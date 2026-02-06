"""Custom evaluators for MolecularAgent using pydantic-evals framework."""

from dataclasses import dataclass
from typing import Any

from pydantic_evals.evaluators import Evaluator, EvaluatorContext

from app.evals.task import MolecularAgentInput, MolecularAgentOutput


@dataclass
class ToolSelectionEvaluator(Evaluator[MolecularAgentInput, MolecularAgentOutput]):
    """Evaluator that checks if the correct tools were called.

    Uses Jaccard similarity (intersection over union) for partial credit.
    """

    async def evaluate(
        self, ctx: EvaluatorContext[MolecularAgentInput, MolecularAgentOutput]
    ) -> float:
        """Evaluate tool selection accuracy.

        Returns:
            Score from 0.0 to 1.0 based on tool selection overlap
        """
        expected = ctx.expected_output
        actual = ctx.output

        if expected is None:
            return 1.0

        expected_tools = set(expected.expected_tools) if expected.expected_tools else set()
        actual_tools = set(actual.tools_called) if actual.tools_called else set()

        # Both empty - correct (no tools needed)
        if not expected_tools and not actual_tools:
            return 1.0

        # Called tools when none expected
        if not expected_tools and actual_tools:
            return 0.0

        # Didn't call tools when some expected
        if expected_tools and not actual_tools:
            return 0.0

        # Calculate Jaccard similarity
        intersection = expected_tools & actual_tools
        union = expected_tools | actual_tools

        return len(intersection) / len(union) if union else 1.0


@dataclass
class ResponseContainsEvaluator(Evaluator[MolecularAgentInput, MolecularAgentOutput]):
    """Evaluator that checks if response contains expected keywords."""

    case_sensitive: bool = False

    async def evaluate(
        self, ctx: EvaluatorContext[MolecularAgentInput, MolecularAgentOutput]
    ) -> float:
        """Evaluate presence of expected keywords in response.

        Returns:
            Score from 0.0 to 1.0 representing proportion of keywords found
        """
        expected = ctx.expected_output
        actual = ctx.output

        if expected is None or not expected.expected_contains:
            return 1.0

        if not actual.response:
            return 0.0

        response = actual.response if self.case_sensitive else actual.response.lower()
        keywords = expected.expected_contains

        found_count = 0
        for keyword in keywords:
            check_keyword = keyword if self.case_sensitive else keyword.lower()
            if check_keyword in response:
                found_count += 1

        return found_count / len(keywords)


@dataclass
class ResponseNotContainsEvaluator(Evaluator[MolecularAgentInput, MolecularAgentOutput]):
    """Evaluator that checks if response does NOT contain forbidden keywords."""

    case_sensitive: bool = False

    async def evaluate(
        self, ctx: EvaluatorContext[MolecularAgentInput, MolecularAgentOutput]
    ) -> float:
        """Evaluate absence of forbidden keywords in response.

        Returns:
            Score from 0.0 to 1.0 representing proportion of keywords NOT found
        """
        expected = ctx.expected_output
        actual = ctx.output

        if expected is None or not expected.expected_not_contains:
            return 1.0

        if not actual.response:
            return 1.0

        response = actual.response if self.case_sensitive else actual.response.lower()
        forbidden = expected.expected_not_contains

        not_found_count = 0
        for keyword in forbidden:
            check_keyword = keyword if self.case_sensitive else keyword.lower()
            if check_keyword not in response:
                not_found_count += 1

        return not_found_count / len(forbidden)


@dataclass
class LLMJudgeEvaluator(Evaluator[MolecularAgentInput, MolecularAgentOutput]):
    """Evaluator that uses an LLM to judge response quality."""

    model: str = "claude-haiku-4-20250514"

    RUBRIC: str = """You are evaluating an AI agent's response to a drug discovery query.

Rate the response on a scale of 1-5:
1 - Poor: Response is wrong, irrelevant, or harmful
2 - Below Average: Partially correct but missing key information
3 - Average: Correct but could be more helpful
4 - Good: Correct, helpful, and well-structured
5 - Excellent: Comprehensive and demonstrates expert understanding

User Query: {query}
Agent Response: {response}
Tools Called: {tools}

Respond with ONLY a number 1-5."""

    async def evaluate(
        self, ctx: EvaluatorContext[MolecularAgentInput, MolecularAgentOutput]
    ) -> float:
        """Use LLM to judge response quality.

        Returns:
            Normalized score from 0.0 to 1.0 (1-5 scale mapped to 0-1)
        """
        import anthropic
        import os

        api_key = os.environ.get("ANTHROPIC_API_KEY")
        if not api_key:
            return 0.5  # Default to middle score if no API key

        actual = ctx.output
        query = ctx.inputs.query

        prompt = self.RUBRIC.format(
            query=query,
            response=actual.response or "",
            tools=", ".join(actual.tools_called) if actual.tools_called else "None",
        )

        try:
            client = anthropic.AsyncAnthropic(api_key=api_key)
            message = await client.messages.create(
                model=self.model,
                max_tokens=10,
                messages=[{"role": "user", "content": prompt}],
            )

            response_text = message.content[0].text.strip()

            # Extract rating (1-5)
            rating = 3  # Default
            for char in response_text:
                if char.isdigit():
                    rating = int(char)
                    break

            rating = max(1, min(5, rating))

            # Normalize to 0-1
            return (rating - 1) / 4

        except Exception:
            return 0.5  # Default on error


@dataclass
class CombinedEvaluator(Evaluator[MolecularAgentInput, MolecularAgentOutput]):
    """Combines multiple evaluators with configurable weights."""

    tool_weight: float = 0.4
    contains_weight: float = 0.3
    not_contains_weight: float = 0.15
    llm_judge_weight: float = 0.15
    use_llm_judge: bool = False

    async def evaluate(
        self, ctx: EvaluatorContext[MolecularAgentInput, MolecularAgentOutput]
    ) -> float:
        """Calculate weighted combined score from all evaluators."""
        tool_eval = ToolSelectionEvaluator()
        contains_eval = ResponseContainsEvaluator()
        not_contains_eval = ResponseNotContainsEvaluator()

        tool_score = await tool_eval.evaluate(ctx)
        contains_score = await contains_eval.evaluate(ctx)
        not_contains_score = await not_contains_eval.evaluate(ctx)

        total_weight = self.tool_weight + self.contains_weight + self.not_contains_weight
        score = (
            tool_score * self.tool_weight
            + contains_score * self.contains_weight
            + not_contains_score * self.not_contains_weight
        )

        if self.use_llm_judge:
            llm_eval = LLMJudgeEvaluator()
            llm_score = await llm_eval.evaluate(ctx)
            score += llm_score * self.llm_judge_weight
            total_weight += self.llm_judge_weight

        return score / total_weight if total_weight > 0 else 0.0
