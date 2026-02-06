"""Evaluation framework for MolecularAgent using pydantic-evals.

This module provides tools for evaluating agent performance:
- Custom evaluators for tool selection, response quality, and LLM-as-judge
- Dataset loading utilities
- Evaluation runner for batch testing

Usage:
    from app.evals import run_evaluation, load_dataset

    # Load and run evaluation
    results = run_evaluation("app/evals/datasets/tool-selection.json", verbose=True)

    # Or use CLI:
    # python scripts/run-evals.py --dataset app/evals/datasets/tool-selection.json
"""

# Lazy imports to avoid circular dependencies and improve startup time
__all__ = [
    "ToolSelectionEvaluator",
    "ResponseContainsEvaluator",
    "ResponseNotContainsEvaluator",
    "LLMJudgeEvaluator",
    "CombinedEvaluator",
    "run_evaluation",
    "load_dataset",
    "load_datasets_from_dir",
    "molecular_agent_task",
    "MolecularAgentInput",
    "MolecularAgentOutput",
    "MolecularAgentExpectedOutput",
]


def __getattr__(name: str):
    """Lazy import to avoid circular dependencies."""
    if name in (
        "ToolSelectionEvaluator",
        "ResponseContainsEvaluator",
        "ResponseNotContainsEvaluator",
        "LLMJudgeEvaluator",
        "CombinedEvaluator",
    ):
        from app.evals.evaluators import (
            CombinedEvaluator,
            LLMJudgeEvaluator,
            ResponseContainsEvaluator,
            ResponseNotContainsEvaluator,
            ToolSelectionEvaluator,
        )

        return {
            "ToolSelectionEvaluator": ToolSelectionEvaluator,
            "ResponseContainsEvaluator": ResponseContainsEvaluator,
            "ResponseNotContainsEvaluator": ResponseNotContainsEvaluator,
            "LLMJudgeEvaluator": LLMJudgeEvaluator,
            "CombinedEvaluator": CombinedEvaluator,
        }[name]

    if name in ("run_evaluation", "load_dataset", "load_datasets_from_dir"):
        from app.evals.runner import load_dataset, load_datasets_from_dir, run_evaluation

        return {
            "run_evaluation": run_evaluation,
            "load_dataset": load_dataset,
            "load_datasets_from_dir": load_datasets_from_dir,
        }[name]

    if name in (
        "molecular_agent_task",
        "MolecularAgentInput",
        "MolecularAgentOutput",
        "MolecularAgentExpectedOutput",
    ):
        from app.evals.task import (
            MolecularAgentExpectedOutput,
            MolecularAgentInput,
            MolecularAgentOutput,
            molecular_agent_task,
        )

        return {
            "molecular_agent_task": molecular_agent_task,
            "MolecularAgentInput": MolecularAgentInput,
            "MolecularAgentOutput": MolecularAgentOutput,
            "MolecularAgentExpectedOutput": MolecularAgentExpectedOutput,
        }[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
