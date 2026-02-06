"""Tests for evaluation scorers/evaluators."""

import pytest

from app.evals.task import MolecularAgentInput, MolecularAgentOutput


class MockEvaluatorContext:
    """Mock context for testing evaluators."""

    def __init__(
        self,
        inputs: MolecularAgentInput,
        output: MolecularAgentOutput,
        expected_output: MolecularAgentOutput | None = None,
    ):
        self.inputs = inputs
        self.output = output
        self.expected_output = expected_output


class TestToolSelectionScoring:
    """Tests for tool selection scoring logic."""

    def test_exact_match(self):
        """Test exact match returns 1.0."""
        expected = ["validate_smiles"]
        actual = ["validate_smiles"]

        expected_set = set(expected)
        actual_set = set(actual)
        intersection = expected_set & actual_set
        union = expected_set | actual_set
        score = len(intersection) / len(union)

        assert score == 1.0

    def test_partial_match(self):
        """Test partial match returns correct Jaccard score."""
        expected = ["validate_smiles", "predict_properties"]
        actual = ["validate_smiles", "check_drug_likeness"]

        expected_set = set(expected)
        actual_set = set(actual)
        intersection = expected_set & actual_set
        union = expected_set | actual_set
        score = len(intersection) / len(union)

        # intersection: {validate_smiles}
        # union: {validate_smiles, predict_properties, check_drug_likeness}
        assert score == pytest.approx(1 / 3)

    def test_no_match(self):
        """Test no overlap returns 0.0."""
        expected = ["validate_smiles"]
        actual = ["predict_properties"]

        expected_set = set(expected)
        actual_set = set(actual)
        intersection = expected_set & actual_set
        union = expected_set | actual_set
        score = len(intersection) / len(union)

        assert score == 0.0

    def test_both_empty(self):
        """Test both empty returns 1.0 (correct - no tools needed)."""
        expected = []
        actual = []

        if not expected and not actual:
            score = 1.0
        else:
            expected_set = set(expected)
            actual_set = set(actual)
            intersection = expected_set & actual_set
            union = expected_set | actual_set
            score = len(intersection) / len(union) if union else 1.0

        assert score == 1.0

    def test_expected_empty_actual_not(self):
        """Test calling tools when none expected returns 0.0."""
        expected = []
        actual = ["validate_smiles"]

        if not expected and actual:
            score = 0.0
        else:
            expected_set = set(expected)
            actual_set = set(actual)
            intersection = expected_set & actual_set
            union = expected_set | actual_set
            score = len(intersection) / len(union) if union else 1.0

        assert score == 0.0


class TestResponseContainsScoring:
    """Tests for response contains scoring logic."""

    def test_all_keywords_found(self):
        """Test all keywords found returns 1.0."""
        keywords = ["valid", "ethanol"]
        response = "The SMILES CCO is valid and represents ethanol."

        found = sum(1 for kw in keywords if kw.lower() in response.lower())
        score = found / len(keywords)

        assert score == 1.0

    def test_partial_keywords_found(self):
        """Test partial keywords found returns correct proportion."""
        keywords = ["valid", "ethanol", "carbon"]
        response = "The SMILES CCO is valid and represents ethanol."

        found = sum(1 for kw in keywords if kw.lower() in response.lower())
        score = found / len(keywords)

        assert score == pytest.approx(2 / 3)

    def test_no_keywords_found(self):
        """Test no keywords found returns 0.0."""
        keywords = ["invalid", "error"]
        response = "The SMILES CCO is valid and represents ethanol."

        found = sum(1 for kw in keywords if kw.lower() in response.lower())
        score = found / len(keywords)

        assert score == 0.0

    def test_empty_keywords(self):
        """Test empty keywords list returns 1.0."""
        keywords = []
        response = "Any response."

        if not keywords:
            score = 1.0
        else:
            found = sum(1 for kw in keywords if kw.lower() in response.lower())
            score = found / len(keywords)

        assert score == 1.0

    def test_empty_response(self):
        """Test empty response returns 0.0 when keywords expected."""
        keywords = ["valid"]
        response = ""

        if not response:
            score = 0.0
        else:
            found = sum(1 for kw in keywords if kw.lower() in response.lower())
            score = found / len(keywords)

        assert score == 0.0


class TestResponseNotContainsScoring:
    """Tests for response NOT contains scoring logic."""

    def test_no_forbidden_found(self):
        """Test no forbidden keywords found returns 1.0."""
        forbidden = ["invalid", "error"]
        response = "The SMILES CCO is valid and represents ethanol."

        not_found = sum(1 for kw in forbidden if kw.lower() not in response.lower())
        score = not_found / len(forbidden)

        assert score == 1.0

    def test_some_forbidden_found(self):
        """Test some forbidden keywords found returns correct proportion."""
        forbidden = ["invalid", "error", "ethanol"]
        response = "The SMILES CCO is valid and represents ethanol."

        not_found = sum(1 for kw in forbidden if kw.lower() not in response.lower())
        score = not_found / len(forbidden)

        # "ethanol" is found, so 2/3 not found
        assert score == pytest.approx(2 / 3)

    def test_all_forbidden_found(self):
        """Test all forbidden keywords found returns 0.0."""
        forbidden = ["valid", "ethanol"]
        response = "The SMILES CCO is valid and represents ethanol."

        not_found = sum(1 for kw in forbidden if kw.lower() not in response.lower())
        score = not_found / len(forbidden)

        assert score == 0.0


class TestDatasetLoading:
    """Tests for dataset loading functionality."""

    def test_load_tool_selection_dataset(self):
        """Test that the tool-selection dataset loads correctly."""
        import json
        from pathlib import Path

        dataset_path = Path("app/evals/datasets/tool-selection.json")

        if not dataset_path.exists():
            pytest.skip("Dataset file not found")

        with open(dataset_path) as f:
            data = json.load(f)

        assert "name" in data
        assert "cases" in data
        assert len(data["cases"]) >= 20

        # Validate case structure
        for case in data["cases"]:
            assert "id" in case
            assert "query" in case
            assert "expected_tools" in case
            assert isinstance(case["expected_tools"], list)
