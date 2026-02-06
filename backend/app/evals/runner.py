"""Evaluation runner utilities for MolecularAgent."""

import json
from datetime import datetime
from pathlib import Path
from typing import Any
from uuid import uuid4

from pydantic_evals import Case, Dataset

from app.evals.evaluators import (
    CombinedEvaluator,
    LLMJudgeEvaluator,
    ResponseContainsEvaluator,
    ResponseNotContainsEvaluator,
    ToolSelectionEvaluator,
)
from app.evals.task import (
    MolecularAgentExpectedOutput,
    MolecularAgentInput,
    MolecularAgentOutput,
    molecular_agent_task,
)


def load_dataset(path: str | Path) -> Dataset[MolecularAgentInput, MolecularAgentOutput]:
    """Load evaluation cases from a JSON file.

    Expected JSON format:
    {
        "name": "dataset-name",
        "description": "Dataset description",
        "cases": [
            {
                "id": "case-001",
                "query": "Is CCO a valid SMILES?",
                "expected_tools": ["validate_smiles"],
                "expected_contains": ["valid", "ethanol"],
                "expected_not_contains": ["invalid"],
                "tags": ["tool-selection", "smiles"]
            }
        ]
    }

    Args:
        path: Path to JSON dataset file

    Returns:
        pydantic_evals Dataset ready for evaluation
    """
    path = Path(path)

    with open(path) as f:
        data = json.load(f)

    cases = []
    for case_data in data.get("cases", []):
        case_input = MolecularAgentInput(query=case_data["query"])

        # Create expected output for evaluators to use
        expected = MolecularAgentOutput(
            response="",  # Not checking exact response
            tools_called=[],
            expected_tools=case_data.get("expected_tools", []),
            expected_contains=case_data.get("expected_contains", []),
            expected_not_contains=case_data.get("expected_not_contains", []),
        )

        case = Case(
            name=case_data.get("id", f"case-{len(cases)}"),
            inputs=case_input,
            expected_output=expected,
            metadata={
                "tags": case_data.get("tags", []),
                "description": case_data.get("description", ""),
            },
        )
        cases.append(case)

    return Dataset(
        cases=cases,
        evaluators=[
            ToolSelectionEvaluator(),
            ResponseContainsEvaluator(),
            ResponseNotContainsEvaluator(),
        ],
    )


def load_datasets_from_dir(
    directory: str | Path,
) -> list[tuple[str, Dataset[MolecularAgentInput, MolecularAgentOutput]]]:
    """Load all dataset JSON files from a directory.

    Args:
        directory: Path to directory containing JSON dataset files

    Returns:
        List of (name, Dataset) tuples
    """
    directory = Path(directory)
    datasets = []

    for json_file in directory.glob("*.json"):
        try:
            dataset = load_dataset(json_file)
            datasets.append((json_file.stem, dataset))
        except Exception as e:
            print(f"Warning: Failed to load {json_file}: {e}")

    return datasets


def run_evaluation(
    dataset_path: str | Path,
    output_dir: str | Path | None = None,
    use_llm_judge: bool = False,
    verbose: bool = False,
) -> dict[str, Any]:
    """Run evaluation on a dataset and optionally save results.

    Args:
        dataset_path: Path to dataset JSON file
        output_dir: Optional directory to save results
        use_llm_judge: Whether to use LLM-as-judge (slower, costs money)
        verbose: Print detailed output for each case

    Returns:
        Dictionary with evaluation results and statistics
    """
    dataset = load_dataset(dataset_path)

    # Add LLM judge if requested
    if use_llm_judge:
        dataset = Dataset(
            cases=dataset.cases,
            evaluators=[
                ToolSelectionEvaluator(),
                ResponseContainsEvaluator(),
                ResponseNotContainsEvaluator(),
                LLMJudgeEvaluator(),
            ],
        )

    # Run evaluation
    report = dataset.evaluate_sync(molecular_agent_task)

    if verbose:
        report.print(include_input=True, include_output=True)

    # Build results summary
    results = {
        "run_id": str(uuid4()),
        "timestamp": datetime.utcnow().isoformat(),
        "dataset_path": str(dataset_path),
        "total_cases": len(dataset.cases),
        "use_llm_judge": use_llm_judge,
        "report": _extract_report_data(report),
    }

    # Save results if output directory specified
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / f"eval-{results['run_id'][:8]}-{datetime.utcnow().strftime('%Y%m%d-%H%M%S')}.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, default=str)

        results["output_file"] = str(output_file)

    return results


def _extract_report_data(report) -> dict[str, Any]:
    """Extract serializable data from pydantic-evals report."""
    # The report structure depends on pydantic-evals version
    # We extract what we can in a version-agnostic way
    try:
        return {
            "summary": str(report),
            "cases": [
                {
                    "name": getattr(case, "name", "unknown"),
                    "passed": getattr(case, "passed", None),
                    "scores": getattr(case, "scores", {}),
                }
                for case in getattr(report, "cases", [])
            ],
        }
    except Exception:
        return {"summary": str(report)}


def compare_reports(
    baseline_path: str | Path,
    current_path: str | Path,
) -> dict[str, Any]:
    """Compare two evaluation reports.

    Args:
        baseline_path: Path to baseline report JSON
        current_path: Path to current report JSON

    Returns:
        Comparison results with regressions and improvements
    """
    with open(baseline_path) as f:
        baseline = json.load(f)
    with open(current_path) as f:
        current = json.load(f)

    # Simple comparison - extend as needed
    comparison = {
        "baseline_run_id": baseline.get("run_id"),
        "current_run_id": current.get("run_id"),
        "baseline_total": baseline.get("total_cases", 0),
        "current_total": current.get("total_cases", 0),
    }

    return comparison
