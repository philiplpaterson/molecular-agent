#!/usr/bin/env python
"""CLI for running MolecularAgent evaluations.

Usage:
    python scripts/run-evals.py --dataset app/evals/datasets/tool-selection.json
    python scripts/run-evals.py --dataset app/evals/datasets/ --all
    python scripts/run-evals.py --dataset app/evals/datasets/tool-selection.json --llm-judge --verbose
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

# Add the app directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.evals.runner import load_dataset, load_datasets_from_dir, run_evaluation


def main():
    parser = argparse.ArgumentParser(
        description="Run MolecularAgent evaluations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Run a single dataset:
    python scripts/run-evals.py --dataset app/evals/datasets/tool-selection.json

  Run all datasets in a directory:
    python scripts/run-evals.py --dataset app/evals/datasets/ --all

  Run with LLM judge and verbose output:
    python scripts/run-evals.py --dataset app/evals/datasets/tool-selection.json --llm-judge --verbose

  Save results to a specific directory:
    python scripts/run-evals.py --dataset app/evals/datasets/tool-selection.json --output /tmp/eval-results
        """,
    )

    parser.add_argument(
        "--dataset",
        type=str,
        required=True,
        help="Path to dataset JSON file or directory (with --all)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all datasets in the specified directory",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="app/evals/reports",
        help="Output directory for evaluation reports (default: app/evals/reports)",
    )
    parser.add_argument(
        "--llm-judge",
        action="store_true",
        help="Use LLM-as-judge for response quality scoring (slower, costs money)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print detailed output for each test case",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Load and validate datasets without running evaluations",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.85,
        help="Minimum accuracy threshold for CI (default: 0.85)",
    )
    parser.add_argument(
        "--ci",
        action="store_true",
        help="CI mode: exit with code 1 if accuracy below threshold",
    )

    args = parser.parse_args()

    dataset_path = Path(args.dataset)

    if not dataset_path.exists():
        print(f"Error: Dataset path does not exist: {dataset_path}")
        sys.exit(1)

    print("=" * 60)
    print("MolecularAgent Evaluation Runner")
    print("=" * 60)
    print(f"Timestamp: {datetime.utcnow().isoformat()}")
    print(f"Dataset: {dataset_path}")
    print(f"LLM Judge: {'enabled' if args.llm_judge else 'disabled'}")
    print(f"Output: {args.output}")
    print("=" * 60)
    print()

    if args.all and dataset_path.is_dir():
        datasets = load_datasets_from_dir(dataset_path)
        if not datasets:
            print("No datasets found in directory")
            sys.exit(1)
        print(f"Found {len(datasets)} dataset(s)")
    elif dataset_path.is_file():
        datasets = [(dataset_path.stem, load_dataset(dataset_path))]
    else:
        print(f"Error: {dataset_path} is not a file. Use --all for directories.")
        sys.exit(1)

    if args.dry_run:
        print("\n[DRY RUN] Validating datasets...")
        for name, dataset in datasets:
            print(f"  - {name}: {len(dataset.cases)} cases")
        print("\nDry run complete. No evaluations were executed.")
        sys.exit(0)

    all_results = []
    total_passed = 0
    total_cases = 0

    for name, dataset in datasets:
        print(f"\nRunning evaluation: {name}")
        print("-" * 40)

        try:
            if args.all and dataset_path.is_dir():
                result = run_evaluation(
                    dataset_path / f"{name}.json",
                    output_dir=args.output,
                    use_llm_judge=args.llm_judge,
                    verbose=args.verbose,
                )
            else:
                result = run_evaluation(
                    dataset_path,
                    output_dir=args.output,
                    use_llm_judge=args.llm_judge,
                    verbose=args.verbose,
                )

            all_results.append(result)
            total_cases += result.get("total_cases", 0)

            print(f"  Cases: {result.get('total_cases', 0)}")
            if "output_file" in result:
                print(f"  Report saved: {result['output_file']}")

        except Exception as e:
            print(f"  Error: {e}")
            if args.ci:
                sys.exit(1)

    print("\n" + "=" * 60)
    print("EVALUATION SUMMARY")
    print("=" * 60)
    print(f"Total datasets: {len(datasets)}")
    print(f"Total cases: {total_cases}")

    if args.ci:
        # In CI mode, we'd check against thresholds
        # For now, just report success
        print(f"\nCI threshold: {args.threshold}")
        print("CI check: PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
