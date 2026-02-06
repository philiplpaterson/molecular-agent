"""API routes for running and viewing evaluations."""

import json
from datetime import datetime
from pathlib import Path
from typing import Any
from uuid import uuid4

from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel, Field

router = APIRouter()

# Store for tracking background eval runs
_eval_runs: dict[str, dict[str, Any]] = {}


class EvalRunRequest(BaseModel):
    """Request to start an evaluation run."""

    dataset_path: str = Field(
        default="app/evals/datasets/tool-selection.json",
        description="Path to dataset JSON file",
    )
    use_llm_judge: bool = Field(
        default=False,
        description="Use LLM-as-judge for quality scoring",
    )


class EvalRunResponse(BaseModel):
    """Response from starting an evaluation run."""

    run_id: str
    status: str
    message: str


class EvalReportSummary(BaseModel):
    """Summary of an evaluation report."""

    run_id: str
    timestamp: str
    dataset_path: str
    total_cases: int
    status: str


async def _run_evaluation_background(
    run_id: str,
    dataset_path: str,
    use_llm_judge: bool,
) -> None:
    """Background task to run evaluation."""
    from app.evals.runner import run_evaluation

    try:
        _eval_runs[run_id]["status"] = "running"
        _eval_runs[run_id]["started_at"] = datetime.utcnow().isoformat()

        result = run_evaluation(
            dataset_path=dataset_path,
            output_dir="app/evals/reports",
            use_llm_judge=use_llm_judge,
            verbose=False,
        )

        _eval_runs[run_id]["status"] = "completed"
        _eval_runs[run_id]["completed_at"] = datetime.utcnow().isoformat()
        _eval_runs[run_id]["result"] = result

    except Exception as e:
        _eval_runs[run_id]["status"] = "failed"
        _eval_runs[run_id]["error"] = str(e)
        _eval_runs[run_id]["completed_at"] = datetime.utcnow().isoformat()


@router.post("/run", response_model=EvalRunResponse)
async def start_eval_run(
    request: EvalRunRequest,
    background_tasks: BackgroundTasks,
) -> EvalRunResponse:
    """Start an evaluation run in the background.

    The evaluation runs asynchronously. Use GET /evals/run/{run_id}
    to check status and get results.
    """
    run_id = str(uuid4())

    # Validate dataset exists
    dataset_path = Path(request.dataset_path)
    if not dataset_path.exists():
        raise HTTPException(
            status_code=400,
            detail=f"Dataset not found: {request.dataset_path}",
        )

    # Initialize run tracking
    _eval_runs[run_id] = {
        "run_id": run_id,
        "status": "pending",
        "dataset_path": request.dataset_path,
        "use_llm_judge": request.use_llm_judge,
        "created_at": datetime.utcnow().isoformat(),
    }

    # Start background task
    background_tasks.add_task(
        _run_evaluation_background,
        run_id,
        request.dataset_path,
        request.use_llm_judge,
    )

    return EvalRunResponse(
        run_id=run_id,
        status="pending",
        message="Evaluation started. Use GET /evals/run/{run_id} to check status.",
    )


@router.get("/run/{run_id}")
async def get_eval_run(run_id: str) -> dict[str, Any]:
    """Get status and results of an evaluation run."""
    if run_id not in _eval_runs:
        raise HTTPException(status_code=404, detail=f"Run not found: {run_id}")

    return _eval_runs[run_id]


@router.get("/reports", response_model=list[EvalReportSummary])
async def list_reports() -> list[EvalReportSummary]:
    """List all saved evaluation reports."""
    reports_dir = Path("app/evals/reports")

    if not reports_dir.exists():
        return []

    summaries = []
    for report_file in reports_dir.glob("eval-*.json"):
        try:
            with open(report_file) as f:
                data = json.load(f)

            summaries.append(
                EvalReportSummary(
                    run_id=data.get("run_id", report_file.stem),
                    timestamp=data.get("timestamp", ""),
                    dataset_path=data.get("dataset_path", ""),
                    total_cases=data.get("total_cases", 0),
                    status="completed",
                )
            )
        except Exception:
            continue

    # Sort by timestamp descending
    summaries.sort(key=lambda x: x.timestamp, reverse=True)
    return summaries


@router.get("/reports/{run_id}")
async def get_report(run_id: str) -> dict[str, Any]:
    """Get a specific evaluation report by run ID."""
    reports_dir = Path("app/evals/reports")

    # Look for report file
    for report_file in reports_dir.glob("eval-*.json"):
        try:
            with open(report_file) as f:
                data = json.load(f)
            if data.get("run_id") == run_id or run_id in report_file.stem:
                return data
        except Exception:
            continue

    raise HTTPException(status_code=404, detail=f"Report not found: {run_id}")


@router.get("/datasets")
async def list_datasets() -> list[dict[str, Any]]:
    """List available evaluation datasets."""
    datasets_dir = Path("app/evals/datasets")

    if not datasets_dir.exists():
        return []

    datasets = []
    for dataset_file in datasets_dir.glob("*.json"):
        try:
            with open(dataset_file) as f:
                data = json.load(f)

            datasets.append(
                {
                    "name": data.get("name", dataset_file.stem),
                    "path": str(dataset_file),
                    "description": data.get("description", ""),
                    "num_cases": len(data.get("cases", [])),
                    "version": data.get("version", ""),
                }
            )
        except Exception:
            continue

    return datasets
