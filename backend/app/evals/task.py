"""Task definition for MolecularAgent evaluation.

This module defines the input/output types and task function
that wraps the molecular agent for evaluation.
"""

import time
from typing import Any

from pydantic import BaseModel, Field


class MolecularAgentInput(BaseModel):
    """Input for a MolecularAgent evaluation case."""

    query: str = Field(..., description="The user query to send to the agent")


class MolecularAgentExpectedOutput(BaseModel):
    """Expected output for validation."""

    expected_tools: list[str] = Field(
        default_factory=list,
        description="List of tool names expected to be called",
    )
    expected_contains: list[str] = Field(
        default_factory=list,
        description="Keywords/phrases expected in the response",
    )
    expected_not_contains: list[str] = Field(
        default_factory=list,
        description="Keywords/phrases that should NOT appear",
    )


class MolecularAgentOutput(BaseModel):
    """Output from a MolecularAgent evaluation run."""

    response: str = Field(default="", description="The agent's response text")
    tools_called: list[str] = Field(
        default_factory=list,
        description="Names of tools that were called",
    )
    tool_call_details: list[dict[str, Any]] = Field(
        default_factory=list,
        description="Full tool call information including arguments and results",
    )
    latency_ms: int = Field(default=0, description="Response latency in milliseconds")
    error: str | None = Field(default=None, description="Error message if failed")

    # For comparison with expected output
    expected_tools: list[str] = Field(default_factory=list)
    expected_contains: list[str] = Field(default_factory=list)
    expected_not_contains: list[str] = Field(default_factory=list)


async def molecular_agent_task(inputs: MolecularAgentInput) -> MolecularAgentOutput:
    """Task function that runs the molecular agent for evaluation.

    This is the function that pydantic-evals will call for each test case.

    Args:
        inputs: The evaluation case input containing the query

    Returns:
        MolecularAgentOutput with response and tool call information
    """
    from sqlmodel import Session

    from app.core.db import engine
    from app.services.agent import MolecularAgent

    start_time = time.perf_counter()

    try:
        with Session(engine) as session:
            agent = MolecularAgent(session)
            result = await agent.process_message(inputs.query)

            latency_ms = int((time.perf_counter() - start_time) * 1000)

            # Extract tool names from tool_calls
            tools_called = [tc["name"] for tc in result.get("tool_calls", [])]

            return MolecularAgentOutput(
                response=result.get("response", ""),
                tools_called=tools_called,
                tool_call_details=result.get("tool_calls", []),
                latency_ms=latency_ms,
            )

    except Exception as e:
        latency_ms = int((time.perf_counter() - start_time) * 1000)
        return MolecularAgentOutput(
            response="",
            tools_called=[],
            tool_call_details=[],
            latency_ms=latency_ms,
            error=str(e),
        )


def molecular_agent_task_sync(inputs: MolecularAgentInput) -> MolecularAgentOutput:
    """Synchronous wrapper for molecular_agent_task."""
    import asyncio

    return asyncio.run(molecular_agent_task(inputs))
