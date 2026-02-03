import uuid

from fastapi import APIRouter
from pydantic import BaseModel

from app.api.deps import SessionDep

router = APIRouter()


class ChatRequest(BaseModel):
    message: str
    conversation_id: uuid.UUID | None = None


class ToolCall(BaseModel):
    name: str
    arguments: dict
    result: dict


class ChatResponse(BaseModel):
    response: str
    tool_calls: list[ToolCall]
    conversation_id: uuid.UUID


@router.post("/chat", response_model=ChatResponse)
async def chat(request: ChatRequest, session: SessionDep):
    from app.services.agent import MolecularAgent

    agent = MolecularAgent(session=session)
    result = await agent.process_message(
        user_message=request.message,
        conversation_id=request.conversation_id,
    )
    return ChatResponse(
        response=result["response"],
        tool_calls=[ToolCall(**tc) for tc in result["tool_calls"]],
        conversation_id=result["conversation_id"],
    )
