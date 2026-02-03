import uuid
from datetime import datetime

from sqlmodel import Field, Relationship, SQLModel


class ConversationBase(SQLModel):
    pass


class Conversation(ConversationBase, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)

    messages: list["Message"] = Relationship(back_populates="conversation")


class ConversationCreate(ConversationBase):
    pass


class ConversationPublic(ConversationBase):
    id: uuid.UUID
    created_at: datetime
    updated_at: datetime


class MessageBase(SQLModel):
    role: str  # "user", "assistant", "tool"
    content: str
    tool_calls: str | None = None  # JSON string of tool calls
    tool_call_id: str | None = None  # For tool responses


class Message(MessageBase, table=True):
    id: uuid.UUID = Field(default_factory=uuid.uuid4, primary_key=True)
    conversation_id: uuid.UUID = Field(foreign_key="conversation.id")
    created_at: datetime = Field(default_factory=datetime.utcnow)

    conversation: Conversation = Relationship(back_populates="messages")


class MessageCreate(MessageBase):
    conversation_id: uuid.UUID


class MessagePublic(MessageBase):
    id: uuid.UUID
    conversation_id: uuid.UUID
    created_at: datetime
