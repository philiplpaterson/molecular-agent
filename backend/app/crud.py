import uuid
from datetime import datetime

from sqlmodel import Session, select

from app.models import Conversation, ConversationCreate, Message, MessageCreate


def create_conversation(session: Session, conversation_in: ConversationCreate) -> Conversation:
    db_conversation = Conversation.model_validate(conversation_in)
    session.add(db_conversation)
    session.commit()
    session.refresh(db_conversation)
    return db_conversation


def get_conversation(session: Session, conversation_id: uuid.UUID) -> Conversation | None:
    return session.get(Conversation, conversation_id)


def get_conversation_messages(session: Session, conversation_id: uuid.UUID) -> list[Message]:
    statement = select(Message).where(Message.conversation_id == conversation_id).order_by(Message.created_at)
    return list(session.exec(statement).all())


def create_message(session: Session, message_in: MessageCreate) -> Message:
    db_message = Message.model_validate(message_in)
    session.add(db_message)

    # Update conversation's updated_at
    conversation = session.get(Conversation, message_in.conversation_id)
    if conversation:
        conversation.updated_at = datetime.utcnow()
        session.add(conversation)

    session.commit()
    session.refresh(db_message)
    return db_message


def get_or_create_conversation(session: Session, conversation_id: uuid.UUID | None) -> Conversation:
    if conversation_id:
        conversation = get_conversation(session, conversation_id)
        if conversation:
            return conversation
    return create_conversation(session, ConversationCreate())
