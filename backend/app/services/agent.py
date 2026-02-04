import json
import uuid
from dataclasses import dataclass, field

from pydantic_ai import Agent, RunContext
from pydantic_ai.messages import ModelRequest, ModelResponse, TextPart, UserPromptPart
from sqlmodel import Session

from app import crud
from app.core.config import settings
from app.models import ConversationCreate, MessageCreate
from app.tools import (
    CompoundSearchTool,
    DrugLikenessTool,
    LiteratureSearchTool,
    MoleculeGeneratorTool,
    PropertyPredictorTool,
    SimilaritySearchTool,
    SMILESValidatorTool,
)


@dataclass
class AgentDeps:
    """Dependencies passed to the agent during execution."""

    session: Session
    conversation_id: uuid.UUID
    tool_calls: list[dict] = field(default_factory=list)


def get_model_name() -> str:
    """Get the model name based on configuration."""
    if settings.LLM_PROVIDER == "anthropic":
        return f"anthropic:{settings.ANTHROPIC_MODEL}"
    elif settings.LLM_PROVIDER == "openai":
        return f"openai:{settings.OPENAI_MODEL}"
    else:
        raise ValueError(f'Error: LLM_PROVIDER must be either "openai" or "anthropic".')


SYSTEM_PROMPT = """You are MolecularAgent, an AI assistant specialized in drug discovery and molecular analysis.

You have access to the following tools:

MOLECULAR ANALYSIS TOOLS:
1. validate_smiles - Validate SMILES strings and convert to canonical form
2. predict_properties - Predict molecular properties (MW, LogP, TPSA, etc.)
3. check_drug_likeness - Check Lipinski's Rule of Five for drug-likeness assessment
4. similarity_search - Calculate molecular similarity between compounds
5. generate_molecules - Generate random drug-like molecules

KNOWLEDGE RETRIEVAL TOOLS (RAG):
6. search_literature - Search scientific literature (PubMed abstracts) for papers and research
7. search_compounds - Search ChEMBL for compounds and clinical candidates

When users ask about molecules, use the appropriate tools to provide accurate information.
Always validate SMILES strings before using them in other tools.
For questions about scientific research, use the search_literature tool.
For questions about compounds and clinical candidates, use the search_compounds tool.
Cite sources when providing information from the literature.
Explain results in a clear, accessible way while maintaining scientific accuracy.
"""

# Create the PydanticAI agent
molecular_agent = Agent(
    get_model_name(),
    deps_type=AgentDeps,
    system_prompt=SYSTEM_PROMPT,
)


@molecular_agent.tool
def validate_smiles(ctx: RunContext[AgentDeps], smiles: str) -> dict:
    """
    Validate a SMILES string and return its canonical form.

    Args:
        smiles: The SMILES string to validate (e.g., 'CCO' for ethanol)
    """
    tool = SMILESValidatorTool()
    result = tool.execute(smiles=smiles)
    ctx.deps.tool_calls.append({"name": "validate_smiles", "arguments": {"smiles": smiles}, "result": result})
    return result


@molecular_agent.tool
def predict_properties(ctx: RunContext[AgentDeps], smiles: str) -> dict:
    """
    Predict molecular properties from a SMILES string.

    Args:
        smiles: The SMILES string of the molecule (must be valid)
    """
    tool = PropertyPredictorTool()
    result = tool.execute(smiles=smiles)
    ctx.deps.tool_calls.append({"name": "predict_properties", "arguments": {"smiles": smiles}, "result": result})
    return result


@molecular_agent.tool
def check_drug_likeness(ctx: RunContext[AgentDeps], smiles: str) -> dict:
    """
    Check if a molecule passes Lipinski's Rule of Five for drug-likeness.

    Args:
        smiles: The SMILES string of the molecule to check
    """
    tool = DrugLikenessTool()
    result = tool.execute(smiles=smiles)
    ctx.deps.tool_calls.append({"name": "check_drug_likeness", "arguments": {"smiles": smiles}, "result": result})
    return result


@molecular_agent.tool
def similarity_search(ctx: RunContext[AgentDeps], query_smiles: str, target_smiles: list[str]) -> dict:
    """
    Calculate Tanimoto similarity between a query molecule and target molecules.

    Args:
        query_smiles: The SMILES string of the query molecule
        target_smiles: List of SMILES strings to compare against
    """
    tool = SimilaritySearchTool()
    result = tool.execute(query_smiles=query_smiles, target_smiles=target_smiles)
    ctx.deps.tool_calls.append({
        "name": "similarity_search",
        "arguments": {"query_smiles": query_smiles, "target_smiles": target_smiles},
        "result": result,
    })
    return result


@molecular_agent.tool
def generate_molecules(ctx: RunContext[AgentDeps], num_molecules: int, seed_smiles: str | None = None) -> dict:
    """
    Generate random drug-like molecules.

    Args:
        num_molecules: Number of molecules to generate (1-20)
        seed_smiles: Optional seed SMILES to generate similar molecules
    """
    tool = MoleculeGeneratorTool()
    result = tool.execute(num_molecules=num_molecules, seed_smiles=seed_smiles)
    ctx.deps.tool_calls.append({
        "name": "generate_molecules",
        "arguments": {"num_molecules": num_molecules, "seed_smiles": seed_smiles},
        "result": result,
    })
    return result


# --- RAG Tools ---


@molecular_agent.tool
def search_literature(ctx: RunContext[AgentDeps], query: str, num_results: int = 5) -> dict:
    """
    Search scientific literature (PubMed abstracts) for relevant papers.

    Args:
        query: The search query (e.g., 'mechanism of aspirin', 'EGFR inhibitors')
        num_results: Number of results to return (1-10, default 5)
    """
    tool = LiteratureSearchTool()
    result = tool.execute(query=query, num_results=num_results)
    ctx.deps.tool_calls.append({
        "name": "search_literature",
        "arguments": {"query": query, "num_results": num_results},
        "result": result,
    })
    return result


@molecular_agent.tool
def search_compounds(ctx: RunContext[AgentDeps], query: str, num_results: int = 5, min_phase: int | None = None) -> dict:
    """
    Search ChEMBL for compounds including approved drugs and clinical candidates.

    Args:
        query: The search query (e.g., 'kinase inhibitor', 'anticancer')
        num_results: Number of results to return (1-10, default 5)
        min_phase: Minimum clinical phase filter (1-4, where 4 is approved)
    """
    tool = CompoundSearchTool()
    result = tool.execute(query=query, num_results=num_results, min_phase=min_phase)
    ctx.deps.tool_calls.append({
        "name": "search_compounds",
        "arguments": {"query": query, "num_results": num_results, "min_phase": min_phase},
        "result": result,
    })
    return result


class MolecularAgent:
    """High-level wrapper for the molecular agent."""

    def __init__(self, session: Session):
        self.session = session

    async def process_message(
        self,
        user_message: str,
        conversation_id: uuid.UUID | None = None,
    ) -> dict:
        """
        Process a user message and return the agent's response.

        Args:
            user_message: The user's input message
            conversation_id: Optional existing conversation ID

        Returns:
            Dictionary with response, tool_calls, and conversation_id
        """
        # Get or create conversation and subsequent messages
        conversation = crud.get_or_create_conversation(self.session, conversation_id)

        history_messages = crud.get_conversation_messages(self.session, conversation.id)

        message_history: list[ModelRequest | ModelResponse] = []
        for msg in history_messages:
            if msg.role == "user":
                message_history.append(ModelRequest(parts=[UserPromptPart(content=msg.content)]))
            elif msg.role == "assistant":
                message_history.append(ModelResponse(parts=[TextPart(content=msg.content)]))

        # Save user message
        crud.create_message(
            self.session,
            MessageCreate(
                conversation_id=conversation.id,
                role="user",
                content=user_message,
            ),
        )

        # Run the agent
        deps = AgentDeps(session=self.session, conversation_id=conversation.id, tool_calls=[])
        result = await molecular_agent.run(user_message, deps=deps, message_history=message_history)

        # Save assistant response
        crud.create_message(
            self.session,
            MessageCreate(
                conversation_id=conversation.id,
                role="assistant",
                content=result.output,
                tool_calls=json.dumps(deps.tool_calls) if deps.tool_calls else None,
            ),
        )

        return {
            "response": result.output,
            "tool_calls": deps.tool_calls,
            "conversation_id": conversation.id,
        }
