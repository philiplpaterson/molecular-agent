from fastapi import APIRouter
from pydantic import BaseModel

from app.services.rag import rag_service

router = APIRouter()


class SearchQuery(BaseModel):
    query: str
    num_results: int = 5


class LiteratureResult(BaseModel):
    title: str
    authors: str
    journal: str
    year: str
    pmid: str
    relevance_score: float
    snippet: str


class LiteratureSearchResponse(BaseModel):
    query: str
    num_results: int
    results: list[LiteratureResult]


class CompoundResult(BaseModel):
    name: str
    chembl_id: str
    smiles: str
    max_phase: str
    molecule_type: str
    therapeutic_areas: str
    relevance_score: float


class CompoundSearchResponse(BaseModel):
    query: str
    num_results: int
    compounds: list[CompoundResult]


class CollectionStats(BaseModel):
    literature: int
    compounds: int


@router.post("/literature/search")
def search_literature(query: SearchQuery) -> dict:
    """Search scientific literature."""
    results = rag_service.search_literature(query.query, n_results=query.num_results)

    formatted_results = []
    for r in results:
        metadata = r.get("metadata", {})
        formatted_results.append({
            "title": metadata.get("title", ""),
            "authors": metadata.get("authors", ""),
            "journal": metadata.get("journal", ""),
            "year": metadata.get("year", ""),
            "pmid": metadata.get("pmid", ""),
            "relevance_score": r.get("similarity", 0),
            "snippet": r.get("content", "")[:500],
        })

    return {
        "query": query.query,
        "num_results": len(formatted_results),
        "results": formatted_results,
    }


class CompoundSearchQuery(BaseModel):
    query: str
    num_results: int = 5
    min_phase: int | None = None


@router.post("/compounds/search")
def search_compounds(query: CompoundSearchQuery) -> dict:
    """Search ChEMBL compounds."""
    results = rag_service.search_compounds(
        query.query,
        n_results=query.num_results,
        max_phase=query.min_phase,
    )

    compounds = []
    for r in results:
        metadata = r.get("metadata", {})
        compounds.append({
            "name": metadata.get("name", ""),
            "chembl_id": metadata.get("chembl_id", ""),
            "smiles": metadata.get("smiles", ""),
            "max_phase": metadata.get("max_phase", ""),
            "molecule_type": metadata.get("molecule_type", ""),
            "therapeutic_areas": metadata.get("therapeutic_areas", ""),
            "relevance_score": r.get("similarity", 0),
        })

    return {
        "query": query.query,
        "num_results": len(compounds),
        "compounds": compounds,
    }


@router.get("/stats")
def get_stats() -> CollectionStats:
    """Get collection statistics."""
    stats = rag_service.get_collection_stats()
    return CollectionStats(
        literature=stats.get("literature", 0),
        compounds=stats.get("compounds", 0),
    )
