from typing import Any

from app.services.rag import rag_service
from app.tools.base import BaseTool


class CompoundSearchTool(BaseTool):
    """Tool to search for compounds using RAG."""

    @property
    def name(self) -> str:
        return "search_compounds"

    @property
    def description(self) -> str:
        return "Search for compounds from ChEMBL database including approved drugs and clinical candidates. Can filter by clinical phase. Returns compound names, SMILES, and therapeutic areas."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "The search query (e.g., 'kinase inhibitor', 'anticancer', 'diabetes treatment')",
                },
                "num_results": {
                    "type": "integer",
                    "description": "Number of results to return (default: 5, max: 10)",
                    "minimum": 1,
                    "maximum": 10,
                },
                "min_phase": {
                    "type": "integer",
                    "description": "Minimum clinical phase (1-4, where 4 is approved)",
                    "minimum": 1,
                    "maximum": 4,
                },
            },
            "required": ["query"],
        }

    def execute(
        self,
        query: str,
        num_results: int = 5,
        min_phase: int | None = None,
        **kwargs: Any,
    ) -> dict:
        """
        Search for compounds.

        Args:
            query: The search query
            num_results: Number of results to return
            min_phase: Minimum clinical phase filter

        Returns:
            Dictionary with compound information
        """
        num_results = max(1, min(10, num_results))

        try:
            results = rag_service.search_compounds(
                query,
                n_results=num_results,
                max_phase=min_phase,
            )

            if not results:
                return {
                    "query": query,
                    "num_results": 0,
                    "compounds": [],
                    "message": "No compounds found. The database may need to be populated with data.",
                }

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
                "query": query,
                "num_results": len(compounds),
                "filter_min_phase": min_phase,
                "compounds": compounds,
            }

        except Exception as e:
            return {
                "query": query,
                "error": f"Error searching compounds: {str(e)}",
            }
