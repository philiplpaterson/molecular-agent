from typing import Any

from app.services.rag import rag_service
from app.tools.base import BaseTool


class LiteratureSearchTool(BaseTool):
    """Tool to search scientific literature using RAG."""

    @property
    def name(self) -> str:
        return "search_literature"

    @property
    def description(self) -> str:
        return "Search scientific literature (PubMed abstracts) for information about drugs, mechanisms, clinical studies, and drug discovery topics. Returns relevant papers with citations."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "The search query (e.g., 'mechanism of aspirin', 'EGFR inhibitors clinical trials')",
                },
                "num_results": {
                    "type": "integer",
                    "description": "Number of results to return (default: 5, max: 10)",
                    "minimum": 1,
                    "maximum": 10,
                },
            },
            "required": ["query"],
        }

    def execute(self, query: str, num_results: int = 5, **kwargs: Any) -> dict:
        """
        Search scientific literature.

        Args:
            query: The search query
            num_results: Number of results to return

        Returns:
            Dictionary with search results and citations
        """
        num_results = max(1, min(10, num_results))

        try:
            results = rag_service.search_literature(query, n_results=num_results)

            if not results:
                return {
                    "query": query,
                    "num_results": 0,
                    "results": [],
                    "message": "No relevant literature found. The database may need to be populated with data.",
                }

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
                    "snippet": r.get("content", "")[:500] + "..." if len(r.get("content", "")) > 500 else r.get("content", ""),
                })

            return {
                "query": query,
                "num_results": len(formatted_results),
                "results": formatted_results,
            }

        except Exception as e:
            return {
                "query": query,
                "error": f"Error searching literature: {str(e)}",
            }
