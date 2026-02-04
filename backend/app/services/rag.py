import hashlib
from typing import Any

import chromadb
from chromadb.config import Settings as ChromaSettings

from app.core.config import settings

# Collection names
LITERATURE_COLLECTION = "literature"
COMPOUNDS_COLLECTION = "compounds"


class RAGService:
    """Service for RAG operations using ChromaDB."""

    def __init__(self):
        self._client = None

    @property
    def client(self) -> chromadb.HttpClient:
        """Lazy initialization of ChromaDB client."""
        if self._client is None:
            self._client = chromadb.HttpClient(
                host=settings.CHROMADB_HOST,
                port=settings.CHROMADB_PORT,
                settings=ChromaSettings(anonymized_telemetry=False),
            )
        return self._client

    def _get_or_create_collection(self, name: str) -> chromadb.Collection:
        """Get or create a collection."""
        return self.client.get_or_create_collection(
            name=name,
            metadata={"hnsw:space": "cosine"},
        )

    def _generate_id(self, text: str, prefix: str = "") -> str:
        """Generate a unique ID for a document."""
        hash_val = hashlib.md5(text.encode()).hexdigest()[:12]
        return f"{prefix}_{hash_val}" if prefix else hash_val

    # --- Literature Operations ---

    def index_literature(
        self,
        documents: list[dict[str, Any]],
    ) -> int:
        """
        Index literature documents (PubMed abstracts, papers).

        Args:
            documents: List of dicts with keys: title, abstract, authors, pmid, year, journal

        Returns:
            Number of documents indexed
        """
        collection = self._get_or_create_collection(LITERATURE_COLLECTION)

        ids = []
        texts = []
        metadatas = []

        seen_ids = set()
        for doc in documents:
            # Combine title and abstract for embedding
            text = f"{doc.get('title', '')} {doc.get('abstract', '')}"
            doc_id = self._generate_id(text, prefix=doc.get("pmid", "lit"))

            # Skip duplicates within batch
            if doc_id in seen_ids:
                continue
            seen_ids.add(doc_id)

            ids.append(doc_id)
            texts.append(text)
            metadatas.append({
                "title": doc.get("title", ""),
                "authors": doc.get("authors", ""),
                "pmid": doc.get("pmid", ""),
                "year": str(doc.get("year", "")),
                "journal": doc.get("journal", ""),
                "source": "pubmed",
            })

        if texts:
            collection.upsert(ids=ids, documents=texts, metadatas=metadatas)

        return len(texts)

    def search_literature(
        self,
        query: str,
        n_results: int = 5,
        year_min: int | None = None,
    ) -> list[dict[str, Any]]:
        """
        Search literature by semantic similarity.

        Args:
            query: Search query
            n_results: Number of results to return
            year_min: Optional minimum year filter

        Returns:
            List of matching documents with scores
        """
        collection = self._get_or_create_collection(LITERATURE_COLLECTION)

        where_filter = None
        if year_min:
            where_filter = {"year": {"$gte": str(year_min)}}

        results = collection.query(
            query_texts=[query],
            n_results=n_results,
            where=where_filter,
            include=["documents", "metadatas", "distances"],
        )

        return self._format_results(results)

    # --- Compound Operations ---

    def index_compounds(
        self,
        compounds: list[dict[str, Any]],
    ) -> int:
        """
        Index compound information (ChEMBL data).

        Args:
            compounds: List of dicts with keys: chembl_id, smiles, name, description,
                       max_phase, molecule_type, therapeutic_areas

        Returns:
            Number of compounds indexed
        """
        collection = self._get_or_create_collection(COMPOUNDS_COLLECTION)

        ids = []
        texts = []
        metadatas = []

        seen_ids = set()
        for compound in compounds:
            # Create searchable text
            text = " ".join([
                compound.get("name", ""),
                compound.get("description", ""),
                compound.get("therapeutic_areas", ""),
            ])

            if not text.strip():
                text = compound.get("smiles", "")

            compound_id = self._generate_id(
                compound.get("smiles", text),
                prefix=compound.get("chembl_id", "cmpd"),
            )

            # Skip duplicates within batch
            if compound_id in seen_ids:
                continue
            seen_ids.add(compound_id)

            ids.append(compound_id)
            texts.append(text)
            metadatas.append({
                "chembl_id": compound.get("chembl_id", ""),
                "smiles": compound.get("smiles", ""),
                "name": compound.get("name", ""),
                "max_phase": int(compound.get("max_phase", 0) or 0),
                "molecule_type": compound.get("molecule_type", ""),
                "therapeutic_areas": compound.get("therapeutic_areas", ""),
                "source": "chembl",
            })

        if texts:
            collection.upsert(ids=ids, documents=texts, metadatas=metadatas)

        return len(texts)

    def search_compounds(
        self,
        query: str,
        n_results: int = 5,
        max_phase: int | None = None,
    ) -> list[dict[str, Any]]:
        """
        Search compounds by semantic similarity.

        Args:
            query: Search query
            n_results: Number of results to return
            max_phase: Optional filter by max clinical phase (1-4)

        Returns:
            List of matching compounds with scores
        """
        collection = self._get_or_create_collection(COMPOUNDS_COLLECTION)

        where_filter = None
        if max_phase:
            where_filter = {"max_phase": {"$gte": max_phase}}

        results = collection.query(
            query_texts=[query],
            n_results=n_results,
            where=where_filter,
            include=["documents", "metadatas", "distances"],
        )

        return self._format_results(results)

    # --- Utility Methods ---

    def _format_results(self, results: dict) -> list[dict[str, Any]]:
        """Format ChromaDB results into a cleaner structure."""
        formatted = []

        if not results or not results.get("ids") or not results["ids"][0]:
            return formatted

        ids = results["ids"][0]
        documents = results.get("documents", [[]])[0]
        metadatas = results.get("metadatas", [[]])[0]
        distances = results.get("distances", [[]])[0]

        for i, doc_id in enumerate(ids):
            # Convert distance to similarity score (cosine distance to similarity)
            similarity = 1 - distances[i] if i < len(distances) else 0

            formatted.append({
                "id": doc_id,
                "content": documents[i] if i < len(documents) else "",
                "metadata": metadatas[i] if i < len(metadatas) else {},
                "similarity": round(similarity, 4),
            })

        return formatted

    def get_collection_stats(self) -> dict[str, int]:
        """Get document counts for all collections."""
        stats = {}
        for name in [LITERATURE_COLLECTION, COMPOUNDS_COLLECTION]:
            try:
                collection = self.client.get_collection(name)
                stats[name] = collection.count()
            except Exception:
                stats[name] = 0
        return stats

    def clear_collection(self, name: str) -> bool:
        """Clear a collection (for re-indexing)."""
        try:
            self.client.delete_collection(name)
            return True
        except Exception:
            return False


# Global instance
rag_service = RAGService()
