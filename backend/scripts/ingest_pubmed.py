#!/usr/bin/env python3
"""
Ingest PubMed abstracts into ChromaDB.

Uses NCBI E-utilities API to fetch abstracts related to drug discovery.
"""

import argparse
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import httpx

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.rag import rag_service

NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Drug discovery related search terms
DEFAULT_SEARCH_TERMS = [
    "drug discovery",
    "molecular docking",
    "QSAR",
    "pharmacokinetics",
    "drug target",
    "lead optimization",
    "ADMET",
    "drug repurposing",
    "structure activity relationship",
    "high throughput screening",
]


def search_pubmed(query: str, max_results: int = 100) -> list[str]:
    """Search PubMed and return PMIDs."""
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "sort": "relevance",
    }

    response = httpx.get(f"{NCBI_BASE_URL}/esearch.fcgi", params=params, timeout=30)
    response.raise_for_status()

    data = response.json()
    return data.get("esearchresult", {}).get("idlist", [])


def fetch_abstracts(pmids: list[str]) -> list[dict]:
    """Fetch abstract details for given PMIDs."""
    if not pmids:
        return []

    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
    }

    response = httpx.get(f"{NCBI_BASE_URL}/efetch.fcgi", params=params, timeout=60)
    response.raise_for_status()

    # Parse XML response
    root = ET.fromstring(response.content)
    articles = []

    for article in root.findall(".//PubmedArticle"):
        try:
            medline = article.find(".//MedlineCitation")
            if medline is None:
                continue

            pmid_elem = medline.find(".//PMID")
            pmid = pmid_elem.text if pmid_elem is not None else ""

            article_elem = medline.find(".//Article")
            if article_elem is None:
                continue

            title_elem = article_elem.find(".//ArticleTitle")
            title = title_elem.text if title_elem is not None else ""

            abstract_elem = article_elem.find(".//Abstract/AbstractText")
            abstract = abstract_elem.text if abstract_elem is not None else ""

            # Get authors
            authors = []
            for author in article_elem.findall(".//Author"):
                lastname = author.find("LastName")
                forename = author.find("ForeName")
                if lastname is not None:
                    name = lastname.text
                    if forename is not None:
                        name = f"{forename.text} {name}"
                    authors.append(name)

            # Get journal and year
            journal_elem = article_elem.find(".//Journal/Title")
            journal = journal_elem.text if journal_elem is not None else ""

            year_elem = article_elem.find(".//Journal/JournalIssue/PubDate/Year")
            year = int(year_elem.text) if year_elem is not None and year_elem.text else 0

            if title and abstract:
                articles.append({
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "authors": ", ".join(authors[:5]),  # First 5 authors
                    "journal": journal,
                    "year": year,
                })

        except Exception as e:
            print(f"Error parsing article: {e}")
            continue

    return articles


def ingest_pubmed_data(
    search_terms: list[str] | None = None,
    max_per_term: int = 100,
    batch_size: int = 50,
) -> int:
    """
    Ingest PubMed abstracts into ChromaDB.

    Args:
        search_terms: List of search terms (uses defaults if None)
        max_per_term: Maximum results per search term
        batch_size: Number of PMIDs to fetch at once

    Returns:
        Total number of documents indexed
    """
    if search_terms is None:
        search_terms = DEFAULT_SEARCH_TERMS

    total_indexed = 0
    all_pmids = set()

    print(f"Searching PubMed with {len(search_terms)} terms...")

    for term in search_terms:
        print(f"  Searching: {term}")
        pmids = search_pubmed(term, max_per_term)
        all_pmids.update(pmids)
        time.sleep(0.5)  # Rate limiting

    print(f"Found {len(all_pmids)} unique PMIDs")

    # Fetch and index in batches
    pmid_list = list(all_pmids)
    for i in range(0, len(pmid_list), batch_size):
        batch = pmid_list[i : i + batch_size]
        print(f"Fetching batch {i // batch_size + 1}/{(len(pmid_list) + batch_size - 1) // batch_size}...")

        articles = fetch_abstracts(batch)
        if articles:
            indexed = rag_service.index_literature(articles)
            total_indexed += indexed
            print(f"  Indexed {indexed} articles")

        time.sleep(1)  # Rate limiting

    return total_indexed


def main():
    parser = argparse.ArgumentParser(description="Ingest PubMed abstracts into ChromaDB")
    parser.add_argument("--max-per-term", type=int, default=100, help="Max results per search term")
    parser.add_argument("--terms", nargs="+", help="Custom search terms")
    parser.add_argument("--clear", action="store_true", help="Clear existing data before ingesting")

    args = parser.parse_args()

    if args.clear:
        print("Clearing existing literature collection...")
        rag_service.clear_collection("literature")

    print("Starting PubMed ingestion...")
    total = ingest_pubmed_data(
        search_terms=args.terms,
        max_per_term=args.max_per_term,
    )

    print(f"\nIngestion complete! Total documents indexed: {total}")


if __name__ == "__main__":
    main()
