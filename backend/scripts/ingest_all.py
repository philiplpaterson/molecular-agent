#!/usr/bin/env python3
"""
Run all data ingestion scripts.

This script can be run as a Docker job to populate the vector database.
"""

import argparse
import sys
import time
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.rag import rag_service


def wait_for_chromadb(max_retries: int = 30, delay: int = 2) -> bool:
    """Wait for ChromaDB to be available."""
    print("Waiting for ChromaDB to be available...")
    for i in range(max_retries):
        try:
            rag_service.client.heartbeat()
            print("ChromaDB is ready!")
            return True
        except Exception as e:
            print(f"  Attempt {i + 1}/{max_retries}: {e}")
            time.sleep(delay)
    return False


def main():
    parser = argparse.ArgumentParser(description="Run all data ingestion scripts")
    parser.add_argument("--clear", action="store_true", help="Clear existing data before ingesting")
    parser.add_argument("--skip-pubmed", action="store_true", help="Skip PubMed ingestion")
    parser.add_argument("--skip-chembl", action="store_true", help="Skip ChEMBL ingestion")
    parser.add_argument("--pubmed-max", type=int, default=100, help="Max PubMed results per term")
    parser.add_argument("--chembl-max", type=int, default=1500, help="Max ChEMBL approved drugs")

    args = parser.parse_args()

    # Wait for ChromaDB
    if not wait_for_chromadb():
        print("ERROR: ChromaDB is not available. Exiting.")
        sys.exit(1)

    # Clear collections if requested
    if args.clear:
        print("\nClearing all collections...")
        rag_service.clear_collection("literature")
        rag_service.clear_collection("compounds")
        print("Collections cleared.")

    total_indexed = 0

    # Ingest PubMed
    if not args.skip_pubmed:
        print("\n" + "=" * 50)
        print("INGESTING PUBMED ABSTRACTS")
        print("=" * 50)
        from scripts.ingest_pubmed import ingest_pubmed_data

        count = ingest_pubmed_data(max_per_term=args.pubmed_max)
        total_indexed += count
        print(f"Indexed {count} literature documents")

    # Ingest ChEMBL
    if not args.skip_chembl:
        print("\n" + "=" * 50)
        print("INGESTING CHEMBL COMPOUNDS")
        print("=" * 50)
        from scripts.ingest_chembl import ingest_chembl_data

        count = ingest_chembl_data(max_approved=args.chembl_max)
        total_indexed += count
        print(f"Indexed {count} compounds")

    # Final summary
    print("\n" + "=" * 50)
    print("INGESTION COMPLETE")
    print("=" * 50)
    print(f"Total documents indexed: {total_indexed}")

    stats = rag_service.get_collection_stats()
    print("\nCollection statistics:")
    for name, count in stats.items():
        print(f"  {name}: {count} documents")


if __name__ == "__main__":
    main()
