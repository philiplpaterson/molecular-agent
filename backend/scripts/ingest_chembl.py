#!/usr/bin/env python3
"""
Ingest ChEMBL compound data into ChromaDB.

Uses ChEMBL web services API to fetch approved drugs and clinical candidates.
"""

import argparse
import sys
import time
from pathlib import Path

import httpx

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.rag import rag_service

CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"


def fetch_approved_drugs(limit: int = 1000, offset: int = 0) -> list[dict]:
    """Fetch approved drugs from ChEMBL."""
    params = {
        "max_phase": 4,  # Approved drugs
        "limit": min(limit, 1000),
        "offset": offset,
        "format": "json",
    }

    response = httpx.get(
        f"{CHEMBL_BASE_URL}/molecule",
        params=params,
        timeout=60,
    )
    response.raise_for_status()

    data = response.json()
    molecules = data.get("molecules", [])

    compounds = []
    for mol in molecules:
        try:
            # Get SMILES from molecule_structures
            structures = mol.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles", "")

            if not smiles:
                continue

            # Get therapeutic areas from indication_class
            indication = mol.get("indication_class", "") or ""

            compounds.append({
                "chembl_id": mol.get("molecule_chembl_id", ""),
                "name": mol.get("pref_name", "") or "",
                "smiles": smiles,
                "description": mol.get("molecule_type", "") or "",
                "max_phase": mol.get("max_phase", 0),
                "molecule_type": mol.get("molecule_type", "") or "",
                "therapeutic_areas": indication,
            })

        except Exception as e:
            print(f"Error parsing molecule: {e}")
            continue

    return compounds


def fetch_clinical_candidates(min_phase: int = 2, limit: int = 1000) -> list[dict]:
    """Fetch clinical candidates from ChEMBL."""
    params = {
        "max_phase__gte": min_phase,
        "max_phase__lt": 4,  # Not yet approved
        "limit": min(limit, 1000),
        "format": "json",
    }

    response = httpx.get(
        f"{CHEMBL_BASE_URL}/molecule",
        params=params,
        timeout=60,
    )
    response.raise_for_status()

    data = response.json()
    molecules = data.get("molecules", [])

    compounds = []
    for mol in molecules:
        try:
            structures = mol.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles", "")

            if not smiles:
                continue

            compounds.append({
                "chembl_id": mol.get("molecule_chembl_id", ""),
                "name": mol.get("pref_name", "") or "",
                "smiles": smiles,
                "description": f"Phase {mol.get('max_phase', '')} clinical candidate",
                "max_phase": mol.get("max_phase", 0),
                "molecule_type": mol.get("molecule_type", "") or "",
                "therapeutic_areas": mol.get("indication_class", "") or "",
            })

        except Exception as e:
            print(f"Error parsing molecule: {e}")
            continue

    return compounds


def ingest_chembl_data(
    max_approved: int = 1000,
    max_clinical: int = 500,
    batch_size: int = 100,
) -> int:
    """
    Ingest ChEMBL compounds into ChromaDB.

    Args:
        max_approved: Maximum approved drugs to fetch
        max_clinical: Maximum clinical candidates to fetch
        batch_size: Number of compounds to index at once

    Returns:
        Total number of compounds indexed
    """
    total_indexed = 0

    # Fetch approved drugs
    print(f"Fetching up to {max_approved} approved drugs from ChEMBL...")
    approved = []
    offset = 0

    while len(approved) < max_approved:
        batch = fetch_approved_drugs(limit=min(1000, max_approved - len(approved)), offset=offset)
        if not batch:
            break
        approved.extend(batch)
        offset += len(batch)
        print(f"  Fetched {len(approved)} approved drugs...")
        time.sleep(0.5)

    print(f"Total approved drugs fetched: {len(approved)}")

    # Index approved drugs
    for i in range(0, len(approved), batch_size):
        batch = approved[i : i + batch_size]
        indexed = rag_service.index_compounds(batch)
        total_indexed += indexed
        print(f"  Indexed batch {i // batch_size + 1}: {indexed} compounds")

    # Fetch clinical candidates
    print(f"\nFetching up to {max_clinical} clinical candidates from ChEMBL...")
    clinical = fetch_clinical_candidates(min_phase=2, limit=max_clinical)
    print(f"Total clinical candidates fetched: {len(clinical)}")

    # Index clinical candidates
    for i in range(0, len(clinical), batch_size):
        batch = clinical[i : i + batch_size]
        indexed = rag_service.index_compounds(batch)
        total_indexed += indexed
        print(f"  Indexed batch {i // batch_size + 1}: {indexed} compounds")

    return total_indexed


def main():
    parser = argparse.ArgumentParser(description="Ingest ChEMBL compounds into ChromaDB")
    parser.add_argument("--max-approved", type=int, default=1000, help="Max approved drugs to fetch")
    parser.add_argument("--max-clinical", type=int, default=500, help="Max clinical candidates to fetch")
    parser.add_argument("--clear", action="store_true", help="Clear existing data before ingesting")

    args = parser.parse_args()

    if args.clear:
        print("Clearing existing compounds collection...")
        rag_service.clear_collection("compounds")

    print("Starting ChEMBL ingestion...")
    total = ingest_chembl_data(
        max_approved=args.max_approved,
        max_clinical=args.max_clinical,
    )

    print(f"\nIngestion complete! Total compounds indexed: {total}")


if __name__ == "__main__":
    main()
