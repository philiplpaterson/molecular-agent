from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from app.tools.base import BaseTool


class SimilaritySearchTool(BaseTool):
    """Tool to calculate molecular similarity using Tanimoto coefficient."""

    @property
    def name(self) -> str:
        return "similarity_search"

    @property
    def description(self) -> str:
        return "Calculate Tanimoto similarity between a query molecule and a list of target molecules. Returns similarity scores from 0 (no similarity) to 1 (identical). Uses Morgan fingerprints (ECFP4-like)."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "query_smiles": {
                    "type": "string",
                    "description": "The SMILES string of the query molecule",
                },
                "target_smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings to compare against",
                },
            },
            "required": ["query_smiles", "target_smiles"],
        }

    def execute(self, query_smiles: str, target_smiles: list[str], **kwargs: Any) -> dict:
        """
        Calculate Tanimoto similarity between query and target molecules.

        Args:
            query_smiles: The SMILES string of the query molecule
            target_smiles: List of SMILES strings to compare against

        Returns:
            Dictionary with similarity scores sorted by similarity (descending)
        """
        query_mol = Chem.MolFromSmiles(query_smiles.strip())
        if query_mol is None:
            return {"error": f"Invalid query SMILES: '{query_smiles}'"}

        try:
            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)

            similarities = []
            for target in target_smiles:
                target_mol = Chem.MolFromSmiles(target.strip())
                if target_mol is None:
                    similarities.append({"smiles": target, "similarity": 0.0, "error": "Invalid SMILES"})
                    continue

                target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, radius=2, nBits=2048)
                similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
                similarities.append({"smiles": target, "similarity": round(similarity, 4)})

            # Sort by similarity descending
            similarities.sort(key=lambda x: x.get("similarity", 0), reverse=True)

            return {
                "query_smiles": query_smiles,
                "num_targets": len(target_smiles),
                "similarities": similarities,
            }
        except Exception as e:
            return {"error": f"Error calculating similarity: {str(e)}"}
