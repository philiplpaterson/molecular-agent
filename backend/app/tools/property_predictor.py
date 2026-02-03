from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from app.tools.base import BaseTool


class PropertyPredictorTool(BaseTool):
    """Tool to predict molecular properties from SMILES."""

    @property
    def name(self) -> str:
        return "predict_properties"

    @property
    def description(self) -> str:
        return "Predict molecular properties from a SMILES string including molecular weight, LogP, topological polar surface area (TPSA), rotatable bonds, and hydrogen bond donors/acceptors."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "string",
                    "description": "The SMILES string of the molecule",
                }
            },
            "required": ["smiles"],
        }

    def execute(self, smiles: str, **kwargs: Any) -> dict:
        """
        Predict molecular properties from SMILES.

        Args:
            smiles: The SMILES string of the molecule

        Returns:
            Dictionary with molecular properties
        """
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return {"error": f"Invalid SMILES string: '{smiles}'"}

        try:
            return {
                "molecular_weight": round(Descriptors.MolWt(mol), 2),
                "logp": round(Descriptors.MolLogP(mol), 2),
                "tpsa": round(Descriptors.TPSA(mol), 2),
                "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                "hbd": rdMolDescriptors.CalcNumHBD(mol),  # Hydrogen bond donors
                "hba": rdMolDescriptors.CalcNumHBA(mol),  # Hydrogen bond acceptors
                "num_atoms": mol.GetNumAtoms(),
                "num_heavy_atoms": mol.GetNumHeavyAtoms(),
                "num_rings": rdMolDescriptors.CalcNumRings(mol),
                "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            }
        except Exception as e:
            return {"error": f"Error calculating properties: {str(e)}"}
