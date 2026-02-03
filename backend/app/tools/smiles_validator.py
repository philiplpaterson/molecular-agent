from typing import Any

from rdkit import Chem

from app.tools.base import BaseTool


class SMILESValidatorTool(BaseTool):
    """Tool to validate SMILES strings and convert to canonical form."""

    @property
    def name(self) -> str:
        return "validate_smiles"

    @property
    def description(self) -> str:
        return "Validate a SMILES string and return its canonical form. SMILES (Simplified Molecular Input Line Entry System) is a notation for representing molecular structures as text."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "string",
                    "description": "The SMILES string to validate",
                }
            },
            "required": ["smiles"],
        }

    def execute(self, smiles: str, **kwargs: Any) -> dict:
        """
        Validate a SMILES string and return canonical form.

        Args:
            smiles: The SMILES string to validate

        Returns:
            Dictionary with validation result, canonical SMILES, and any error
        """
        if not smiles or not smiles.strip():
            return {
                "valid": False,
                "canonical_smiles": None,
                "error": "Empty SMILES string provided",
            }

        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            if mol is None:
                return {
                    "valid": False,
                    "canonical_smiles": None,
                    "error": f"Invalid SMILES string: '{smiles}'",
                }

            canonical = Chem.MolToSmiles(mol, canonical=True)
            return {
                "valid": True,
                "canonical_smiles": canonical,
                "error": None,
            }
        except Exception as e:
            return {
                "valid": False,
                "canonical_smiles": None,
                "error": f"Error processing SMILES: {str(e)}",
            }
