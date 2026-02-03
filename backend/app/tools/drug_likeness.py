from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from app.tools.base import BaseTool


class DrugLikenessTool(BaseTool):
    """Tool to check Lipinski's Rule of Five for drug-likeness."""

    @property
    def name(self) -> str:
        return "check_drug_likeness"

    @property
    def description(self) -> str:
        return "Check if a molecule passes Lipinski's Rule of Five for drug-likeness. The rules are: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10. Molecules violating more than one rule are unlikely to be orally active drugs."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "string",
                    "description": "The SMILES string of the molecule to check",
                }
            },
            "required": ["smiles"],
        }

    def execute(self, smiles: str, **kwargs: Any) -> dict:
        """
        Check Lipinski's Rule of Five.

        Args:
            smiles: The SMILES string of the molecule

        Returns:
            Dictionary with drug-likeness assessment
        """
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return {"error": f"Invalid SMILES string: '{smiles}'"}

        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)

            violations = []

            if mw > 500:
                violations.append(f"Molecular weight ({mw:.1f}) > 500")
            if logp > 5:
                violations.append(f"LogP ({logp:.2f}) > 5")
            if hbd > 5:
                violations.append(f"Hydrogen bond donors ({hbd}) > 5")
            if hba > 10:
                violations.append(f"Hydrogen bond acceptors ({hba}) > 10")

            passes = len(violations) <= 1

            return {
                "passes_lipinski": passes,
                "violations": violations,
                "num_violations": len(violations),
                "properties": {
                    "molecular_weight": round(mw, 2),
                    "logp": round(logp, 2),
                    "hbd": hbd,
                    "hba": hba,
                },
                "assessment": "Drug-like" if passes else "Not drug-like (too many violations)",
            }
        except Exception as e:
            return {"error": f"Error checking drug-likeness: {str(e)}"}
