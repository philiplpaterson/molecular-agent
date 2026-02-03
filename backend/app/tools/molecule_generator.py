import random
from typing import Any

from rdkit import Chem
from rdkit.Chem import AllChem

from app.tools.base import BaseTool

# Common drug-like molecular fragments for generation
FRAGMENTS = [
    "c1ccccc1",  # benzene
    "c1ccncc1",  # pyridine
    "c1ccc2[nH]ccc2c1",  # indole
    "c1cnc2ccccc2n1",  # quinazoline
    "C1CCNCC1",  # piperidine
    "C1CCOCC1",  # tetrahydropyran
    "C1CCNC1",  # pyrrolidine
    "c1ccc2c(c1)cccc2",  # naphthalene
    "c1csc2ccccc12",  # benzothiophene
    "c1coc2ccccc12",  # benzofuran
]

LINKERS = ["C", "CC", "CCC", "O", "N", "CO", "CN", "CCO", "CCN", "C(=O)N", "C(=O)O"]

SUBSTITUENTS = ["C", "CC", "CCC", "O", "N", "F", "Cl", "Br", "OC", "NC", "C(=O)C", "C#N", "C(F)(F)F"]


class MoleculeGeneratorTool(BaseTool):
    """Tool to generate random drug-like molecules."""

    @property
    def name(self) -> str:
        return "generate_molecules"

    @property
    def description(self) -> str:
        return "Generate random drug-like molecules. Can generate molecules similar to a seed structure or completely random molecules using common drug-like fragments."

    @property
    def parameters(self) -> dict:
        return {
            "type": "object",
            "properties": {
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of molecules to generate (1-20)",
                    "minimum": 1,
                    "maximum": 20,
                },
                "seed_smiles": {
                    "type": "string",
                    "description": "Optional seed SMILES to generate similar molecules",
                },
            },
            "required": ["num_molecules"],
        }

    def _generate_random_molecule(self) -> str | None:
        """Generate a random molecule from fragments."""
        # Select 1-2 core fragments
        num_cores = random.randint(1, 2)
        cores = random.sample(FRAGMENTS, min(num_cores, len(FRAGMENTS)))

        # Build molecule
        if len(cores) == 1:
            smiles = cores[0]
        else:
            linker = random.choice(LINKERS)
            smiles = f"{cores[0]}{linker}{cores[1]}"

        # Add 0-2 substituents
        num_subs = random.randint(0, 2)
        for _ in range(num_subs):
            sub = random.choice(SUBSTITUENTS)
            smiles = f"{smiles}{sub}"

        # Validate the generated SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True)
        return None

    def _mutate_molecule(self, seed_smiles: str) -> str | None:
        """Generate a molecule similar to the seed by mutation."""
        mol = Chem.MolFromSmiles(seed_smiles)
        if mol is None:
            return None

        try:
            # Try to add a random substituent or fragment
            sub = random.choice(SUBSTITUENTS + FRAGMENTS[:3])
            new_smiles = f"{seed_smiles}{sub}"
            new_mol = Chem.MolFromSmiles(new_smiles)
            if new_mol is not None:
                return Chem.MolToSmiles(new_mol, canonical=True)
        except Exception:
            pass

        # Fallback: return a random molecule
        return self._generate_random_molecule()

    def execute(self, num_molecules: int = 5, seed_smiles: str | None = None, **kwargs: Any) -> dict:
        """
        Generate random drug-like molecules.

        Args:
            num_molecules: Number of molecules to generate (1-20)
            seed_smiles: Optional seed SMILES to generate similar molecules

        Returns:
            Dictionary with list of generated SMILES
        """
        num_molecules = max(1, min(20, num_molecules))

        generated = []
        attempts = 0
        max_attempts = num_molecules * 10  # Limit attempts to avoid infinite loops

        while len(generated) < num_molecules and attempts < max_attempts:
            attempts += 1

            if seed_smiles:
                smiles = self._mutate_molecule(seed_smiles)
            else:
                smiles = self._generate_random_molecule()

            if smiles and smiles not in generated:
                generated.append(smiles)

        return {
            "num_generated": len(generated),
            "seed_smiles": seed_smiles,
            "generated_smiles": generated,
        }
