import pytest
from rdkit import Chem

from app.tools.molecule_generator import MoleculeGeneratorTool


class TestMoleculeGenerator:
    @pytest.fixture
    def tool(self):
        return MoleculeGeneratorTool()

    def test_generate_default_count(self, tool):
        result = tool.execute(num_molecules=5)
        assert "error" not in result
        assert result["num_generated"] <= 5
        assert len(result["generated_smiles"]) == result["num_generated"]

    def test_generated_molecules_are_valid(self, tool):
        result = tool.execute(num_molecules=10)
        for smiles in result["generated_smiles"]:
            mol = Chem.MolFromSmiles(smiles)
            assert mol is not None, f"Invalid SMILES generated: {smiles}"

    def test_generate_with_seed(self, tool):
        result = tool.execute(num_molecules=3, seed_smiles="c1ccccc1")
        assert "error" not in result
        assert result["seed_smiles"] == "c1ccccc1"
        assert result["num_generated"] <= 3

    def test_limit_num_molecules(self, tool):
        # Test that max limit is enforced
        result = tool.execute(num_molecules=100)
        assert result["num_generated"] <= 20

    def test_unique_molecules(self, tool):
        result = tool.execute(num_molecules=10)
        # Check that all generated molecules are unique
        assert len(result["generated_smiles"]) == len(set(result["generated_smiles"]))

    def test_tool_properties(self, tool):
        assert tool.name == "generate_molecules"
        assert "generate" in tool.description.lower()
