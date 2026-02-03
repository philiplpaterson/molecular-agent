import pytest

from app.tools.drug_likeness import DrugLikenessTool


class TestDrugLikeness:
    @pytest.fixture
    def tool(self):
        return DrugLikenessTool()

    def test_drug_like_aspirin(self, tool):
        # Aspirin passes Lipinski's rules
        result = tool.execute(smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert "error" not in result
        assert result["passes_lipinski"] is True
        assert result["num_violations"] == 0

    def test_drug_like_ethanol(self, tool):
        result = tool.execute(smiles="CCO")
        assert "error" not in result
        assert result["passes_lipinski"] is True

    def test_properties_returned(self, tool):
        result = tool.execute(smiles="CCO")
        assert "properties" in result
        props = result["properties"]
        assert "molecular_weight" in props
        assert "logp" in props
        assert "hbd" in props
        assert "hba" in props

    def test_invalid_smiles(self, tool):
        result = tool.execute(smiles="invalid")
        assert "error" in result

    def test_tool_properties(self, tool):
        assert tool.name == "check_drug_likeness"
        assert "lipinski" in tool.description.lower()
