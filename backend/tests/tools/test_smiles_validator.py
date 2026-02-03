import pytest

from app.tools.smiles_validator import SMILESValidatorTool


class TestSMILESValidator:
    @pytest.fixture
    def tool(self):
        return SMILESValidatorTool()

    def test_valid_smiles_ethanol(self, tool):
        result = tool.execute(smiles="CCO")
        assert result["valid"] is True
        assert result["canonical_smiles"] == "CCO"
        assert result["error"] is None

    def test_valid_smiles_aspirin(self, tool):
        result = tool.execute(smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert result["valid"] is True
        assert result["canonical_smiles"] is not None
        assert result["error"] is None

    def test_valid_smiles_benzene(self, tool):
        result = tool.execute(smiles="c1ccccc1")
        assert result["valid"] is True
        assert result["canonical_smiles"] == "c1ccccc1"

    def test_invalid_smiles(self, tool):
        result = tool.execute(smiles="invalid_smiles")
        assert result["valid"] is False
        assert result["canonical_smiles"] is None
        assert "Invalid SMILES" in result["error"]

    def test_empty_smiles(self, tool):
        result = tool.execute(smiles="")
        assert result["valid"] is False
        assert result["canonical_smiles"] is None
        assert "Empty SMILES" in result["error"]

    def test_whitespace_smiles(self, tool):
        result = tool.execute(smiles="  CCO  ")
        assert result["valid"] is True
        assert result["canonical_smiles"] == "CCO"

    def test_tool_properties(self, tool):
        assert tool.name == "validate_smiles"
        assert "validate" in tool.description.lower()
        assert "smiles" in tool.parameters["properties"]
