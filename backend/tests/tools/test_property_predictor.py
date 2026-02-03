import pytest

from app.tools.property_predictor import PropertyPredictorTool


class TestPropertyPredictor:
    @pytest.fixture
    def tool(self):
        return PropertyPredictorTool()

    def test_ethanol_properties(self, tool):
        result = tool.execute(smiles="CCO")
        assert "error" not in result
        assert result["molecular_weight"] == pytest.approx(46.07, rel=0.01)
        assert result["hbd"] == 1  # One OH group
        assert result["hba"] == 1  # One oxygen

    def test_aspirin_properties(self, tool):
        # Aspirin: CC(=O)Oc1ccccc1C(=O)O
        result = tool.execute(smiles="CC(=O)Oc1ccccc1C(=O)O")
        assert "error" not in result
        assert result["molecular_weight"] == pytest.approx(180.16, rel=0.01)
        assert result["num_aromatic_rings"] == 1

    def test_invalid_smiles(self, tool):
        result = tool.execute(smiles="invalid")
        assert "error" in result

    def test_all_properties_present(self, tool):
        result = tool.execute(smiles="CCO")
        expected_keys = [
            "molecular_weight",
            "logp",
            "tpsa",
            "rotatable_bonds",
            "hbd",
            "hba",
            "num_atoms",
            "num_heavy_atoms",
            "num_rings",
            "num_aromatic_rings",
        ]
        for key in expected_keys:
            assert key in result

    def test_tool_properties(self, tool):
        assert tool.name == "predict_properties"
        assert "properties" in tool.description.lower()
