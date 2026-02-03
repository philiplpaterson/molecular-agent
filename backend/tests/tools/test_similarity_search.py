import pytest

from app.tools.similarity_search import SimilaritySearchTool


class TestSimilaritySearch:
    @pytest.fixture
    def tool(self):
        return SimilaritySearchTool()

    def test_identical_molecules(self, tool):
        result = tool.execute(query_smiles="CCO", target_smiles=["CCO"])
        assert "error" not in result
        assert len(result["similarities"]) == 1
        assert result["similarities"][0]["similarity"] == 1.0

    def test_different_molecules(self, tool):
        result = tool.execute(
            query_smiles="CCO",  # ethanol
            target_smiles=["c1ccccc1"],  # benzene
        )
        assert "error" not in result
        assert len(result["similarities"]) == 1
        assert result["similarities"][0]["similarity"] < 1.0

    def test_multiple_targets(self, tool):
        result = tool.execute(
            query_smiles="CCO",
            target_smiles=["CCO", "CCCO", "c1ccccc1"],
        )
        assert "error" not in result
        assert len(result["similarities"]) == 3
        # Results should be sorted by similarity descending
        similarities = [r["similarity"] for r in result["similarities"]]
        assert similarities == sorted(similarities, reverse=True)

    def test_invalid_query_smiles(self, tool):
        result = tool.execute(query_smiles="invalid", target_smiles=["CCO"])
        assert "error" in result

    def test_invalid_target_smiles(self, tool):
        result = tool.execute(query_smiles="CCO", target_smiles=["invalid"])
        assert "error" not in result
        assert result["similarities"][0]["error"] == "Invalid SMILES"

    def test_tool_properties(self, tool):
        assert tool.name == "similarity_search"
        assert "similarity" in tool.description.lower()
