from app.tools.base import BaseTool
from app.tools.compound_search import CompoundSearchTool
from app.tools.drug_likeness import DrugLikenessTool
from app.tools.literature_search import LiteratureSearchTool
from app.tools.molecule_generator import MoleculeGeneratorTool
from app.tools.property_predictor import PropertyPredictorTool
from app.tools.similarity_search import SimilaritySearchTool
from app.tools.smiles_validator import SMILESValidatorTool

__all__ = [
    "BaseTool",
    "SMILESValidatorTool",
    "PropertyPredictorTool",
    "DrugLikenessTool",
    "SimilaritySearchTool",
    "MoleculeGeneratorTool",
    # RAG tools
    "LiteratureSearchTool",
    "CompoundSearchTool",
]

AVAILABLE_TOOLS: dict[str, type[BaseTool]] = {
    "validate_smiles": SMILESValidatorTool,
    "predict_properties": PropertyPredictorTool,
    "check_drug_likeness": DrugLikenessTool,
    "similarity_search": SimilaritySearchTool,
    "generate_molecules": MoleculeGeneratorTool,
    # RAG tools
    "search_literature": LiteratureSearchTool,
    "search_compounds": CompoundSearchTool,
}
