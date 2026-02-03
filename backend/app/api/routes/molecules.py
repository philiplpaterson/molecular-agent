from fastapi import APIRouter
from pydantic import BaseModel

router = APIRouter()


class SMILESInput(BaseModel):
    smiles: str


class SMILESListInput(BaseModel):
    query_smiles: str
    target_smiles: list[str]


class ValidationResult(BaseModel):
    valid: bool
    canonical_smiles: str | None = None
    error: str | None = None


class PropertyResult(BaseModel):
    molecular_weight: float
    logp: float
    tpsa: float
    rotatable_bonds: int
    hbd: int
    hba: int


class DrugLikenessResult(BaseModel):
    passes_lipinski: bool
    violations: list[str]
    properties: dict


class SimilarityResult(BaseModel):
    smiles: str
    similarity: float


@router.post("/validate", response_model=ValidationResult)
def validate_smiles(data: SMILESInput):
    from app.tools.smiles_validator import SMILESValidatorTool

    tool = SMILESValidatorTool()
    result = tool.execute(smiles=data.smiles)
    return ValidationResult(**result)


@router.post("/properties", response_model=PropertyResult)
def get_properties(data: SMILESInput):
    from app.tools.property_predictor import PropertyPredictorTool

    tool = PropertyPredictorTool()
    result = tool.execute(smiles=data.smiles)
    return PropertyResult(**result)


@router.post("/drug-likeness", response_model=DrugLikenessResult)
def check_drug_likeness(data: SMILESInput):
    from app.tools.drug_likeness import DrugLikenessTool

    tool = DrugLikenessTool()
    result = tool.execute(smiles=data.smiles)
    return DrugLikenessResult(**result)


@router.post("/similarity", response_model=list[SimilarityResult])
def calculate_similarity(data: SMILESListInput):
    from app.tools.similarity_search import SimilaritySearchTool

    tool = SimilaritySearchTool()
    result = tool.execute(query_smiles=data.query_smiles, target_smiles=data.target_smiles)
    return [SimilarityResult(**r) for r in result["similarities"]]
