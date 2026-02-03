from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_validate_valid_smiles():
    response = client.post("/api/v1/molecules/validate", json={"smiles": "CCO"})
    assert response.status_code == 200
    data = response.json()
    assert data["valid"] is True
    assert data["canonical_smiles"] == "CCO"


def test_validate_invalid_smiles():
    response = client.post("/api/v1/molecules/validate", json={"smiles": "invalid"})
    assert response.status_code == 200
    data = response.json()
    assert data["valid"] is False


def test_get_properties():
    response = client.post("/api/v1/molecules/properties", json={"smiles": "CCO"})
    assert response.status_code == 200
    data = response.json()
    assert "molecular_weight" in data
    assert "logp" in data
    assert "tpsa" in data


def test_drug_likeness():
    response = client.post("/api/v1/molecules/drug-likeness", json={"smiles": "CCO"})
    assert response.status_code == 200
    data = response.json()
    assert "passes_lipinski" in data
    assert "violations" in data
    assert "properties" in data


def test_similarity():
    response = client.post(
        "/api/v1/molecules/similarity",
        json={"query_smiles": "CCO", "target_smiles": ["CCO", "CCCO"]},
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    assert data[0]["similarity"] == 1.0  # Identical molecule
