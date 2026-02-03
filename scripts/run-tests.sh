#!/bin/bash
set -e

echo "🧪 Running MolecularAgent tests..."

# Run backend tests
echo "Running backend tests..."
docker compose exec -T backend pytest -v --cov=app --cov-report=term-missing

echo ""
echo "✅ All tests passed!"
