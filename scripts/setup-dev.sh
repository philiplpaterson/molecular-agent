#!/bin/bash
set -e

echo "🧬 Setting up MolecularAgent development environment..."

# Check if .env exists, if not copy from example
if [ ! -f .env ]; then
    echo "📋 Creating .env from .env.example..."
    cp .env.example .env
    echo "⚠️  Please update .env with your API keys before running the agent."
fi

# Build containers
echo "🔨 Building Docker containers..."
docker compose build

# Start services
echo "🚀 Starting services..."
docker compose up -d

# Wait for database to be ready
echo "⏳ Waiting for database to be ready..."
sleep 5

# Run migrations
echo "📊 Running database migrations..."
docker compose exec -T backend alembic upgrade head || echo "Note: Migrations may need to be generated first"

echo ""
echo "✅ Setup complete!"
echo ""
echo "🌐 Frontend: http://localhost:5173"
echo "🔧 Backend API: http://localhost:8000"
echo "📚 API Docs: http://localhost:8000/docs"
echo ""
echo "To view logs: docker compose logs -f"
echo "To stop: docker compose down"
