from fastapi import APIRouter

from app.api.routes import agent, health, molecules

api_router = APIRouter()
api_router.include_router(health.router, prefix="/health", tags=["health"])
api_router.include_router(molecules.router, prefix="/molecules", tags=["molecules"])
api_router.include_router(agent.router, prefix="/agent", tags=["agent"])
