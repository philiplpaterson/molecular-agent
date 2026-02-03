from fastapi import APIRouter

router = APIRouter()


@router.get("")
def health_check():
    return {"status": "healthy"}


@router.get("/ready")
def readiness_check():
    # TODO: Add database and Redis health checks
    return {"status": "ready", "database": "connected", "redis": "connected"}
