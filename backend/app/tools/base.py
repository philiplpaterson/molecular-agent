from abc import ABC, abstractmethod
from typing import Any


class BaseTool(ABC):
    """Base class for all molecular tools."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Tool name used for identification."""
        pass

    @property
    @abstractmethod
    def description(self) -> str:
        """Description of what the tool does."""
        pass

    @property
    @abstractmethod
    def parameters(self) -> dict:
        """JSON schema for tool parameters."""
        pass

    @abstractmethod
    def execute(self, **kwargs: Any) -> dict:
        """Execute the tool with given parameters."""
        pass
