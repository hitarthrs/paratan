"""
Component registry for the simple mirror machine.

Adding a new component type to Paratan 2 means:
  1. Create a new class that inherits from MachineComponent
  2. Decorate it with @register_component("your_name")
  3. Reference "your_name" in your machine config

No changes to SimpleMachineBuilder are needed.
"""

from __future__ import annotations
from paratan2.models.components.base import MachineComponent

_REGISTRY: dict[str, type[MachineComponent]] = {}


def register_component(name: str):
    """Decorator to register a MachineComponent class under a string key."""
    def decorator(cls: type[MachineComponent]) -> type[MachineComponent]:
        if name in _REGISTRY:
            raise ValueError(f"Component '{name}' is already registered.")
        _REGISTRY[name] = cls
        return cls
    return decorator


def get_component_class(name: str) -> type[MachineComponent]:
    if name not in _REGISTRY:
        raise KeyError(
            f"No component registered under '{name}'. "
            f"Available: {sorted(_REGISTRY.keys())}"
        )
    return _REGISTRY[name]


def list_components() -> list[str]:
    return sorted(_REGISTRY.keys())
