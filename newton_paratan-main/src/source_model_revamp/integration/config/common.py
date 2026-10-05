"""Shared validation and coercion helpers for source model configuration mappings"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from typing import Any
import numpy as np

def canonical_model_name(value: str | None, aliases: Mapping[str, str], default: str) -> str:
    """
    Return the canonical model identifier for one configured name
    
    Alias matching is case insensitive and treats ASCII hyphen and underscore spellings as equivalent
    """
    if value is None:
        return default
    key = str(value).strip()
    lowered = key.lower().replace("-", "_")
    for alias, canonical in aliases.items():
        if lowered == alias.lower().replace("-", "_"):
            return canonical
        
    raise ValueError(f"Unsupported model name {value!r}, supported values are {sorted(set(aliases.values()))}")

def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    """Return a mapping value or an empty mapping when the input is absent"""
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    
    return value

def _check_unknown(mapping: Mapping[str, Any], allowed: set[str], name: str, strict: bool) -> None:
    """Reject configuration keys outside the allowed set when strict parsing is enabled"""
    unknown = set(mapping) - allowed
    if strict and unknown:
        raise ValueError(f"Unknown keys in {name}: {sorted(unknown)}")

def _number(value: Any, name: str, *, positive: bool = False, nonnegative: bool = False, default: float | None = None) -> float:
    """Coerce one finite floating point value and apply optional sign constraints"""
    if value is None:
        if default is None:
            raise ValueError(f"{name} is required")
        value = default
    scalar = float(value)
    if not np.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    if positive and scalar <= 0.0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and scalar < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    
    return scalar

def _integer(value: Any, name: str, *, minimum: int = 1, default: int | None = None) -> int:
    """Coerce one integer value and enforce the configured minimum"""
    if value is None:
        if default is None:
            raise ValueError(f"{name} is required")
        value = default
    integer = int(value)
    if integer < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
 
    return integer

def _bool(value: Any, default: bool = False) -> bool:
    """Coerce common textual boolean values before falling back to Python truth conversion"""
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        key = value.strip().lower()
        if key in {"true", "yes", "1", "on"}:
            return True
        if key in {"false", "no", "0", "off"}:
            return False
  
    return bool(value)

def _array(value: Any, name: str, *, default: Sequence[float] | None = None, positive: bool = False, nonnegative: bool = False) -> tuple[float, ...]:
    """Coerce one finite one dimensional numeric sequence and apply optional sign constraints"""
    if value is None:
        if default is None:
            raise ValueError(f"{name} is required")
        value = default
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D sequence")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    if positive and np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain only positive values")
    if nonnegative and np.any(arr < 0.0):
        raise ValueError(f"{name} must contain only nonnegative values")
  
    return tuple(float(x) for x in arr)

def _point3(value: Any, name: str) -> tuple[float, float, float]:
    """Return exactly three finite coordinates as a tuple"""
    arr = _array(value, name)
    if len(arr) != 3:
        raise ValueError(f"{name} must contain exactly 3 coordinates")
 
    return (arr[0], arr[1], arr[2])

def _sum_thickness_cm(items: Any) -> float:
    """Sum `thickness` values from mapping or sequence style geometry layer definitions in cm"""
    if items is None:
        return 0.0
    if isinstance(items, Mapping):
        return float(sum(float(v.get("thickness", 0.0)) for v in items.values() if isinstance(v, Mapping)))
    if isinstance(items, Sequence) and not isinstance(items, (str, bytes)):
        return float(sum(float(v.get("thickness", 0.0)) for v in items if isinstance(v, Mapping)))
   
    return 0.0

def _coerce_float_list(value: Any, default: tuple[float, ...] = ()) -> tuple[float, ...]:
    """Return a scalar or sequence input as a tuple of floats"""
    if value is None:
        return tuple(default)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return tuple(float(x) for x in value)
   
    return (float(value),)