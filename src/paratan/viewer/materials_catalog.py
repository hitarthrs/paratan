"""Source-agnostic material catalog for the viewer materials panel.

Architecture
------------
* ``MaterialRecord`` / ``CompositionLine`` are plain, JSON-friendly DTOs.
  The UI and Trame state never hold ``openmc.Material`` objects.
* ``MaterialCatalog`` is the protocol every backend must implement.
* ``OpenMCMaterialsModuleCatalog`` is today's adapter: it introspects
  ``paratan.materials.material`` (``materials_list`` + module attributes).
  A later YAML / XML / simulation-manifest backend can implement the same
  protocol without touching the UI.

Lookup keys are stable strings (``slug``). Display names come from the
OpenMC material ``name`` when present.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Iterable, Protocol, runtime_checkable


@dataclass(frozen=True)
class CompositionLine:
    """One nuclide or element entry in a material definition."""

    symbol: str
    fraction: float
    fraction_type: str  # OpenMC percent_type: ao, wo, ...


@dataclass(frozen=True)
class MaterialRecord:
    """Normalized material description, independent of the data source."""

    key: str
    display_name: str
    density: float | None
    density_units: str | None
    openmc_id: int | None
    composition: tuple[CompositionLine, ...]
    source: str
    aliases: tuple[str, ...] = ()
    notes: str = ""

    def to_summary(self) -> dict[str, Any]:
        dens = None
        if self.density is not None:
            unit = self.density_units or ""
            dens = f"{self.density:g} {unit}".strip()
        return {
            "key": self.key,
            "label": self.display_name,
            "density": dens or "—",
            "n_lines": len(self.composition),
            "openmc_id": self.openmc_id,
        }

    def to_detail(self) -> dict[str, Any]:
        lines = sorted(self.composition, key=lambda c: c.fraction, reverse=True)
        return {
            "key": self.key,
            "label": self.display_name,
            "source": self.source,
            "density": self.density,
            "density_units": self.density_units or "",
            "openmc_id": self.openmc_id,
            "notes": self.notes,
            "composition": [
                {
                    "symbol": line.symbol,
                    "fraction": line.fraction,
                    "fraction_type": line.fraction_type,
                    "percent": 100.0 * line.fraction,
                }
                for line in lines
            ],
        }


@runtime_checkable
class MaterialCatalog(Protocol):
    """Read-only material library used by the viewer."""

    @property
    def source_label(self) -> str: ...

    def ensure_loaded(self) -> None:
        """Eagerly build the index (may import heavy modules)."""

    def list_records(self) -> list[MaterialRecord]: ...

    def get(self, key: str) -> MaterialRecord | None: ...

    def resolve(self, name: str) -> MaterialRecord | None:
        """Best-effort match by display name, alias, or key."""


def _slug(text: str) -> str:
    text = text.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "material"


def record_from_openmc_material(
    mat: Any,
    *,
    key: str,
    source: str,
    aliases: Iterable[str] = (),
) -> MaterialRecord:
    """Convert an ``openmc.Material`` into a ``MaterialRecord``."""
    composition: list[CompositionLine] = []
    for entry in getattr(mat, "nuclides", ()) or ():
        composition.append(
            CompositionLine(
                symbol=str(entry.name),
                fraction=float(entry.percent),
                fraction_type=str(entry.percent_type),
            )
        )
    density = getattr(mat, "density", None)
    try:
        density_f = float(density) if density is not None else None
    except (TypeError, ValueError):
        density_f = None
    mid = getattr(mat, "id", None)
    try:
        mid_i = int(mid) if mid is not None else None
    except (TypeError, ValueError):
        mid_i = None
    display = str(getattr(mat, "name", None) or key)
    alias_set = {_slug(a) for a in aliases if a} | {_slug(display), _slug(key)}
    return MaterialRecord(
        key=key,
        display_name=display,
        density=density_f,
        density_units=str(getattr(mat, "density_units", "") or "") or None,
        openmc_id=mid_i,
        composition=tuple(composition),
        source=source,
        aliases=tuple(sorted(a for a in alias_set if a and a != key)),
    )


class OpenMCMaterialsModuleCatalog:
    """Catalog backed by ``src.paratan.materials.material`` module objects."""

    SOURCE = "paratan.materials.material"

    def __init__(self, module: Any | None = None) -> None:
        self._module = module
        self._by_key: dict[str, MaterialRecord] = {}
        self._alias_index: dict[str, str] = {}
        self._loaded = False

    @property
    def source_label(self) -> str:
        return self.SOURCE

    def ensure_loaded(self) -> None:
        if self._loaded:
            return
        module = self._module
        if module is None:
            from src.paratan.materials import material as module  # lazy: openmc + matplotlib
            self._module = module
        records: dict[str, MaterialRecord] = {}
        # Prefer the curated list; fall back to scanning module attributes.
        materials = list(getattr(module, "materials_list", []) or [])
        attr_names: dict[int, str] = {}
        for name, value in vars(module).items():
            if self._is_openmc_material(value):
                attr_names[id(value)] = name
        if not materials:
            materials = [v for v in vars(module).values() if self._is_openmc_material(v)]
        seen: set[int] = set()
        for mat in materials:
            if id(mat) in seen:
                continue
            seen.add(id(mat))
            attr = attr_names.get(id(mat))
            display = str(getattr(mat, "name", None) or attr or "material")
            key = _slug(attr or display)
            # Disambiguate duplicate slugs.
            base, n = key, 2
            while key in records:
                key = f"{base}_{n}"
                n += 1
            aliases = [attr, display] if attr else [display]
            records[key] = record_from_openmc_material(
                mat, key=key, source=self.SOURCE, aliases=aliases
            )
        self._by_key = dict(sorted(records.items(), key=lambda kv: kv[1].display_name.lower()))
        alias_index: dict[str, str] = {}
        for key, rec in self._by_key.items():
            alias_index[key] = key
            for alias in rec.aliases:
                alias_index.setdefault(alias, key)
            alias_index.setdefault(_slug(rec.display_name), key)
        self._alias_index = alias_index
        self._loaded = True

    @staticmethod
    def _is_openmc_material(value: Any) -> bool:
        cls = type(value)
        return cls.__module__.startswith("openmc") and cls.__name__ == "Material"

    def list_records(self) -> list[MaterialRecord]:
        self.ensure_loaded()
        return list(self._by_key.values())

    def get(self, key: str) -> MaterialRecord | None:
        self.ensure_loaded()
        return self._by_key.get(key)

    def resolve(self, name: str) -> MaterialRecord | None:
        self.ensure_loaded()
        if not name:
            return None
        key = self._alias_index.get(_slug(name))
        return self._by_key.get(key) if key else None


def default_catalog() -> MaterialCatalog:
    """Factory used by the viewer; swap here when a new backend lands."""
    return OpenMCMaterialsModuleCatalog()
