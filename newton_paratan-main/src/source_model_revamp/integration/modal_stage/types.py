"""Internal modal stage density, basis cache, and reusable basis state"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.modal.basis import ModalBasisConvergenceAssessment
from source_model_revamp.fbis.modal.density_weighting import Eq42DensityProfile
from source_model_revamp.fbis.modal import ModalFBISBasis, ModalFBISResult

def electron_density_input_is_initial_guess_only(modal_density_closure_model: str) -> bool:
    """Return whether the configured electron density is an initial guess rather than the benchmark prescribed density authority"""
    return str(modal_density_closure_model).strip().lower() != "egedal_beam_plasma_quasineutral"

@dataclass(frozen=True)
class OperatingPointDensityState:
    """
    Solved confined density composition for one operating point
    
    Cell and node profiles are in m⁻³ and include fast D, fast T, total positive charge, electron density, the Eq 70 support mask, and scalar closure diagnostics
    """
    zeta_cells: np.ndarray
    zeta_nodes: np.ndarray
    cell_volumes_m3: np.ndarray
    fast_deuterium_cell_density_m3: np.ndarray
    fast_deuterium_node_density_m3: np.ndarray
    total_positive_charge_cell_density_m3: np.ndarray
    total_positive_charge_node_density_m3: np.ndarray
    electron_cell_density_m3: np.ndarray
    electron_node_density_m3: np.ndarray
    fast_deuterium_midplane_density_m3: float
    electron_midplane_density_m3: float
    fast_deuterium_confined_volume_average_density_m3: float
    electron_confined_volume_average_density_m3: float
    electron_collision_density_m3: float
    electron_parent_maxwellian_n0_m3: float
    profile_support_mask: np.ndarray
    eq70_residual_support_mask: np.ndarray
    profile_scope: str
    electron_profile_scope: str
    expander_electron_profile_available: bool
    scalar_convergence_history: tuple[dict[str, object], ...]
    exact_midplane_quasineutrality_relative_error: float
    confined_inventory_relative_error: float
    fast_tritium_cell_density_m3: np.ndarray | None = None
    fast_tritium_node_density_m3: np.ndarray | None = None
    fast_tritium_midplane_density_m3: float = 0.0
    fast_tritium_confined_volume_average_density_m3: float = 0.0

    def __post_init__(self) -> None:
        """Fill absent fast T profiles with exact zero arrays matching the fast D cell and node grids"""
        if self.fast_tritium_cell_density_m3 is None:
            object.__setattr__(self, "fast_tritium_cell_density_m3", np.zeros_like(np.asarray(self.fast_deuterium_cell_density_m3, dtype=float)))
        if self.fast_tritium_node_density_m3 is None:
            object.__setattr__(self, "fast_tritium_node_density_m3", np.zeros_like(np.asarray(self.fast_deuterium_node_density_m3, dtype=float)))

    def as_metadata(self) -> dict[str, object]:
        """Return density profiles, scalar density values, support masks, and quasineutral inventory diagnostics as metadata"""
        return {'fast_D_midplane_density_m3': self.fast_deuterium_midplane_density_m3, 'fast_D_confined_volume_average_density_m3': self.fast_deuterium_confined_volume_average_density_m3, 'fast_T_midplane_density_m3': self.fast_tritium_midplane_density_m3, 'fast_T_confined_volume_average_density_m3': self.fast_tritium_confined_volume_average_density_m3, 'fast_D_density_profile_m3': self.fast_deuterium_cell_density_m3, 'fast_T_density_profile_m3': self.fast_tritium_cell_density_m3, 'electron_midplane_density_m3': self.electron_midplane_density_m3, 'electron_confined_volume_average_density_m3': self.electron_confined_volume_average_density_m3, 'electron_collision_density_m3': self.electron_collision_density_m3, 'electron_parent_maxwellian_n0_m3': self.electron_parent_maxwellian_n0_m3, 'electron_density_profile_m3': self.electron_cell_density_m3, 'total_positive_charge_profile_m3': self.total_positive_charge_cell_density_m3, 'electron_density_profile_support_mask': self.profile_support_mask, 'eq70_residual_support_mask': self.eq70_residual_support_mask, 'electron_profile_scope': self.electron_profile_scope, 'expander_electron_profile_available': self.expander_electron_profile_available, 'exact_midplane_quasineutrality_relative_error': self.exact_midplane_quasineutrality_relative_error, 'electron_ion_confined_inventory_relative_error': self.confined_inventory_relative_error, 'composition_closure_history': list(self.scalar_convergence_history)}

@dataclass(frozen=True)
class _ModalDensityClosureResult:
    """Single species modal result together with its collision state, solved density closure state, active Eq 42 basis, and convergence history"""
    modal: ModalFBISResult
    collision_state: FBISCollisionParameterState
    electron_midplane_density_m3: float
    electron_collision_density_m3: float
    ion_densities_m3: np.ndarray
    ion_charge_numbers: np.ndarray
    ion_masses_kg: np.ndarray
    converged: bool
    iterations: int
    relative_error: float
    history: list[dict[str, object]]
    density_state: OperatingPointDensityState | None = None
    modal_basis: "_ReusableModalBasisBundle | None" = None
    eq42_density_profile: Eq42DensityProfile | None = None
    eq42_density_basis_converged: bool = False
    eq42_density_basis_iterations: int = 0
    eq42_density_basis_history: list[dict[str, object]] | None = None

@dataclass(frozen=True)
class _ModalBasisCacheEntry:
    """One immutable cached basis result with optional grid convergence assessment and reference first eigenvalue"""
    basis: ModalFBISBasis
    convergence_assessment: ModalBasisConvergenceAssessment | None
    published_reference_lambda1: float | None

@dataclass
class ModalBasisCache:
    """Bounded run local cache for density dependent modal basis results"""
    max_entries: int = 256
    entries: dict[tuple[object, ...], _ModalBasisCacheEntry] = field(default_factory=dict)
    hits: int = 0
    misses: int = 0
    bypasses: int = 0
    evictions: int = 0

    def lookup(self, key: tuple[object, ...]) -> _ModalBasisCacheEntry | None:
        """Return a cached basis entry and update hit or miss counters"""
        entry = self.entries.get(key)
        if entry is None:
            self.misses += 1
        else:
            self.hits += 1
        return entry

    def store(self, key: tuple[object, ...], entry: _ModalBasisCacheEntry) -> None:
        """Store one basis entry and evict the oldest retained entry when the cache is full"""
        if key not in self.entries and len(self.entries) >= self.max_entries:
            self.entries.pop(next(iter(self.entries)))
            self.evictions += 1
        self.entries[key] = entry

    def record_bypass(self) -> None:
        """Record one intentional cache bypass"""
        self.bypasses += 1

    def as_metadata(self) -> dict[str, int]:
        """Return cache hit, miss, bypass, eviction, size, and capacity counters"""
        return {
            "modal_basis_cache_hits": int(self.hits),
            "modal_basis_cache_misses": int(self.misses),
            "modal_basis_cache_bypasses": int(self.bypasses),
            "modal_basis_cache_evictions": int(self.evictions),
            "modal_basis_cache_entries": int(len(self.entries)),
            "modal_basis_cache_max_entries": int(self.max_entries),
        }

@dataclass(frozen=True)
class _ReusableModalBasisBundle:
    """Selected modal basis plus Eq 42 profile, convergence assessment, reference first eigenvalue, and run local cache"""

    basis: ModalFBISBasis
    convergence_assessment: ModalBasisConvergenceAssessment | None
    published_reference_lambda1: float | None = None
    eq42_density_profile: Eq42DensityProfile | None = None
    basis_cache: ModalBasisCache | None = None

    def __getattr__(self, name: str):
        """Forward unknown attributes to the contained ModalFBISBasis for compatibility with earlier callers"""
        return getattr(self.basis, name)
