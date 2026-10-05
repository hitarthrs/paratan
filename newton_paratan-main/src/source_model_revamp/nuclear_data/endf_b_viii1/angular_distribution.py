"""Load and evaluate conditional DD and DT neutron angular laws from ENDF B VIII.1"""
from __future__ import annotations
from dataclasses import dataclass, field, replace
from importlib import resources
from typing import Any
import numpy as np
from numpy.polynomial.legendre import legval
from numpy.typing import ArrayLike
from source_model_revamp.nuclear_data.endf_b_viii1.interpolation import bracket_linear_energy, linear_blend, linear_probability_density, validate_linear_energy_interpolation

_PACKAGE = "source_model_revamp.nuclear_data.endf_b_viii1"
_REACTION_FILES = {"dd_n": "ddn_angular_data.npz", "dt_n": "dtn_angular_data.npz"}
_REACTION_ALIASES = {"dd": "dd_n", "dd_n": "dd_n", "d(d,n)3he": "dd_n", "dt": "dt_n", "dt_n": "dt_n", "t(d,n)4he": "dt_n"}

@dataclass(frozen=True)
class EvaluatedFusionAngularDistribution:
    """Immutable evaluated angular law for one neutron producing fusion reaction

    `incident_energy_eV` has shape `(n_energy,)` and uses equivalent deuteron
    laboratory kinetic energy on a stationary target

    `mu = cos(theta_cm)` is measured from the incident deuteron direction in the
    center of mass frame and `p(mu | E)` is normalized so that
    `integral from −1 to 1 of p(mu | E) dmu = 1`

    Representation code `12` stores tabulated `p(mu | E)` values while code `0`
    stores Legendre coefficients for the same conditional probability density
    """
    reaction_key: str
    incident_energy_eV: np.ndarray
    energy_interpolation_breakpoints: np.ndarray
    energy_interpolation_laws: np.ndarray
    representation_code: int
    mat: int
    mf: int
    mt: int
    lct: int
    law: int
    lang: int
    q_value_eV: float
    target_awr: float
    projectile_awr: float
    neutron_awp: float
    residual_awp: float
    residual_zap: int
    angular_offsets: np.ndarray | None
    angular_mu: np.ndarray | None
    angular_probability_density_per_mu: np.ndarray | None
    angular_item_count: np.ndarray | None
    legendre_coefficients: np.ndarray | None
    legendre_order_count: np.ndarray | None
    runtime_validation: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate energy interpolation metadata, reaction identity, and representation fields"""
        energy = np.asarray(self.incident_energy_eV, dtype=float)
        if energy.ndim != 1 or energy.size < 2 or np.any(np.diff(energy) <= 0.0):
            raise ValueError("incident_energy_eV must be a strictly increasing 1D grid")
        validate_linear_energy_interpolation(self.energy_interpolation_breakpoints, self.energy_interpolation_laws, energy.size)
        expected_identity = { "dd_n": (128, 6, 50, 2, 2, 12, 2003), "dt_n": (131, 6, 50, 2, 2, 0, 2004)}
        if self.reaction_key not in expected_identity:
            raise ValueError("fusion angular data have an unsupported reaction identity")
        identity = (self.mat, self.mf, self.mt, self.lct, self.law, self.lang, self.residual_zap)
        if identity != expected_identity[self.reaction_key]:
            raise ValueError("fusion angular dataset reaction identity is inconsistent")
        if self.representation_code == 12:
            if any(value is None for value in (self.angular_offsets, self.angular_mu, self.angular_probability_density_per_mu, self.angular_item_count)):
                raise ValueError("LANG=12 data are incomplete")
        elif self.representation_code == 0:
            if self.legendre_coefficients is None or self.legendre_order_count is None:
                raise ValueError("LANG=0 data are incomplete")
        else:
            raise ValueError(f"unsupported LAW=2 representation {self.representation_code}")

    @property
    def incident_energy_bounds_eV(self) -> tuple[float, float]:
        """Return the inclusive equivalent deuteron laboratory energy bounds in eV"""
        return float(self.incident_energy_eV[0]), float(self.incident_energy_eV[-1])

    @property
    def representation(self) -> str:
        """Return the canonical angular representation name"""
        return ("tabulated_probability_linear_in_mu" if self.representation_code == 12 else "legendre_coefficients")

    def _tabulated_knot(self, index: int) -> tuple[np.ndarray, np.ndarray]:
        """Return the tabulated `mu` grid and `p(mu)` values for one incident energy knot"""
        if self.angular_offsets is None or self.angular_mu is None:
            raise ValueError("tabulated angular data are unavailable")
        if self.angular_probability_density_per_mu is None:
            raise ValueError("tabulated angular probability data are unavailable")
        start = int(self.angular_offsets[index])
        stop = int(self.angular_offsets[index + 1])
        return (self.angular_mu[start:stop], self.angular_probability_density_per_mu[start:stop],)

    def _tabulated_probability_at_knot(self, index: int, mu: np.ndarray) -> np.ndarray:
        """Evaluate one tabulated angular knot by linear interpolation in `mu`"""
        knot_mu, knot_probability = self._tabulated_knot(index)
        return linear_probability_density(mu, knot_mu, knot_probability)

    def _legendre_probability_at_coefficients(self, coefficients: np.ndarray, mu: np.ndarray) -> np.ndarray:
        """Evaluate `p(mu) = 1/2 + sum_l [(2l + 1)/2] a_l P_l(mu)`"""
        # ENDF coefficients `a_l` omit the normalized `l = 0` probability term
        polynomial_coefficients = np.zeros(coefficients.size + 1, dtype=float)
        polynomial_coefficients[0] = 0.5
        orders = np.arange(1, coefficients.size + 1, dtype=float)
        polynomial_coefficients[1:] = 0.5 * (2.0 * orders + 1.0) * coefficients

        return np.asarray(legval(mu, polynomial_coefficients), dtype=float)

    def legendre_coefficients_at_energy(self, incident_energy_eV: float) -> np.ndarray:
        """Interpolate the padded Legendre coefficient row at one incident energy

        The returned shape is `(n_order_max,)` and zero padded orders remain zero
        """
        if self.legendre_coefficients is None:
            raise ValueError("this angular distribution is not represented by Legendre coefficients")
        energy = float(incident_energy_eV)
        lower, upper, fraction = bracket_linear_energy(np.asarray(energy), self.incident_energy_eV)
        lower_index = int(np.asarray(lower))
        upper_index = int(np.asarray(upper))
        weight = float(np.asarray(fraction))

        return linear_blend(self.legendre_coefficients[lower_index], self.legendre_coefficients[upper_index], weight)

    def probability_density_mu(self, incident_energy_eV: float, mu: ArrayLike) -> np.ndarray | float:
        """Evaluate the normalized conditional density `p(mu | E)` per unit `mu`

        Scalar `mu` input returns a float and array input preserves the query shape
        """
        query = np.asarray(mu, dtype=float)
        scalar = query.ndim == 0
        if np.any(~np.isfinite(query)) or np.any(query < -1.0) or np.any(query > 1.0):
            raise ValueError("mu must be finite and lie inside [-1, 1]")
        energy = float(incident_energy_eV)
        lower, upper, fraction = bracket_linear_energy(np.asarray(energy), self.incident_energy_eV)
        lower_index = int(np.asarray(lower))
        upper_index = int(np.asarray(upper))
        weight = float(np.asarray(fraction))
        if self.representation_code == 12:
            lower_probability = self._tabulated_probability_at_knot(lower_index, query)
            upper_probability = self._tabulated_probability_at_knot(upper_index, query)
            probability = linear_blend(lower_probability, upper_probability, weight)
        else:
            coefficients = self.legendre_coefficients_at_energy(energy)
            probability = self._legendre_probability_at_coefficients(coefficients, query)
        if np.any(probability < -1.0e-12):
            raise ValueError("interpolated ENDF angular probability is materially negative")
        probability = np.maximum(probability, 0.0)
        if scalar:
            return float(probability)
        
        return probability

    def probability_density_solid_angle_sr(self, incident_energy_eV: float, mu: ArrayLike) -> np.ndarray | float:
        """Return `p_Omega = p_mu / (2 pi)` in `sr^−1` for uniform center of mass azimuth"""
        probability = self.probability_density_mu(incident_energy_eV, mu)
        if np.ndim(probability) == 0:
            return float(probability) / (2.0 * np.pi)
        
        return np.asarray(probability, dtype=float) / (2.0 * np.pi)

    def normalization_integral(self, incident_energy_eV: float) -> float:
        """Integrate `p(mu | E)` over `−1 <= mu <= 1` for a normalization check"""
        energy = float(incident_energy_eV)
        if self.representation_code == 12:
            lower, upper, _ = bracket_linear_energy(np.asarray(energy), self.incident_energy_eV)
            lower_mu, _ = self._tabulated_knot(int(np.asarray(lower)))
            upper_mu, _ = self._tabulated_knot(int(np.asarray(upper)))
            mu = np.unique(np.concatenate((lower_mu, upper_mu)))
            probability = np.asarray(self.probability_density_mu(energy, mu), dtype=float)
            return float(np.trapezoid(probability, mu))
        nodes, weights = np.polynomial.legendre.leggauss(96)
        probability = np.asarray(self.probability_density_mu(energy, nodes), dtype=float)

        return float(np.dot(weights, probability))

def _angular_probability_validation(distribution: EvaluatedFusionAngularDistribution) -> dict[str, Any]:
    """Validate knot normalization and nonnegative angular probability over the evaluated law"""
    normalization_errors: list[float] = []
    minimum_probability = np.inf
    if distribution.representation_code == 12:
        for index in range(distribution.incident_energy_eV.size):
            mu, probability = distribution._tabulated_knot(index)
            if (mu.size < 2 or mu[0] != -1.0 or mu[-1] != 1.0 or np.any(np.diff(mu) <= 0.0) or np.any(~np.isfinite(probability)) or np.any(probability < 0.0)):
                raise ValueError("tabulated ENDF angular PDF is invalid")
            integral = float(np.trapezoid(probability, mu))
            normalization_errors.append(abs(integral - 1.0))
            minimum_probability = min(minimum_probability, float(np.min(probability)))
    else:
        coefficients = np.asarray(distribution.legendre_coefficients, dtype=float)
        order_counts = np.asarray(distribution.legendre_order_count, dtype=np.int64)
        if (coefficients.ndim != 2 or order_counts.shape != (coefficients.shape[0],) or coefficients.shape[0] != distribution.incident_energy_eV.size or np.any(~np.isfinite(coefficients)) or np.any(order_counts < 0) or np.any(order_counts > coefficients.shape[1])):
            raise ValueError("Legendre ENDF angular coefficients are invalid")
        for index, order_count in enumerate(order_counts):
            count = int(order_count)
            probability_coefficients = np.zeros(count + 1, dtype=float)
            probability_coefficients[0] = 0.5
            if count:
                orders = np.arange(1, count + 1, dtype=float)
                probability_coefficients[1:] = (0.5 * (2.0 * orders + 1.0) * coefficients[index, :count])
            polynomial = np.polynomial.legendre.Legendre(probability_coefficients)
            derivative_roots = polynomial.deriv().roots()
            real_roots = np.real(derivative_roots[np.isreal(derivative_roots)])
            candidates = np.concatenate((np.asarray([-1.0, 1.0]), real_roots[(real_roots > -1.0) & (real_roots < 1.0)]))
            values = np.asarray(polynomial(candidates), dtype=float)
            minimum_probability = min(minimum_probability, float(np.min(values)))
            normalization_errors.append(abs(float(polynomial.integ()(1.0) - polynomial.integ()(-1.0)) - 1.0))

    maximum_normalization_error = max(normalization_errors, default=np.inf)
    if maximum_normalization_error > 1.0e-5:
        raise ValueError("ENDF angular PDF normalization is invalid")
    if not np.isfinite(minimum_probability) or minimum_probability < -1.0e-10:
        raise ValueError("ENDF angular PDF is materially negative")
    
    return {
        "angular_law_validation_passed": True,
        "angular_normalization_measured": True,
        "maximum_angular_normalization_error": maximum_normalization_error,
        "angular_normalization_tolerance": 1.0e-5,
        "minimum_angular_probability_density_per_mu": minimum_probability,
        "minimum_angular_probability_tolerance": -1.0e-10,
        "incident_energy_min_eV": distribution.incident_energy_bounds_eV[0],
        "incident_energy_max_eV": distribution.incident_energy_bounds_eV[1],
        "incident_energy_point_count": int(distribution.incident_energy_eV.size),
        "reference_frame": "center_of_mass",
        "incident_energy_frame": "deuteron_lab_on_stationary_target",
    }

def _normalize_reaction_key(reaction_key: str) -> str:
    """Map accepted DD and DT aliases to the canonical neutron reaction key"""
    normalized = str(reaction_key).strip().lower()
    try:
        return _REACTION_ALIASES[normalized]
    except KeyError as exc:
        raise ValueError(f"unsupported fusion angular reaction {reaction_key!r}") from exc

def load_evaluated_fusion_angular_distribution(reaction_key: str) -> EvaluatedFusionAngularDistribution:
    """Load one packaged angular dataset and attach measured runtime validation results"""
    key = _normalize_reaction_key(reaction_key)
    resource = resources.files(_PACKAGE).joinpath(_REACTION_FILES[key])
    with resources.as_file(resource) as path:
        with np.load(path, allow_pickle=False) as archive:
            arrays = {name: np.asarray(archive[name]).copy() for name in archive.files}
    representation_code = int(arrays["representation_code"].item())
    distribution = EvaluatedFusionAngularDistribution(
        reaction_key=key,
        incident_energy_eV=arrays["incident_energy_eV"],
        energy_interpolation_breakpoints=arrays["energy_interpolation_breakpoints"],
        energy_interpolation_laws=arrays["energy_interpolation_laws"],
        representation_code=representation_code,
        mat=int(arrays["mat"].item()),
        mf=int(arrays["mf"].item()),
        mt=int(arrays["mt"].item()),
        lct=int(arrays["lct"].item()),
        law=int(arrays["law"].item()),
        lang=int(arrays["lang"].item()),
        q_value_eV=float(arrays["q_value_eV"].item()),
        target_awr=float(arrays["target_awr"].item()),
        projectile_awr=float(arrays["projectile_awr"].item()),
        neutron_awp=float(arrays["neutron_awp"].item()),
        residual_awp=float(arrays["residual_awp"].item()),
        residual_zap=int(arrays["residual_zap"].item()),
        angular_offsets=arrays.get("angular_offsets"),
        angular_mu=arrays.get("angular_mu"),
        angular_probability_density_per_mu=arrays.get("angular_probability_density_per_mu"),
        angular_item_count=arrays.get("angular_item_count"),
        legendre_coefficients=arrays.get("legendre_coefficients"),
        legendre_order_count=arrays.get("legendre_order_count"),
    )
    angular_validation = _angular_probability_validation(distribution)

    return replace(distribution, runtime_validation={**angular_validation,},)

def load_dd_neutron_angular_distribution() -> EvaluatedFusionAngularDistribution:
    """Load the evaluated `D(d,n)3He` neutron angular distribution"""
    return load_evaluated_fusion_angular_distribution("dd_n")

def load_dt_neutron_angular_distribution() -> EvaluatedFusionAngularDistribution:
    """Load the evaluated `T(d,n)4He` neutron angular distribution"""
    return load_evaluated_fusion_angular_distribution("dt_n")