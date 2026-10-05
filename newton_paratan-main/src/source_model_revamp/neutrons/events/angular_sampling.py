"""
Sampling from evaluated ENDF B VIII.1 center of mass angular distributions

The tabulated path samples a piecewise linear probability density in `mu = cos(theta_cm)`
The Legendre path linearly interpolates angular coefficients in incident energy and samples the resulting normalized law by rejection
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from source_model_revamp.nuclear_data.endf_b_viii1.angular_distribution import EvaluatedFusionAngularDistribution
from source_model_revamp.nuclear_data.endf_b_viii1.interpolation import bracket_linear_energy

@dataclass(frozen=True)
class _PiecewiseLinearKnotSampler:
    """
    Inverse CDF sampler for one tabulated angular knot
    
    `mu`, `probability_density`, and `cumulative_area` are one dimensional arrays spanning `mu = −1` through `mu = 1`
    """
    mu: np.ndarray
    probability_density: np.ndarray
    cumulative_area: np.ndarray

    @classmethod
    def from_arrays(cls, mu: np.ndarray, probability_density: np.ndarray,) -> "_PiecewiseLinearKnotSampler":
        """Build a normalized piecewise linear knot sampler from tabulated `p(mu)` values"""
        cosine = np.asarray(mu, dtype=float)
        probability = np.asarray(probability_density, dtype=float)
        if cosine.ndim != 1 or probability.shape != cosine.shape:
            raise ValueError("tabulated angular data must be equal length 1D arrays")
        if cosine.size < 2 or cosine[0] != -1.0 or cosine[-1] != 1.0:
            raise ValueError("tabulated angular data must span minus one to one")
        widths = np.diff(cosine)
        segment_area = 0.5 * (probability[:-1] + probability[1:]) * widths
        if np.any(segment_area < 0.0) or not np.all(np.isfinite(segment_area)):
            raise ValueError("tabulated angular segment areas must be finite and nonnegative")
        total = float(np.sum(segment_area))
        if total <= 0.0:
            raise ValueError("tabulated angular distribution has zero normalization")
        normalized_probability = probability / total
        normalized_area = segment_area / total
        cumulative = np.concatenate(([0.0], np.cumsum(normalized_area)))
        cumulative[-1] = 1.0

        return cls(mu=cosine, probability_density=normalized_probability, cumulative_area=cumulative)

    def sample(self, random_values: np.ndarray) -> np.ndarray:
        """Sample `mu` values by analytically inverting the integrated linear density inside each knot interval"""
        values = np.asarray(random_values, dtype=float)
        if np.any(values < 0.0) or np.any(values >= 1.0):
            raise ValueError("random_values must lie inside [0, 1)")
        segment = np.searchsorted(self.cumulative_area, values, side="right") - 1
        segment = np.clip(segment, 0, self.mu.size - 2)
        local_area = values - self.cumulative_area[segment]
        x0 = self.mu[segment]
        width = self.mu[segment + 1] - x0
        y0 = self.probability_density[segment]
        y1 = self.probability_density[segment + 1]
        slope = (y1 - y0) / width
        nearly_constant = np.abs(slope) <= 32.0 * np.finfo(float).eps * np.maximum(np.abs(y0) / np.maximum(width, np.finfo(float).tiny), 1.0)
        displacement = np.zeros_like(values)
        if np.any(nearly_constant):
            constant_density = np.maximum(y0[nearly_constant], np.finfo(float).tiny)
            displacement[nearly_constant] = local_area[nearly_constant] / constant_density
        varying = ~nearly_constant
        if np.any(varying):
            discriminant = np.maximum(y0[varying] ** 2 + 2.0 * slope[varying] * local_area[varying], 0.0)
            denominator = y0[varying] + np.sqrt(discriminant)
            fallback = np.abs(denominator) <= np.finfo(float).tiny
            value = np.empty(np.count_nonzero(varying), dtype=float)
            value[~fallback] = 2.0 * local_area[varying][~fallback] / denominator[~fallback]
            if np.any(fallback):
                value[fallback] = (-y0[varying][fallback] + np.sqrt(discriminant[fallback])) / slope[varying][fallback]
            displacement[varying] = value

        return np.clip(x0 + displacement, -1.0, 1.0)

class EvaluatedFusionAngularSampler:
    """
    Random sampler for one normalized evaluated fusion angular law
    
    Representation code `12` uses tabulated probability density in `mu`
    Representation code `0` uses Legendre coefficients
    """
    def __init__(self, distribution: EvaluatedFusionAngularDistribution):
        """Prepare the representation specific samplers for one evaluated angular distribution"""
        self.distribution = distribution
        self._tabulated_samplers: tuple[_PiecewiseLinearKnotSampler, ...] = ()
        if distribution.representation_code == 12:
            samplers = []
            for index in range(distribution.incident_energy_eV.size):
                mu, probability = distribution._tabulated_knot(index)
                samplers.append(_PiecewiseLinearKnotSampler.from_arrays(mu, probability))
            self._tabulated_samplers = tuple(samplers)

    def _sample_tabulated(self, incident_energy_eV: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """
        Sample a tabulated angular law at requested incident energies
        
        Linear energy interpolation is represented as a stochastic mixture of the lower and upper energy knot distributions
        """
        lower, upper, fraction = bracket_linear_energy(incident_energy_eV, self.distribution.incident_energy_eV,)
        # A lower or upper knot draw reproduces the linearly interpolated angular density in expectation
        choose_upper = rng.random(incident_energy_eV.size) < fraction
        chosen = np.where(choose_upper, upper, lower)
        random_values = rng.random(incident_energy_eV.size)
        result = np.empty(incident_energy_eV.size, dtype=float)
        for knot_index in np.unique(chosen):
            selected = chosen == knot_index
            result[selected] = self._tabulated_samplers[int(knot_index)].sample(random_values[selected])

        return result

    @staticmethod
    def _legendre_probability_rows(coefficients: np.ndarray, mu: np.ndarray,) -> np.ndarray:
        """Evaluate row specific Legendre angular probability densities at supplied `mu` values"""
        count, order_count = coefficients.shape
        result = np.full(count, 0.5, dtype=float)
        if order_count == 0:
            return result
        p_previous = np.ones(count, dtype=float)
        p_current = mu.copy()
        for order in range(1, order_count + 1):
            if order == 1:
                polynomial = p_current
            else:
                polynomial = ((2.0 * order - 1.0) * mu * p_current - (order - 1.0) * p_previous) / order
                p_previous, p_current = p_current, polynomial
            result += (0.5 * (2.0 * order + 1.0) * coefficients[:, order - 1] * polynomial)

        return result

    def _sample_legendre(self, incident_energy_eV: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """
        Sample interpolated Legendre angular laws by rejection on `mu` in `[−1, 1]`
        
        The envelope uses `|P_l(mu)| <= 1` to bound each interpolated row
        """
        if self.distribution.legendre_coefficients is None:
            raise ValueError("Legendre coefficients are unavailable")
        lower, upper, fraction = bracket_linear_energy(incident_energy_eV, self.distribution.incident_energy_eV,)
        coefficients = (self.distribution.legendre_coefficients[lower] + fraction[:, None] * (self.distribution.legendre_coefficients[upper] - self.distribution.legendre_coefficients[lower]))
        orders = np.arange(1, coefficients.shape[1] + 1, dtype=float)
        # The Legendre bound |P_l| <= 1 gives a valid row specific rejection envelope
        envelope = 0.5 + np.sum(0.5 * (2.0 * orders + 1.0)[None, :] * np.abs(coefficients), axis=1,)
        if np.any(~np.isfinite(envelope)) or np.any(envelope <= 0.0):
            raise ValueError("Legendre rejection envelope is invalid")
        result = np.empty(incident_energy_eV.size, dtype=float)
        pending = np.ones(incident_energy_eV.size, dtype=bool)
        attempts = 0
        while np.any(pending):
            attempts += 1
            if attempts > 10000:
                raise RuntimeError("Legendre angular rejection sampling did not converge")
            indices = np.flatnonzero(pending)
            candidates = rng.uniform(-1.0, 1.0, indices.size)
            probability = self._legendre_probability_rows(coefficients[indices], candidates)
            if np.any(probability < -1.0e-12):
                raise ValueError("interpolated ENDF angular probability is materially negative")
            probability = np.maximum(probability, 0.0)
            accepted = rng.random(indices.size) * envelope[indices] <= probability
            if np.any(accepted):
                accepted_indices = indices[accepted]
                result[accepted_indices] = candidates[accepted]
                pending[accepted_indices] = False

        return result

    def sample_mu(self, incident_energy_eV: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """Sample center of mass emission cosine for incident energies inside the evaluated energy range"""
        energy = np.asarray(incident_energy_eV, dtype=float)
        if energy.ndim != 1 or np.any(~np.isfinite(energy)):
            raise ValueError("incident_energy_eV must be a finite 1D array")
        lower, upper = self.distribution.incident_energy_bounds_eV
        if np.any(energy < lower) or np.any(energy > upper):
            raise ValueError(f"incident energy must lie inside [{lower}, {upper}] eV")
        if self.distribution.representation_code == 12:
            return self._sample_tabulated(energy, rng)
        
        return self._sample_legendre(energy, rng)