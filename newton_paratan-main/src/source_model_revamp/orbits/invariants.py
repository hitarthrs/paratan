"""
Single particle invariants and magnetic orbit accessibility
    no electrostatic potential included here

Main variables
    E_J: particle KE in joules
    mass_kg: particle mass in kg
    v_m_s: particle speed in m/s
    B_tilde: normalized magnetic field B / B0
    mirror_ratio: R_M = B_M / B0
    mu_J_per_T: magnetic moment μ = E_perp / B
    Λ: Egedal pitch invariant Λ = μ * B0 / E

Physics equations
    E = 0.5 * m * v^2
    μ = m * v_perp^2 / (2B)
    v_perp^2 / v^2 = Λ * B_tilde
    v_parallel^2 / v^2 =  1 - Λ * B_tilde
    local accessibility
        Λ * B_tilde <= 1
    magnetic trapped/passing boundary
        Λ_M = 1 / R_M
    turning point field
        B_tilde(z_b) = 1 / Λ
"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike

def speed_from_energy_m_s(energy_J: ArrayLike, mass_kg: float):
    """v = sqrt(2E / m)"""
    return np.sqrt(2.0 * np.asarray(energy_J, dtype=float) / mass_kg)

def kinetic_energy_J(speed_m_s: ArrayLike, mass_kg: float):
    """E = 0.5 * m * v^2"""
    v = np.asarray(speed_m_s, dtype=float)
    
    return 0.5 * mass_kg * v**2

def magnetic_moment_J_per_T(perpendicular_speed_m_s: ArrayLike, local_B_T: ArrayLike, mass_kg: float):
    """μ = m * v_perp^2 / (2B)"""
    v_perp = np.asarray(perpendicular_speed_m_s, dtype=float)
    B = np.asarray(local_B_T, dtype=float)

    return 0.5 * mass_kg * v_perp**2 / B

def lambda_from_mu(magnetic_moment_J_per_T: ArrayLike, B0_T: float, energy_J: ArrayLike):
    """Λ = μ * B0 / E"""
    mu = np.asarray(magnetic_moment_J_per_T, dtype=float)
    E = np.asarray(energy_J, dtype=float)

    return mu * B0_T / E

def lambda_from_pitch_angle(pitch_angle_rad: ArrayLike, B_tilde: ArrayLike = 1.0):
    """
    Λ = sin^2(θ) / B_tilde

    At the midplane, B_tilde = 1, so
        Λ = sin^2(θ_midplane)
    """
    θ = np.asarray(pitch_angle_rad, dtype=float)
    B = np.asarray(B_tilde, dtype=float)

    return np.sin(θ) ** 2 / B

def lambda_from_xi(xi: ArrayLike, B_tilde: ArrayLike = 1.0):
    """
    ξ = v_parallel / v
    ξ^2 = 1 - Λ B_tilde
    So
        Λ = (1 - ξ^2) / B_tilde
    """
    xi = np.asarray(xi, dtype=float)
    B = np.asarray(B_tilde, dtype=float)

    return (1.0 - xi**2) / B

def xi_squared(Lambda: ArrayLike, B_tilde: ArrayLike):
    """ξ^2 = v_parallel^2 / v^2 = 1 - Λ * B_tilde"""
    Lam = np.asarray(Lambda, dtype=float)
    B = np.asarray(B_tilde, dtype=float)

    return 1.0 - Lam * B

def perpendicular_velocity_fraction_squared(Lambda: ArrayLike, B_tilde: ArrayLike):
    """v_perp^2 / v^2 = Λ * B_tilde"""
    Lam = np.asarray(Lambda, dtype=float)
    B = np.asarray(B_tilde, dtype=float)

    return Lam * B

def parallel_velocity_fraction_squared(Lambda: ArrayLike, B_tilde: ArrayLike):
    """v_parallel^2 / v^2 = 1 - Λ * B_tilde"""
    return xi_squared(Lambda=Lambda, B_tilde=B_tilde)

def perpendicular_energy_J(energy_J: ArrayLike, Lambda: ArrayLike, B_tilde: ArrayLike):
    """E_perp = E * Λ * B_tilde"""
    E = np.asarray(energy_J, dtype=float)
    Lam = np.asarray(Lambda, dtype=float)
    B = np.asarray(B_tilde, dtype=float)

    return E * Lam * B

def parallel_energy_J(energy_J: ArrayLike, Lambda: ArrayLike, B_tilde: ArrayLike):
    """E_parallel = E (1 - Λ B_tilde)"""
    E = np.asarray(energy_J, dtype=float)

    return E * parallel_velocity_fraction_squared(Lambda=Lambda, B_tilde=B_tilde)

def perpendicular_speed_m_s(speed_m_s: ArrayLike, Lambda: ArrayLike, B_tilde: ArrayLike):
    """v_perp = v * sqrt(Λ * B_tilde)"""
    v = np.asarray(speed_m_s, dtype=float)

    return v * np.sqrt(perpendicular_velocity_fraction_squared(Lambda=Lambda, B_tilde=B_tilde,))

def parallel_speed_m_s(speed_m_s: ArrayLike, Lambda: ArrayLike, B_tilde: ArrayLike):
    """|v_parallel| = v * sqrt(1 - Λ * B_tilde)"""
    v = np.asarray(speed_m_s, dtype=float)

    return v * np.sqrt(parallel_velocity_fraction_squared(Lambda=Lambda, B_tilde=B_tilde))

def signed_parallel_speed_m_s(speed_m_s: ArrayLike, Lambda: ArrayLike, B_tilde: ArrayLike, direction_sign: ArrayLike):
    """v_parallel = s * v * sqrt(1 - Λ * B_tilde)"""
    sign = np.asarray(direction_sign, dtype=float)

    return sign * parallel_speed_m_s(speed_m_s=speed_m_s, Lambda=Lambda,B_tilde=B_tilde)

def lambda_max_accessible( B_tilde: ArrayLike):
    """
    Local accessibility condition
        Λ B_tilde <= 1
    so
        Λ_max(z) = 1 / B_tilde(z)
    """
    B = np.asarray(B_tilde, dtype=float)
    return 1.0 / B

def is_accessible(Lambda: ArrayLike, B_tilde: ArrayLike):
    """
    Magnetic only local accessibility
        Λ * B_tilde <= 1
    """
    Lam = np.asarray(Lambda, dtype=float)
    B = np.asarray(B_tilde, dtype=float)
    return Lam * B <= 1.0

def magnetic_loss_boundary_lambda(mirror_ratio: float) -> float:
    """
    Magnetic trapped/passing boundary
        Λ_M = 1 / R_M
    """
    return 1.0 / mirror_ratio

def is_magnetically_trapped(Lambda: ArrayLike, mirror_ratio: float):
    """
    Magnetic trapped condition
        Λ > Λ_M = 1 / R_M
    """
    Lam = np.asarray(Lambda, dtype=float)
    return Lam > magnetic_loss_boundary_lambda(mirror_ratio)

def is_magnetically_passing(Lambda: ArrayLike, mirror_ratio: float):
    """
    Magnetic passing/loss-cone condition
        Λ < Λ_M = 1 / R_M
    """
    Lam = np.asarray(Lambda, dtype=float)
    return Lam < magnetic_loss_boundary_lambda(mirror_ratio)

def is_magnetic_boundary(Lambda: ArrayLike, mirror_ratio: float):
    """
    Magnetic trapped/passing boundary condition
        Λ = Λ_M = 1 / R_M
    """
    Lam = np.asarray(Lambda, dtype=float)
    return Lam == magnetic_loss_boundary_lambda(mirror_ratio)

def xi_trapped_passing_boundary(mirror_ratio: float):
    """
    Square-mirror trapped/passing boundary in ξ
        ξ_TP = sqrt(1 - 1 / R_M)
    """
    return np.sqrt(1.0 - 1.0 / mirror_ratio)

def turning_point_B_tilde(Lambda: ArrayLike):
    """
    Magnetic turning-point condition
        B_tilde(z_b) = 1 / Λ
    Gives required B_tilde at the turning point
    """
    Lam = np.asarray(Lambda, dtype=float)

    return 1.0 / Lam

def local_pitch_angle_rad(Lambda: ArrayLike, B_tilde: ArrayLike):
    """sin^2(θ) = Λ * B_tilde"""
    Lam = np.asarray(Lambda, dtype=float)
    B = np.asarray(B_tilde, dtype=float)

    return np.arcsin(np.sqrt(Lam * B))