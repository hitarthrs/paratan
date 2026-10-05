"""Shared physical constants and unit conversions for the source model

Read physical constants from SciPy and expose them under shared names
Conversion factors multiply values in the named source unit
"""
from __future__ import annotations
from scipy.constants import c, elementary_charge, epsilon_0, electron_mass, physical_constants, pi, mu_0

SPEED_OF_LIGHT_M_S = float(c)

# Positive elementary charge magnitude in coulombs
ELECTRON_CHARGE_C = float(elementary_charge)
VACUUM_PERMITTIVITY_F_PER_M = float(epsilon_0)

ELECTRON_MASS_KG = float(electron_mass)

# Nuclear and neutron masses in kilograms
# The table index selects the value from the value, unit, uncertainty tuple
PROTON_MASS_KG = float(physical_constants["proton mass"][0])
DEUTERON_MASS_KG = float(physical_constants["deuteron mass"][0])
TRITON_MASS_KG = float(physical_constants["triton mass"][0])
NEUTRON_MASS_KG = float(physical_constants["neutron mass"][0])

# Helion and alpha denote the helium 3 and helium 4 nuclei
HELION_MASS_KG = float(physical_constants["helion mass"][0])
ALPHA_PARTICLE_MASS_KG = float(physical_constants["alpha particle mass"][0])

# Reduced Planck constant in joule seconds
HBAR_J_S = float(physical_constants["Planck constant over 2 pi"][0])

# Convert electronvolt energy units to joules
EV_TO_J = ELECTRON_CHARGE_C
KEV_TO_J = 1.0e3 * EV_TO_J
MEV_TO_J = 1.0e6 * EV_TO_J

# Convert joules to electronvolt energy units
J_TO_EV = 1.0 / EV_TO_J
J_TO_KEV = 1.0 / KEV_TO_J
J_TO_MEV = 1.0 / MEV_TO_J

# Cross section area conversions
BARN_TO_M2 = 1.0e-28
M2_TO_BARN = 1.0e28

# Volume conversions
# Number density conversions use the reciprocal volume factors
CM3_TO_M3 = 1.0e-6
M3_TO_CM3 = 1.0e6

# Length conversion from centimeters to meters
CM_TO_M = 0.01

# Expose vacuum permeability and pi under their existing names
mu_0 = mu_0
pi=pi