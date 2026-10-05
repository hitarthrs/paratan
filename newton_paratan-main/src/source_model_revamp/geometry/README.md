# Geometry 

This module supplies magnetic fields, axial domains, plasma facing surfaces, and seed profiles for the source calculation.
Beam, orbit, and expander calculations use these geometric quantities.

## Files

__init__.py:              Exposes shared geometry types and selected construction functions.
axial_grid_geometry.py:   Calculates axial cell sizes, flux tube areas, cell volumes, and volume integrals.
background_profiles.py:   Evaluates and validates D and T seed or reference shapes, normalizes their midplane densities, and calculates inventories.
coil_fields.py:           Evaluates on axis fields and gradients from effective thin circular coils and reports sampled mirror metrics.
device_domains.py:        Converts ParaTAN dimensions into component locations and constructs the named physics domains.
egedal_mirror_profile.py: Evaluates the analytic field from Egedal Eq 56 and its central mirror gradient.
magnetic_topology.py:     Locates independent left and right magnetic throats and checks the central magnetic well.
plasma_boundaries.py:     Builds material interfaces, applies boundary settings, and finds exhaust path restrictions.

## How it fits into the model

integration/config/geometry.py derives effective coils from ParaTAN inputs and fits their strengths.
integration/pipeline_stages/geometry_stage.py combines the field, throat locations, material boundaries, grids, and seed profiles into GeometryStageResult.
mirror_paratan/geometry/ constructs the OpenMC regions. This package supplies the reduced geometry used by the source physics.

## Coordinates and arrays
| Quantity | Convention |

| Machine positions and radii    | Meters, including any configured midplane offset |
| ParaTAN layout inputs          | Centimeters, converted by derive_paratan_layout |
| Analytic coordinate zeta       | Position relative to the midplane divided by mirror half length |
| B_tilde                        | Magnetic field divided by the reference field B0 |
| Axial cell faces and profiles  | N cells have N + 1 faces and N cell values |

The analytic derivative helper is intended for the central interval abs(zeta) <= 1. 
Flux tube areas use A = A0 / B_tilde. 
Cell volumes use that area at the cell center multiplied by the cell width.

## References

The analytic magnetic field follows J Egedal et al, *Nuclear Fusion* 62 (2022) 126053, Eq 56 and Figure 10.