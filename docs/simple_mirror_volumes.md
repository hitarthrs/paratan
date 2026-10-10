# Analytical simple-mirror component volumes

After constructing the port-enabled `SimpleMachineBuilder`, call
`builder.get_component_volumes()` or
`builder.get_component_volumes(cell_ids=[2000, 2001])`.

Records are keyed by OpenMC cell ID and contain `name`, `material_id`, `volume`,
`units='cm3'`, and `method='analytical'`. No geometry or material inventory is
modified. Room air is excluded. HF magnets, casing layers and shields have
separate named left/right cells and per-side volumes. This is a dimension-based API: editing cell regions
manually after building is not supported.

Functions live in `src/paratan/models/simple_mirror_volumes.py`:

- `cylinder_volume(radius, length)`
- `annulus_volume(inner_radius, outer_radius, length)`
- `port_volume(inner_radius, outer_radius, x_limit)`
- `cylinder_intersection_volume(radius, port_radius, angle=pi/4)`
- `crossing_ports_volume(radius)`

All dimensions are centimeters. Production evaluation uses no numerical
quadrature, meshes, or Monte Carlo. Oblique cylinder intersections have the
closed-form expression

    V = 2*pi*a²*b/sin(angle) * 2F1(-1/2, 1/2; 2; a²/b²)

where `a` and `b` are the smaller and larger radii. SciPy's special-function
evaluator supplies `2F1`; this is equivalent to a complete elliptic-integral
expression. The two perpendicular equal-radius port overlap is `16*r³/3`.
Subtract the port union once, not each port independently. For annular layers
whose bores contain that overlap, the overlap cancels between outer and inner
volumes.

Vessel cylinders/frusta and their overlaps with end cells are evaluated with
piecewise polynomial antiderivatives. LF coils that cannot intersect the ports
retain their exact annular/shell volumes.

Supported port cuts are completely traversed cylindrical domains with the
entire crossing-port overlap contained in the bore/host. The supplied
`simple_parametric_input_new.yaml` is supported for both conical and
perpendicular vessels, including all 29 component cells. If changing dimensions
makes ports meet conical transitions, clip within a host, intersect LF rings
partially, or leave the crossing overlap outside the assumed host, the API
raises `ValueError` rather than substituting an approximate volume. Additional
analytical intersection branches would be needed for those configurations.

Existing geometry is deliberately preserved: OpenMC Intersection.__iand__
mutates regions in place, so the port builder subtracts vessel vacuum from its
stored port regions. Port volumes exclude that intersection analytically.
HF and end cells do not receive port subtractions.
These records do not validate the overall geometry as overlap-free.

Tests:

    PYTHONDONTWRITEBYTECODE=1 /home/hrshah3/miniconda3/bin/python -m unittest discover -s tests -p 'test_*.py' -v

Reference quadrature appears only in tests, checking the oblique-cylinder
expression and end-cell volumes directly from axisymmetric OpenMC CSG.

Component regression coverage includes:

- Separate left/right HF magnet volumes, every casing layer and the shield,
  compared to exact interval sums from actual OpenMC surfaces and membership.
- Zero, one and three HF casing layers, plus volume conservation across the
  union of magnet/casing/shield regions.
- The 35–45 cm magnet, 50 cm long, with 5 cm radial and axial shielding:
  each separate shield record is 56,000*pi cm³.
- Every vessel and central-cell layer for both vessel styles, with independent
  frustum and port-cut references.
- LF magnet/shield references, all port layers on both sides, and end-cell
  vessel exclusions for both vessel styles.
- Cubic scaling of components when dimensions double. End cells are checked
  against their actual scaled CSG instead: the builder's fixed +10 cm placement
  offset means their vessel overlap is not geometrically similar under scaling.
- Explicit rejection of unsupported clipping, transition and overlap cases.
