import openmc
import matplotlib.pyplot as plt
from geometry_lib import redefined_vacuum_vessel_region  # adjust if it's in a different file

# ---------------- Parameters ---------------- #
outer_axial_length = 100
central_axial_length = 350
central_radius = 75
bottleneck_radius = 30
left_bottleneck_length = 50
right_bottleneck_length = 50
axial_midplane = 0.0

# ---------------- Build Region ---------------- #
vv_region, components = redefined_vacuum_vessel_region(
    outer_axial_length,
    central_axial_length,
    central_radius,
    bottleneck_radius,
    left_bottleneck_length,
    right_bottleneck_length,
    axial_midplane
)

# ---------------- Create Cell & Universe ---------------- #
vv_cell = openmc.Cell(region=vv_region)
# vv_cell.fill = openmc.Material(name='vacuum')  # just a placeholder
universe = openmc.Universe(cells=[vv_cell])

# ---------------- Plot ---------------- #
plot = vv_region.plot(
    basis='xz',
    width=(200, 500),  # adjust as needed
    pixels=(1000, 300)
)
plt.title("Vacuum Vessel Test")
plt.savefig("vv_region_test.png", bbox_inches="tight")
print("✅ Plot saved as vv_region_test.png")
