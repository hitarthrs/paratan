#!/openmc_venv/bin/python
#
import openmc
import bohmPrepare
#
# create an openmc model of a simple cylindrical tokamak with FW armor only
#
# useful functions
def get_material(materials, material):
    """
    Search the materials object for a material with a matching name. Openmc
    allows duplicate names, and names are not required, be advised.

    Arguments:
        materials (OpenMC Materials Object): material library to search
        material (string): name of material to be returned

    Returns:
        mat (OpenMC material object): material object with matching name
    """
    for mat in materials:
        if mat.name == material:
            return mat
    # if this returns none, openmc will just assign vacuum to any cell
    # using this material
    raise ValueError(f"no material name {material} was found in the library")
#
# materials section
# using an existing materials.xml file so no need to create
# my_materials=openmc.Materials.from_xml(path='materials_orig.xml') # initialize and load the model
# print(get_material(my_materials,"tungsten")) # test if loading works using Edgar's function

tungsten = openmc.Material(31, name="tungsten")
tungsten.set_density("g/cm3", 19.25)
tungsten.add_element("W", 100)

my_materials = openmc.Materials([tungsten])


my_materials.export_to_xml(path='materials.xml') # write the materials.xml file
#
# geometry
#
# create surfaces for FW armor (Hitarth model cellid=1001) at:
#        central cell cyl (r=75 cm id=9 z-cyl) surface for testing
#        central cell cyl (r=22.3 cm id=14 z-cyl) surface for testing
#        transition cones (are the top and bottom separate surfaces? id=15,17 at +/- 407.87)
#
# surfaces
model_top=openmc.ZPlane(z0=1500, boundary_type='vacuum') 
model_bottom=openmc.ZPlane(z0=-1500, boundary_type='vacuum') 
model_side_outer_radius=openmc.ZCylinder(r=200.0, boundary_type='vacuum') 
#
centralfwarmor_top=openmc.ZPlane(z0=127.9675)
centralfwarmor_bottom=openmc.ZPlane(z0=-127.9675)
#
upperendfwarmor_top=openmc.ZPlane(z0=597.5)
upperendfwarmor_bottom=openmc.ZPlane(z0=325)
#
lowerendfwarmor_top=openmc.ZPlane(z0=-325)
lowerendfwarmor_bottom=openmc.ZPlane(z0=-597.5)
#
centralplasma_outer_radius=openmc.ZCylinder(r=75.0)
centralfwarmor_outer_radius=openmc.ZCylinder(r=75.2)
#
endplasma_outer_radius=openmc.ZCylinder(r=22.3)
endfwarmor_outer_radius=openmc.ZCylinder(r=22.5)
#
lowerconeplasma_outer_radius=openmc.model.ZConeOneSided(x0=0,y0=0,z0=-407.8713105676658,r2=0.07179676972449082,up=True)
lowerconefwarmor_outer_radius=openmc.model.ZConeOneSided(x0=0,y0=0,z0=-408.6177207291796,r2=0.07179676972449082,up=True)
#
upperconeplasma_outer_radius=openmc.model.ZConeOneSided(x0=0,y0=0,z0=407.8713105676658,r2=0.07179676972449082,up=False)
upperconefwarmor_outer_radius=openmc.model.ZConeOneSided(x0=0,y0=0,z0=408.6177207291796,r2=0.07179676972449082,up=False)
#
#
# regions
# plasma regions
lowerendplasma_region=-endplasma_outer_radius & +lowerendfwarmor_bottom & -lowerendfwarmor_top
lowerconeplasma_region=-lowerconeplasma_outer_radius & -centralfwarmor_bottom & +lowerendfwarmor_top
centralplasma_region=-centralplasma_outer_radius & -centralfwarmor_top & +centralfwarmor_bottom
upperconeplasma_region=-upperconeplasma_outer_radius & +centralfwarmor_top & -upperendfwarmor_bottom
upperendplasma_region=-endplasma_outer_radius & +upperendfwarmor_bottom & -upperendfwarmor_top
# fwarmor regions
lowerendfwarmor_region=+endplasma_outer_radius & -endfwarmor_outer_radius & +lowerendfwarmor_bottom & -lowerendfwarmor_top
lowerconefwarmor_region=+lowerconeplasma_outer_radius & -lowerconefwarmor_outer_radius & -centralfwarmor_bottom & +lowerendfwarmor_top
centralfwarmor_region=+centralplasma_outer_radius & -centralfwarmor_outer_radius & -centralfwarmor_top & +centralfwarmor_bottom
upperconefwarmor_region=+upperconeplasma_outer_radius & -upperconefwarmor_outer_radius & +centralfwarmor_top & -upperendfwarmor_bottom
upperendfwarmor_region=+endplasma_outer_radius & -endfwarmor_outer_radius & +upperendfwarmor_bottom & -upperendfwarmor_top
# outside regions
lowerdisk_region=+model_bottom & -lowerendfwarmor_bottom & -model_side_outer_radius
lowerendoutside_region=+endfwarmor_outer_radius & -model_side_outer_radius & +lowerendfwarmor_bottom & -lowerendfwarmor_top
lowerconeoutside_region=+lowerconefwarmor_outer_radius & -model_side_outer_radius & -centralfwarmor_bottom & +lowerendfwarmor_top
centraloutside_region=+centralfwarmor_outer_radius & -model_side_outer_radius & -centralfwarmor_top & +centralfwarmor_bottom
upperconeoutside_region=+upperconefwarmor_outer_radius & -model_side_outer_radius & +centralfwarmor_top & -upperendfwarmor_bottom
upperendoutside_region=+endfwarmor_outer_radius & -model_side_outer_radius & +upperendfwarmor_bottom & -upperendfwarmor_top
upperdisk_region=-model_top & +upperendfwarmor_top & -model_side_outer_radius
#
# create cell and assign material
lowerendplasma_cell=openmc.Cell(region=lowerendplasma_region,fill=None)
lowerconeplasma_cell=openmc.Cell(region=lowerconeplasma_region,fill=None)
centralplasma_cell=openmc.Cell(region=centralplasma_region,fill=None)
upperconeplasma_cell=openmc.Cell(region=upperconeplasma_region,fill=None)
upperendplasma_cell=openmc.Cell(region=upperendplasma_region,fill=None)
#
lowerendfwarmor_cell=openmc.Cell(region=lowerendfwarmor_region,fill=None)
lowerconefwarmor_cell=openmc.Cell(region=lowerconefwarmor_region,fill=None)
centralfwarmor_cell=openmc.Cell(region=centralfwarmor_region,fill=None)
upperconefwarmor_cell=openmc.Cell(region=upperconefwarmor_region,fill=None)
upperendfwarmor_cell=openmc.Cell(region=upperendfwarmor_region,fill=None)
#
lowerdisk_cell=openmc.Cell(region=lowerdisk_region,fill=None)
lowerendoutside_cell=openmc.Cell(region=lowerendoutside_region,fill=None)
lowerconeoutside_cell=openmc.Cell(region=lowerconeoutside_region,fill=None)
centraloutside_cell=openmc.Cell(region=centraloutside_region,fill=None)
upperconeoutside_cell=openmc.Cell(region=upperconeoutside_region,fill=None)
upperendoutside_cell=openmc.Cell(region=upperendoutside_region,fill=None)
upperdisk_cell=openmc.Cell(region=upperdisk_region,fill=None)
#
# make a universe to contain all the cells
geometry = openmc.Geometry([lowerendplasma_cell, lowerconeplasma_cell, centralplasma_cell, upperconeplasma_cell, upperendplasma_cell,
                            lowerendfwarmor_cell, lowerconefwarmor_cell, centralfwarmor_cell, upperconefwarmor_cell, upperendfwarmor_cell,
                            lowerdisk_cell,
                            lowerendoutside_cell, lowerconeoutside_cell, centraloutside_cell, upperconeoutside_cell, upperendoutside_cell,
                            upperdisk_cell])
#
geometry.export_to_xml() # writes the geometry.xml file
#
# plot section
#
# generate some plot slices
plot_xz_overall= geometry.plot(basis='xz', color_by='cell', origin=(0,0,0), width=(1200,1200), pixels=(1000,1000))
plot_xz_overall.figure.savefig('xz-cell_overall.png')
#
plot_xz_lower = geometry.plot(basis='xz', color_by='cell', origin=(0,0,-300), width=(100,100), pixels=(1000,1000))
plot_xz_lower.figure.savefig('xz-cell_lower.png')
#
plot_xz_upper = geometry.plot(basis='xz', color_by='cell', origin=(0,0,300), width=(100,100), pixels=(1000,1000))
plot_xz_upper.figure.savefig('xz-cell_upper.png')
#
#
# source and settings section
# create source from function which uses anvil 2D data
# (may need to modify source strength for 1.638e18 n/sec)
#
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.source = bohmPrepare.sources_2d
settings.batches = 10
settings.particles = 100000
#
# write surface source for FW armor (cellids=6 to 10) to enable "tally" for NWL

settings.surf_source_write = {
    'surface_ids': [12,14,10,18],
    'max_particles': 1000000
} # test writing lower end (actually writes lower and upper ends), lower cone, central, upper cone 

settings.export_to_xml()
#
# tallies section
cell_filter = openmc.CellFilter([6,7,8,9,10])
mytally = openmc.Tally(14) # set tally number to 14 (arbitrary)
mytally.filters = [cell_filter]
mytally.scores = ['flux']
#
tallies = openmc.Tallies([mytally])
tallies.export_to_xml()
#
#
# now run
openmc.run()
#
# post process
