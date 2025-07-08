# %%

## Author: Mason Yu (Email ID: myu233@wisc.edu)
## Contributors: Hitarth Shah (Email ID: hitarth.shah@wisc.edu) and Tim Bohm

import openmc
import matplotlib.pyplot as plt
from pkg_resources import ZipProvider
import openmc.deplete

length = 60
length_cyl = int(length)

vacuum = openmc.Material(0, name='vacuum')
vacuum.set_density('g/cm3', 1e-40)
vacuum.add_nuclide('N14', 1)

air = openmc.Material(4, name='air')
air.set_density('g/cm3', 0.001205)
air.add_nuclide('C12', 1.48583129e-04)
air.add_nuclide('C13', 1.60703476e-06)
air.add_nuclide('N14', 7.81575372e-01)
air.add_nuclide('N15', 2.85532775e-03)
air.add_nuclide('O16', 2.10235865e-01)
air.add_nuclide('O17', 8.00842335e-05)
air.add_nuclide('O18', 4.32033365e-04)
air.add_nuclide('Ar36', 1.55828794e-05)
air.add_nuclide('Ar38', 2.93813883e-06)
air.add_nuclide('Ar40', 4.65260590e-03)

deuterium = openmc.Material(1, name='deuterium')
deuterium.set_density('kg/m3', 0.001)
deuterium.add_nuclide('H2', 1)

stainless = openmc.Material(298, name = 'stainless 304')
stainless.set_density('g/cm3', 8.0)
stainless.add_nuclide('C12', 1.81016539e-03)
stainless.add_nuclide('C13', 1.95782570e-05)
stainless.add_nuclide('Si28', 9.02042466e-03)
stainless.add_nuclide('Si29', 4.58244576e-04)
stainless.add_nuclide('Si30', 3.02431639e-04)
stainless.add_nuclide('P31', 4.07975189e-04)
stainless.add_nuclide('S32', 2.44140964e-04)
stainless.add_nuclide('S33', 1.92763157e-06)
stainless.add_nuclide('S34', 1.09232456e-05)
stainless.add_nuclide('S36', 2.57017543e-08)
stainless.add_nuclide('Cr50', 8.72312748e-03)
stainless.add_nuclide('Cr52', 1.68216831e-01)
stainless.add_nuclide('Cr53', 1.90744383e-02)
stainless.add_nuclide('Cr54', 4.74803142e-03)
stainless.add_nuclide('Mn55', 1.00006144e-02)
stainless.add_nuclide('Fe54', 4.03523669e-02)
stainless.add_nuclide('Fe56', 6.33445864e-01)
stainless.add_nuclide('Fe57', 1.46290275e-02)
stainless.add_nuclide('Fe58', 1.94685500e-03)
stainless.add_nuclide('Ni58', 5.89458369e-02)
stainless.add_nuclide('Ni60', 2.27057109e-02)
stainless.add_nuclide('Ni61', 9.87005296e-04)
stainless.add_nuclide('Ni62', 3.14709137e-03)
stainless.add_nuclide('Ni64', 8.01362752e-04)
#stainless.add_s_alpha_beta('c_Fe56')

aluminum_6061 = openmc.Material(13, name='aluminum_6061')
aluminum_6061.set_density('g/cm3', 2.7)
aluminum_6061.add_nuclide('Mg24', 8.81687922e-03)
aluminum_6061.add_nuclide('Mg25', 1.11620195e-03)
aluminum_6061.add_nuclide('Mg26', 1.22893835e-03)
aluminum_6061.add_nuclide('Al27', 9.77324713e-01)
aluminum_6061.add_nuclide('Si28', 5.34499966e-03)
aluminum_6061.add_nuclide('Si29', 2.71530132e-04)
aluminum_6061.add_nuclide('Si30', 1.79204092e-04)
aluminum_6061.add_nuclide('Ti46', 4.11473671e-05)
aluminum_6061.add_nuclide('Ti47', 3.71074438e-05)
aluminum_6061.add_nuclide('Ti48', 3.67682897e-04)
aluminum_6061.add_nuclide('Ti49', 2.69826977e-05)
aluminum_6061.add_nuclide('Ti50', 2.58355590e-05)
aluminum_6061.add_nuclide('Cr50', 4.42071668e-05)
aluminum_6061.add_nuclide('Cr52', 8.52491208e-04)
aluminum_6061.add_nuclide('Cr53', 9.66656598e-05)
aluminum_6061.add_nuclide('Cr54', 2.40621288e-05)
aluminum_6061.add_nuclide('Mn55', 4.34559057e-04)
aluminum_6061.add_nuclide('Fe54', 1.16134627e-04)
aluminum_6061.add_nuclide('Fe56', 1.82306529e-03)
aluminum_6061.add_nuclide('Fe57', 4.21025279e-05)
aluminum_6061.add_nuclide('Fe58', 5.60307355e-06)
aluminum_6061.add_nuclide('Cu63', 8.11849846e-04)
aluminum_6061.add_nuclide('Cu65', 3.62191869e-04)
aluminum_6061.add_nuclide('Zn64', 2.97894307e-04)
aluminum_6061.add_nuclide('Zn66', 1.68000999e-04)
aluminum_6061.add_nuclide('Zn67', 2.44761643e-05)
aluminum_6061.add_nuclide('Zn68', 1.11778523e-04)
aluminum_6061.add_nuclide('Zn70', 3.69565848e-06)
aluminum_6061.add_s_alpha_beta('c_Al27')

# This is basically YBCO for now
rebco = openmc.Material(name="REBCO tape")
rebco.set_density("g/cm3", 6.3)
rebco.add_element('Y', 7.6923076)
rebco.add_element('Ba', 15.3846153)
rebco.add_element('Cu', 23.0769230)
rebco.add_element('O', 53.8461538)

# Hastelloy C-276 subtrate, composition from Haynes International
hastelloy = openmc.Material(name="Hastelloy C-276")
hastelloy.set_density("g/cm3", 8.89)
hastelloy.add_element("Ni", 55, 'wo')
hastelloy.add_element("Co", 2.5, 'wo')
hastelloy.add_element("Cr", 16, 'wo')
hastelloy.add_element("Mo", 16, 'wo')
hastelloy.add_element("Fe", 5, 'wo')
hastelloy.add_element("W", 4, 'wo')
hastelloy.add_element("Mn", 1, "wo")
hastelloy.add_element("V", 0.35, "wo")
hastelloy.add_element("Cu", 0.15, "wo")

copper = openmc.Material(name="Copper")
copper.set_density("g/cm3", 8.96)
copper.add_element("Cu", 100)

silver = openmc.Material(name='Silver')
silver.set_density('g/cm3', 10.49)
silver.add_element("Ag", 100)

magnet = openmc.Material.mix_materials([hastelloy, copper, rebco], [0.55, 0.43, 0.02], 'vo')

tungsten = openmc.Material(31, name="tungsten")
tungsten.set_density("g/cm3", 19.25)
tungsten.add_element("W", 100)

crispy = openmc.Material(81, name="crispy")
crispy.set_density('g/cm3', 1.4)
crispy.add_nuclide('H1', 6.72038450e-02)
crispy.add_nuclide('H2', 7.72933104e-06)
crispy.add_nuclide('B10', 1.36943828e-01)
crispy.add_nuclide('B11', 5.51216112e-01)
crispy.add_nuclide('C12', 2.26052782e-01)
crispy.add_nuclide('C13', 2.44492547e-03)
crispy.add_nuclide('O16', 1.34096500e-02)
crispy.add_nuclide('O17', 5.10807965e-06)
crispy.add_nuclide('O18', 2.75567455e-05)
crispy.add_nuclide('Cl35', 2.75567455e-05)
crispy.add_nuclide('Cl37', 6.51683424e-04)

# This water does not have the S_ab scattering kernels, meant for mixing
water = openmc.Material(101, name="water")
water.set_density('g/cm3', 1)
water.add_element('H', 66.666)
water.add_element('O', 33.333)

cooled_tungsten = openmc.Material.mix_materials([tungsten, water], [0.8, 0.2], 'vo')
cooled_tungsten.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

tungsten_carbide = openmc.Material(22, name = 'tungsten carbide')
tungsten_carbide.add_elements_from_formula('WC')
tungsten_carbide.set_density('g/cm3', 15.63)

# LaMnO3 (Lanthanum Manganite)
# Source for Density: NIST (https://www.ctcms.nist.gov/~knc6/jsmol/JVASP-11795.html) 
lamno3 = openmc.Material(name='LaMnO3')
lamno3.add_elements_from_formula('LaMnO3')  # Lanthanum Manganite
lamno3.set_density('g/cm3', 6.67)


# Y2O3 (Yttrium Oxide)
yttrium_oxide = openmc.Material(name='Yttrium Oxide: Y2O3')
yttrium_oxide.add_elements_from_formula('Y2O3')
yttrium_oxide.set_density('g/cm3', 5.01)

# Al2O3 (Aluminum Oxide or Alumina)
alumina = openmc.Material(name='Alumina: Al2O3')
alumina.add_elements_from_formula('Al2O3')   # Aluminum Oxide formula
alumina.set_density('g/cm3', 3.99)


# rafm_steel = openmc.Material(900, name="EUROFER 97 RAFM Steel")
# rafm_steel.add_element('Fe', 90)
# rafm_steel.add_element('Cr', 9.21)
# rafm_steel.add_element('C', 0.104)
# rafm_steel.add_element('Mn', 0.502)
# rafm_steel.add_element('V', 0.204)
# rafm_steel.add_element('W', 1.148)
# rafm_steel.add_element('Ta', 0.14 )
# rafm_steel.add_nuclide('N14', 0.0234)
# rafm_steel.add_element('O', 0.001)
# rafm_steel.add_element('P', 0.04)
# rafm_steel.add_element('S', 0.004)
# rafm_steel.add_element('B', 0.01)
# rafm_steel.add_element('Ti', 0.004)
# rafm_steel.add_element('Nb', 0.0012)
# rafm_steel.add_element('Mo', 0.008)
# rafm_steel.add_element('Ni', 0.0214)
# rafm_steel.set_density('g/cm3', 7.798)


## According to Tim's Materials Library
rafm_steel = openmc.Material(name='Eurofer97')
rafm_steel.set_density('g/cm3', 7.75)
rafm_steel.add_nuclide('C12', 1.1545e-03, 'wo')
rafm_steel.add_nuclide('C13', 1.3531e-05, 'wo')
rafm_steel.add_nuclide('V50', 4.9454e-06, 'wo')
rafm_steel.add_nuclide('V51', 2.0126e-03, 'wo')
rafm_steel.add_nuclide('Cr52', 7.9100e-02, 'wo')
rafm_steel.add_nuclide('Cr53', 9.1420e-03, 'wo')
rafm_steel.add_nuclide('Cr54', 2.3186e-03, 'wo')
rafm_steel.add_nuclide('Fe56', 8.6988e-01, 'wo')
rafm_steel.add_nuclide('Fe57', 2.0449e-02, 'wo')
rafm_steel.add_nuclide('Fe58', 2.7690e-03, 'wo')
rafm_steel.add_nuclide('Ta180', 1.7755e-07, 'wo')
rafm_steel.add_nuclide('Ta181', 1.4864e-03, 'wo')
rafm_steel.add_nuclide('W182', 3.0634e-03, 'wo')
rafm_steel.add_nuclide('W183', 1.6634e-03, 'wo')
rafm_steel.add_nuclide('W184', 3.5810e-03, 'wo')
rafm_steel.add_nuclide('W186', 3.3589e-03, 'wo')


cooled_rafm_steel = openmc.Material.mix_materials([rafm_steel, water], [0.8, 0.20], 'vo')
cooled_rafm_steel.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

enrichment_li6 = 0.15
enrichment_li7 = 1- enrichment_li6

LiPb_density = 9965
pb_at_percent = 84.3

LiPb_breeder = openmc.Material(1000, name = "lead-lithium-eutectic-breeder")
# Should double check this number... Varies with temperature
LiPb_breeder.set_density('kg/m3', LiPb_density)
LiPb_breeder.add_element('Li', 100-pb_at_percent, percent_type='ao', enrichment=enrichment_li6*100,
                          enrichment_target='Li6', enrichment_type='ao')
#LiPb_breeder.add_element('Li', 17, percent_type='ao')
LiPb_breeder.add_element('Pb', pb_at_percent)
LiPb_breeder.depletable = True

struct_material = 0.05 
# Shield material with structural support and coolant mixed in
cooled_tungsten_carbide = openmc.Material.mix_materials([tungsten_carbide, water, rafm_steel], [0.5-struct_material, 0.5, struct_material], 'vo')
cooled_tungsten_carbide.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.5)

water_cooled_wc = openmc.Material.mix_materials([tungsten_carbide, water], [0.8, 0.2], 'vo')
water_cooled_wc.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

# Shield material with structural support and coolant mixed in
cooled_stainless = openmc.Material.mix_materials([stainless, water], [0.8, 0.2], 'vo')
cooled_stainless.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

helium_8mpa = openmc.Material(1100, name="helium gas 8 MPa, 450C")
helium_8mpa.set_density("kg/m3", 5.323)
helium_8mpa.add_element("He", 100)

# Portland Concrete
# Source: https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=144

# Create a material for concrete
concrete = openmc.Material(name='Concrete')

# Add elements by atomic number (Z) and weight fraction (wo)
concrete.add_element('H', 0.010000, percent_type='wo')   # Z = 1, 1.00%
concrete.add_element('C', 0.001000, percent_type='wo')   # Z = 6, 0.10%
concrete.add_element('O', 0.529107, percent_type='wo')   # Z = 8, 52.91%
concrete.add_element('Na', 0.016000, percent_type='wo')  # Z = 11, 1.60%
concrete.add_element('Mg', 0.002000, percent_type='wo')  # Z = 12, 0.20%
concrete.add_element('Al', 0.033872, percent_type='wo')  # Z = 13, 3.39%
concrete.add_element('Si', 0.337021, percent_type='wo')  # Z = 14, 33.70%
concrete.add_element('K', 0.013000, percent_type='wo')   # Z = 19, 1.30%
concrete.add_element('Ca', 0.044000, percent_type='wo')  # Z = 20, 4.40%
concrete.add_element('Fe', 0.014000, percent_type='wo')  # Z = 26, 1.40%

# Set the density (in g/cm^3)
concrete.set_density('g/cm3', 2.3)

he_cooled_rafm = openmc.Material.mix_materials([helium_8mpa, rafm_steel], [0.6, 0.4], 'vo')
#print(he_cooled_rafm.density)

he_cooled_tungsten_carbide = openmc.Material.mix_materials([tungsten_carbide, helium_8mpa], [0.6, 0.4], 'vo')

rings = openmc.Material.mix_materials([stainless, vacuum], [0.1, 0.9], 'vo')

tungsten_boride = openmc.Material(1200, name="tungsten boride WB")
tungsten_boride.set_density("g/cm3", 15.43)
tungsten_boride.add_elements_from_formula('WB')

WB2 = openmc.Material(1202, name="tungsten boride WB2")
WB2.set_density("g/cm3", 12.42)
WB2.add_elements_from_formula("WB2")

w2b5 = openmc.Material(1201, name="tungsten boride W2B5")
w2b5.set_density("g/cm3", 12.91)
w2b5.add_elements_from_formula("W2B5")

cooled_w2b5 = openmc.Material.mix_materials([w2b5, water], [0.8, 0.2], 'vo')
cooled_w2b5.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

cooled_wb2 = openmc.Material.mix_materials([WB2, water], [0.8, 0.2], 'vo')
cooled_wb2.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

TiH2 = openmc.Material(1300, name="titanium hydride")
TiH2.set_density("g/cm3", 3.75)
TiH2.add_elements_from_formula("TiH2")

zirconium_hydride = openmc.Material(1400, name="zirconium hydride ZrH2")
zirconium_hydride.set_density("g/cm3", 5.56)
zirconium_hydride.add_elements_from_formula("ZrH2")

tantalum = openmc.Material(name="pure tantalum")
tantalum.set_density("g/cm3", 16.69)
tantalum.add_element("Tantalum", 100)

cooled_TiH2 = openmc.Material.mix_materials([TiH2, water], [0.95, 0.05], 'vo')

flibe = openmc.Material(name="Pure FLiBe molten salt")
flibe.set_density("kg/m3", 1940)
flibe.add_elements_from_formula("Li2BeF4")

beryllium = openmc.Material(name="Pure Beryllium")
beryllium.set_density("g/cm3", 1.85)
beryllium.add_element("Be", 100)

# Define 1050 Aluminum material according to the Aluminum Assocation. Follow link: https://www.aluminum.org/sites/default/files/2021-10/Teal%20Sheet.pdf
aluminum_1050 = openmc.Material(name='Aluminum 1050')
aluminum_1050.add_element('Al', 99.5, percent_type='wo')  # weight percent
aluminum_1050.add_element('Si', 0.25, percent_type='wo')
aluminum_1050.add_element('Fe', 0.25, percent_type='wo')

tantalum_hydride_55 = openmc.Material(name='tantalum hydride, TaH0.55')
tantalum_hydride_55.set_density("g/cm3", 16.69)
tantalum_hydride_55.add_element("Tantalum", 1/1.55)
tantalum_hydride_55.add_element("Hydrogen", 0.55/1.55)
tantalum_hydride_55.temperature = 500.0

tantalum_hydride_46 = openmc.Material(name='tantalum hydride, TaH0.46')
tantalum_hydride_46.set_density("g/cm3", 16.69)
tantalum_hydride_46.add_element("Tantalum", 1/1.46)
tantalum_hydride_46.add_element("Hydrogen", 0.46/1.46)
tantalum_hydride_46.temperature = 500.0

tantalum_hydride_39 = openmc.Material(name='tantalum hydride, TaH0.39')
tantalum_hydride_39.set_density("g/cm3", 16.69)
tantalum_hydride_39.add_element("Tantalum", 1/1.39)
tantalum_hydride_39.add_element("Hydrogen", 0.39/1.39)
tantalum_hydride_39.temperature = 500.0

tantalum_hydride_30 = openmc.Material(name='tantalum hydride, TaH0.30')
tantalum_hydride_30.set_density("g/cm3", 16.69)
tantalum_hydride_30.add_element("Tantalum", 1/1.3)
tantalum_hydride_30.add_element("Hydrogen", 0.3/1.3)
#tantalum_hydride_30.temperature = 500.0

Nak_77 = openmc.Material(name="NaK eutectic, 550C")
Nak_77.set_density("g/cm3", 0.749)
Nak_77.add_element("Na", 23, 'wo')
Nak_77.add_element("K", 77, 'wo')
Nak_77.temperature = 900
Nak_77.depletable = True

potassium = openmc.Material(name="K molten, 550C")
potassium.set_density("g/cm3", 0.82948 )
potassium.add_element("K", 100, 'ao')
potassium.temperature = 900
potassium.depletable = True

KCl = openmc.Material(name="KCl molten, 977C")
KCl.set_density("g/cm3", 1.527)
KCl.add_elements_from_formula("KCl", enrichment=90, enrichment_target="Cl37")
#KCl.add_elements_from_formula("KCl")
KCl.temperature = 1200
KCl.depletable = True

iron = openmc.Material(name='pure iron')
iron.set_density("g/cm3", 7.874)
iron.add_element("Fe", 100, 'ao')

LiH = openmc.Material(name='Lithium Hydride LiH breeder')
LiH.set_density("g/cm3", 0.775)
LiH.add_element("Li", enrichment=4.85, percent=50, enrichment_target='Li6', enrichment_type='ao')
LiH.add_element("H", 50, 'ao')
LiH.temperature = 900

LiD = openmc.Material(name='Lithium Deuteride LiD breeder')
LiD.set_density("g/cm3", 0.885)
LiD.add_element("Li", enrichment=4.85, percent=50, enrichment_target='Li6', enrichment_type='ao')
LiD.add_element("H", enrichment=100, percent=50, enrichment_target='H2', enrichment_type='ao')
LiD.temperature = 900

lead = openmc.Material(name='pure lead')
lead.set_density("g/cm3", 11.34)
lead.add_element("Pb", 100)
lead.temperature = 900

lithium = openmc.Material(name='molten lithium')
lithium.set_density("g/cm3", 0.512)
lithium.add_element("Li", 100)
lithium.temperature = 900


li6 = openmc.Nuclide("Li6")
li7 = openmc.Nuclide("Li7")
enriched_li_density = 0.512*(enrichment_li6*6.015+enrichment_li7*6.941)/6.8716

enriched_lithium = openmc.Material(name='enriched lithium')
enriched_lithium.set_density("g/cm3", enriched_li_density)
enriched_lithium.add_element('Li', 1.0, enrichment=enrichment_li6*100, enrichment_target='Li6')
#enriched_lithium.add_nuclide(li6, enrichment_li6, 'ao')
#enriched_lithium.add_nuclide(li7, 1-enrichment_li6, 'ao')
enriched_lithium.temperature = 900

HfH2 = openmc.Material(name="Hafnium hydride HfH2")
HfH2.set_density("g/cm3", 11.36)
HfH2.add_elements_from_formula("HfH2")

titanium = openmc.Material(name="pure titanium")
titanium.set_density("g/cm3", 4.506)
titanium.add_element("Ti", 100)

MgO = openmc.Material(name="Magnesium oxide MgO")
MgO.set_density("g/cm3", 3.58)
MgO.add_elements_from_formula("MgO")

packing_fraction_pebbles = 0.7

lithium_metatitanate = openmc.Material(name="Lithium metatitanate")
lithium_metatitanate.set_density("g/cm3", packing_fraction_pebbles*3.43)
lithium_metatitanate.add_nuclide(li6, 2 * enrichment_li6 , 'ao')
lithium_metatitanate.add_nuclide(li7, 2 * enrichment_li7 , 'ao')
lithium_metatitanate.add_element("Ti", 1, 'ao')
lithium_metatitanate.add_element('O', 3, 'ao')
#lithium_metatitanate.temperature = 500

lithium_orthosilicate = openmc.Material(name="Lithium orthosilicate")
lithium_orthosilicate.set_density("g/cm3", packing_fraction_pebbles*2.39)
lithium_orthosilicate.add_nuclide(li6, 4 * enrichment_li6 , 'ao')
lithium_orthosilicate.add_nuclide(li7, 4 * enrichment_li7 , 'ao')
lithium_orthosilicate.add_element('Si', 1, 'ao')
lithium_orthosilicate.add_element('O', 4, 'ao')
#lithium_orthosilicate.temperature = 400



# Define SS316LN material from Fusion materails database by Tim Bohhm
ss316ln = openmc.Material(name='SS316LN')

# Set the material density
ss316ln.set_density('g/cm3', 7.93)

# Add nuclides with their atomic fractions (ao)
ss316ln.add_nuclide('B10', 3.0414e-04, percent_type='ao')
ss316ln.add_nuclide('B11', 1.2242e-03, percent_type='ao')
ss316ln.add_nuclide('C12', 1.3609e-03, percent_type='ao')
ss316ln.add_nuclide('C13', 1.4720e-05, percent_type='ao')
ss316ln.add_nuclide('N14', 6.2685e-03, percent_type='ao')
ss316ln.add_nuclide('N15', 2.2901e-05, percent_type='ao')
ss316ln.add_nuclide('Si28', 1.8085e-02, percent_type='ao')
ss316ln.add_nuclide('Si29', 9.1873e-04, percent_type='ao')
ss316ln.add_nuclide('Si30', 6.0634e-04, percent_type='ao')
ss316ln.add_nuclide('P31', 5.3344e-04, percent_type='ao')
ss316ln.add_nuclide('S32', 3.2632e-04, percent_type='ao')
ss316ln.add_nuclide('S33', 2.5765e-06, percent_type='ao')
ss316ln.add_nuclide('S34', 1.4600e-05, percent_type='ao')
ss316ln.add_nuclide('S36', 3.4353e-08, percent_type='ao')
ss316ln.add_nuclide('Cr50', 7.9390e-03, percent_type='ao')
ss316ln.add_nuclide('Cr52', 1.5310e-01, percent_type='ao')
ss316ln.add_nuclide('Cr53', 1.7360e-02, percent_type='ao')
ss316ln.add_nuclide('Cr54', 4.3212e-03, percent_type='ao')
ss316ln.add_nuclide('Mn55', 2.0050e-02, percent_type='ao')
ss316ln.add_nuclide('Fe54', 3.7371e-02, percent_type='ao')
ss316ln.add_nuclide('Fe56', 5.8665e-01, percent_type='ao')
ss316ln.add_nuclide('Fe57', 1.3548e-02, percent_type='ao')
ss316ln.add_nuclide('Fe58', 1.8030e-03, percent_type='ao')
ss316ln.add_nuclide('Co59', 9.3454e-04, percent_type='ao')
ss316ln.add_nuclide('Ni58', 7.6657e-02, percent_type='ao')
ss316ln.add_nuclide('Ni60', 2.9528e-02, percent_type='ao')
ss316ln.add_nuclide('Ni61', 1.2836e-03, percent_type='ao')
ss316ln.add_nuclide('Ni62', 4.0927e-03, percent_type='ao')
ss316ln.add_nuclide('Ni64', 1.0421e-03, percent_type='ao')
ss316ln.add_nuclide('Nb93', 2.9640e-04, percent_type='ao')
ss316ln.add_nuclide('Mo92', 2.0849e-03, percent_type='ao')
ss316ln.add_nuclide('Mo94', 1.3129e-03, percent_type='ao')
ss316ln.add_nuclide('Mo95', 2.2728e-03, percent_type='ao')
ss316ln.add_nuclide('Mo96', 2.3919e-03, percent_type='ao')
ss316ln.add_nuclide('Mo97', 1.3775e-03, percent_type='ao')
ss316ln.add_nuclide('Mo98', 3.4996e-03, percent_type='ao')
ss316ln.add_nuclide('Mo100', 1.4090e-03, percent_type='ao')





MgO_HfH2 = openmc.Material.mix_materials([MgO, HfH2, vacuum], [0.4, 0.4, 0.2], 'vo', name="MgO 40% vo, HfH2 40% vo, void 20% vo, Snead alloy 1")

Fe_HfH2_WB2 = openmc.Material.mix_materials([iron, HfH2, WB2, vacuum], [0.54, 0.34, 0.07, 0.05], 'vo', name="Fe 54 vol % HfH2 34 vol % WB2 7 vol % remainder void, Snead alloy 2")

Ti_HfH2 = openmc.Material.mix_materials([titanium, HfH2, vacuum], [0.5, 0.48, 0.02], 'vo', name="Ti 50 vol % HfH2 48vol % remainder void, Snead alloy 3")

B4C = openmc.Material(name="Boron carbide B4C")
B4C.set_density("g/cm3", 2.5)
B4C.add_elements_from_formula("B4C")

cooled_tungsten_boride = openmc.Material.mix_materials([tungsten_boride, water], [0.8, 0.2], 'vo')
cooled_tungsten_boride.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)

cooled_B4C = openmc.Material.mix_materials([B4C, water], [0.8, 0.2], 'vo')
cooled_B4C.add_s_alpha_beta('c_H_in_H2O', 0.66666*0.2)


HT_Shield_filler = openmc.Material.mix_materials([tungsten_carbide, helium_8mpa, rafm_steel], [0.7, 0.25, 0.05], 'vo', name="HT Shield Filler")
LT_Shield_filler = openmc.Material.mix_materials([tungsten_carbide, helium_8mpa, rafm_steel], [0.35, 0.35, 0.3], 'vo', name="LT Shield Filler")
VV_Filler = openmc.Material.mix_materials([tungsten_carbide, helium_8mpa, rafm_steel], [0.85, 0.10, 0.05], 'vo', name="VV Filler")
Magnet_Winding_Pack = openmc.Material.mix_materials([rebco, hastelloy, copper], [0.03, 0.47, 0.50], 'vo', name = "Winding Pack")

breeder_struct = 0.05
eutectic_breeding_material = openmc.Material.mix_materials([he_cooled_rafm, LiPb_breeder], [breeder_struct, 1-breeder_struct], 'vo', name= 'PbLi Breeder with Structure')

FW_FNSF = openmc.Material.mix_materials([helium_8mpa, rafm_steel], [0.6, 0.4], 'vo', name = "First Wall")
BW_FNSF = openmc.Material.mix_materials([helium_8mpa, rafm_steel], [0.2, 0.8], 'vo', name = "Back Wall")

he_manifold = openmc.Material.mix_materials([helium_8mpa, rafm_steel], [0.7, 0.3], 'vo', name = "Helium Manifold")

Magnet_Winding_Pack_2 = openmc.Material.mix_materials([copper, silver, rebco, lamno3, MgO, yttrium_oxide, alumina, hastelloy], [0.1312508, 0.05250033, 0.02625016, 0.000525, 0.00131251, 0.00013125, 0.000525, 0.7875049], 'vo', name='Magnet Winding Pack 2.0')

materials_list = [vacuum, air, deuterium, aluminum_6061, stainless, beryllium, lead, LiH, LiD,
                  rebco, magnet, tungsten, crispy, water, he_cooled_rafm, iron, lithium,
                  cooled_tungsten, tungsten_carbide, cooled_tungsten_carbide,
                  rafm_steel, LiPb_breeder, rings, tungsten_boride, WB2, w2b5, cooled_w2b5, cooled_wb2,
                  TiH2, cooled_TiH2, zirconium_hydride, water_cooled_wc, copper, hastelloy, flibe, tantalum,
                  tantalum_hydride_55, tantalum_hydride_30, cooled_rafm_steel, Nak_77, potassium, KCl,
                  HfH2, titanium, MgO, MgO_HfH2, Fe_HfH2_WB2, Ti_HfH2, B4C, enriched_lithium, lithium_metatitanate, lithium_orthosilicate, he_cooled_tungsten_carbide
                  ,VV_Filler, concrete, lamno3, alumina, yttrium_oxide, silver, HT_Shield_filler, LT_Shield_filler, Magnet_Winding_Pack, FW_FNSF, BW_FNSF, he_manifold, eutectic_breeding_material
                  ,cooled_stainless, Magnet_Winding_Pack_2, aluminum_1050, ss316ln, cooled_tungsten_boride, cooled_B4C]

for material in materials_list:
    material.depletable = True

materials = openmc.Materials(materials_list)


# %%
"""
fig=openmc.plot_xs(tungsten, ['capture'])
openmc.plot_xs(tungsten_carbide, ["capture"], axis=fig.get_axes()[0])
openmc.plot_xs(TiH2, ['capture'], axis=fig.get_axes()[0])
openmc.plot_xs(LiD, ["capture"], axis=fig.get_axes()[0])
openmc.plot_xs(LiH, ["capture"], axis=fig.get_axes()[0])
#openmc.plot_xs(tantalum_hydride, ["damage"], axis=fig.get_axes()[0])
plt.title("Neutron capture Macroscopic Cross Section")
plt.xlabel("Energy (eV)")
plt.xlim([1e5, 20e6])
plt.ylim([1e-5, 1e0])
plt.legend(["tungsten", "tungsten carbide", "TiH2", 'LiD', 'LiH'], loc='upper left')
#plt.xscale("linear")
#plt.ylim([1e-1, 5])
#plt.yscale("linear")
fig.savefig("./plots/W vs hydrides capture")
plt.show()
"""
# %%
