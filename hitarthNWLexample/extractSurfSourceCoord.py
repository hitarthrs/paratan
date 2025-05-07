#!/openmc_venv/bin/python
#
#
import h5py
import numpy as np
import sys
from scipy.optimize import direct
#import openmc # do not need this
#
#print("Everything imported now")
#
def extract_coords(source_file):
    """Extracts Cartesian coordinates of particle surface crossings given an
    OpenMC surface source file.

    Arguments:
        source_file (str): path to OpenMC surface source file.

    Returns:
        coords (array of array of float): Cartesian coordinates of all particle
            surface crossings.
    """
    # Load source file
    file = h5py.File(source_file, "r")
    # Extract source information
    dataset = file["source_bank"]["r"]
#    dataset2= file["source_bank"]["E"]
    # Construct matrix of particle crossing coordinates
    coords = np.empty((len(dataset), 3))    
#    coords = np.empty((len(dataset), 4))
    coords[:, 0] = dataset["x"]
    coords[:, 1] = dataset["y"]
    coords[:, 2] = dataset["z"]
#    coords[:, 3] = dataset2 # for particle energy

    return coords
#
#
source_file="surface_source.h5"
#
coords = extract_coords(source_file)
#
#print("Printing coords array data type")
#print(coords.dtype)
#print("Printing coords array shape")
#print(coords.shape)
#print(coords)
# write to csv file (doesn't work the way I like)
#print("Writing to CSV file:")
#coords.tofile('surface_source.csv', sep = ',')
#
# just loop over array to print (works but contains [] brackets)
#for xyz in coords:
#    print(xyz)
#
np.savetxt(sys.stdout, coords, delimiter=" ", fmt="%.4e")
#
#print("All done")
