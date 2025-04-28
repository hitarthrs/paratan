import openmc
import numpy as np
import matplotlib.pyplot as plt
import material as m
from geometry_lib import *
from source_lib import *
import yaml


with open('input_files/tandem_parametric_input.yaml', 'r') as f:
    input_data = yaml.safe_load(f)


################### Universe that contains the machine #####################

universe_machine = openmc.Universe(786)
