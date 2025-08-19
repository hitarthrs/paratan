# src/paratan/runner.py

import openmc
import yaml
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
#from paratan.models.simple import build_simple_mirror
from src.paratan.models.tandem_model_builder import build_tandem_model_from_input
from src.paratan.models.simple_model_builder import build_simple_model_from_input
import src.paratan.materials.material as m

# def run_simple_model(input_file):
#     with open(input_file, "r") as f:
#         input_data = yaml.safe_load(f)

#     model = build_simple_mirror(input_data)
#     model.run()

def build_simple_model_only(input_file, output_directory='output'):
    """Build the simple mirror model and create XML files/plots, but don't run the simulation."""
    with open(input_file, "r") as f:
        input_data = yaml.safe_load(f)

    model = build_simple_model_from_input(input_data, output_directory)
    # Don't call model.run() - just return the built model
    return model

def run_simple_model_modular(input_file, output_directory='output'):
    """Run the simple mirror model using the modular builder."""
    with open(input_file, "r") as f:
        input_data = yaml.safe_load(f)

    model = build_simple_model_from_input(input_data, output_directory)
    model.run(cwd=output_directory)

def run_tandem_model(input_file, output_directory = 'output'):
    """Run the tandem model using the modular builder."""
    with open(input_file, "r") as f:
        input_data = yaml.safe_load(f)

    model = build_tandem_model_from_input(input_data, output_directory)
    model.run(cwd = output_directory)
