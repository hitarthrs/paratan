# src/paratan/runner.py

import openmc
import yaml
#from paratan.models.simple import build_simple_mirror
from src.paratan.models.tandem_model_builder import build_tandem_model_from_input

# def run_simple_model(input_file):
#     with open(input_file, "r") as f:
#         input_data = yaml.safe_load(f)

#     model = build_simple_mirror(input_data)
#     model.run()

def run_tandem_model(input_file, output_directory = 'output'):
    with open(input_file, "r") as f:
        input_data = yaml.safe_load(f)

    model = build_tandem_model_from_input(input_data, output_directory)
    model.run(cwd = output_directory)
