from paratan2.config.models import SimpleMachineConfig
from paratan2.runner import build

config = SimpleMachineConfig.from_yaml(
    "/home/hrshah3/paratan/paratan/compare_hf_study_runs/hpc_run/simple_parametric_input_new.yaml"
)
model = build(config, output_dir="/tmp/paratan2_test")
