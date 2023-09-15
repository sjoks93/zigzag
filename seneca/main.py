import os
import sys
import argparse
import re

sys.path.insert(0, os.getcwd())
from zigzag.classes.stages import *


from zigzag.visualization.results.plot_cme import (
    bar_plot_cost_model_evaluations_breakdown,
)

# Get the onnx model, the mapping and accelerator arguments
parser = argparse.ArgumentParser(description="Setup zigzag inputs")
parser.add_argument(
    "--model",
    metavar="path",
    required=True,
    help="path to onnx model, e.g. inputs/examples/my_onnx_model.onnx",
)
parser.add_argument(
    "--mapping",
    metavar="path",
    required=True,
    help="path to mapping file, e.g., inputs.examples.my_mapping",
    default= "lab3.inputs.mapping.mapping_seneca",
)
parser.add_argument(
    "--accelerator",
    metavar="path",
    required=True,
    help="module path to the accelerator, e.g. inputs.examples.accelerator1",
    default="lab3.inputs.hardware.seneca"
)
args = parser.parse_args()

# Initialize the logger
import logging as _logging

_logging_level = _logging.INFO
# _logging_format = "%(asctime)s - %(funcName)s +%(lineno)s - %(levelname)s - %(message)s"
_logging_format = "%(asctime)s - %(levelname)s - %(message)s"
_logging.basicConfig(level=_logging_level, format=_logging_format)

hw_name = args.accelerator.split(".")[-1]
wl_name = re.split(r"/|\.", args.model)[-1]
if wl_name == "onnx":
    wl_name = re.split(r"/|\.", args.model)[-2]
experiment_id = f"{hw_name}-{wl_name}"

# Initialize the MainStage which will start execution.
# The first argument of this init is the list of stages that will be executed in sequence.
# The second argument of this init are the arguments required for these different stages.
mainstage = MainStage(
    [  # Initializes the MainStage as entry point
        WorkloadParserStage, #File workload parser
        #ONNXModelParserStage,  # Parses the ONNX Model into the workload
        AcceleratorParserStage,  # Parses the accelerator
        CompleteSaveStage,  # Saves all received CMEs information to a json
        WorkloadStage,  # Iterates through the different layers in the workload
        SpatialMappingGeneratorStage,  # Generates multiple spatial mappings (SM)
        MinimalLatencyStage,  # Reduces all CMEs, returning minimal latency one
        LomaStage,  # Converts defined temporal_ordering to temporal mapping
        CostModelStage,  # Evaluates generated SM and TM through cost model
    ],
    accelerator=args.accelerator,  # required by AcceleratorParserStage
    workload=args.model,  # required by ONNXModelParserStage
    mapping=args.mapping,  # required by ONNXModelParserStage
    dump_filename_pattern=f"lab3/outputs/seneca_alexnet/{experiment_id}-?.json",  # output file save pattern, ? will be replaced
    loma_lpf_limit=6,  # required by LomaStage
    loma_show_progress_bar=True,  # shows a progress bar while iterating over temporal mappings
)
# mainstage = MainStage([  # Initializes the MainStage as entry point
#     ONNXModelParserStage,  # Parses the ONNX Model into the workload
#     AcceleratorParserStage,  # Parses the accelerator
#     SimpleSaveStage,  # Saves all received CMEs information to a json
#     WorkloadStage,  # Iterates through the different layers in the workload
#     SpatialMappingGeneratorStage,  # Generates multiple spatial mappings (SM)
#     MinimalLatencyStage,  # Reduces all CMEs, returning minimal latency one
#     SalsaStage,  # Find pseudo-optimal temporal mapping
#     CostModelStage  # Evaluates generated SM and TM through cost model
# ],
#     accelerator=args.accelerator,  # required by AcceleratorParserStage
#     workload=args.model,  # required by ONNXModelParserStage
#     mapping=args.mapping,  # required by ONNXModelParserStage
#     dump_filename_pattern=f"outputs/{experiment_id}-layer_?.json",  # output file save pattern
#     loma_lpf_limit=6,  # required by LomaStage
#     loma_show_progress_bar=True,  # shows a progress bar while iterating over temporal mappings
#     salsa_iteration_number=1000,
#     salsa_start_temperature=0.05,
#     salsa_opt_criterion="latency",
#     salsa_number_of_core=8
# )

# Launch the MainStage
answers = mainstage.run()
# Plot the energy and latency breakdown of our cost model evaluation
# cme = answers[0][0]
# save_path = "lab3/outputs/seneca_resnet/breakdown"

# from zigzag.visualization.results.print_mapping import print_mapping
# for i in range (len(answers)):
#     print_mapping(answers[i][0])
#     bar_plot_cost_model_evaluations_breakdown([answers[i][0]], save_path=save_path+str(i)+".png", xtick_rotation=0)

