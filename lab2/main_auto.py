import re
import argparse

from zigzag.classes.stages import *
from zigzag.visualization.results.print_mapping import print_mapping

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
)
parser.add_argument(
    "--accelerator",
    metavar="path",
    required=True,
    help="module path to the accelerator, e.g. inputs.examples.accelerator1",
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
        # WorkloadParserStage,  # Parses a user-defined workload
        ONNXModelParserStage,  # Parses the ONNX Model into the workload
        AcceleratorParserStage,  # Parses the accelerator
        SimpleSaveStage,  # Saves all received CMEs information to a json
        WorkloadStage,  # Iterates through the different layers in the workload
        SpatialMappingGeneratorStage,  # Generates multiple spatial mappings (SM)
        MinimalLatencyStage,  # Reduces all CMEs, returning minimal latency one
        LomaStage,  # Generates multiple temporal mappings (TM)
        CostModelStage,  # Evaluates generated SM and TM through cost model
    ],
    accelerator=args.accelerator,  # required by AcceleratorParserStage
    workload=args.model,  # required by ONNXModelParserStage
    mapping=args.mapping,  # required by ONNXModelParserStage
    dump_filename_pattern=f"lab2/outputs/{experiment_id}-?-auto.json",  # output file save pattern, ? will be replaced
    loma_lpf_limit=6,  # required by LomaStage
    loma_show_progress_bar=True,  # shows a progress bar while iterating over temporal mappings
)

# Launch the MainStage
results = mainstage.run()
cme = results[0][0]
print_mapping(cme)  # print the temporal mapping
