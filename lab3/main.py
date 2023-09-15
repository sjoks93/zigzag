import os
import sys
import argparse
import re

sys.path.insert(0, os.getcwd())

from zigzag.api import get_hardware_performance_zigzag
from zigzag.visualization.results.plot_cme import bar_plot_cost_model_evaluations_total

from inputs.hardware.c_k import accelerator as accelerator_c_k
from inputs.hardware.ox_k import accelerator as accelerator_ox_k
from inputs.hardware.ox_fx_fy import accelerator as accelerator_ox_fx_fy
from inputs.mapping.mapping_c_k import mapping as mapping_c_k
from inputs.mapping.mapping_ox_k import mapping as mapping_ox_k
from inputs.mapping.mapping_ox_fx_fy import mapping as mapping_ox_fx_fy
from zigzag.classes.stages import *

parser = argparse.ArgumentParser(description="Setup zigzag inputs")

parser.add_argument(
    "--model",
    metavar="path",
    required=True,
    help="path to onnx model, e.g. inputs/examples/my_onnx_model.onnx",
)
# Path to the workload onnx model
# onnx_model_path = "zigzag/inputs/examples/workload/resnet18.onnx"
args = parser.parse_args()
onnx_model_path = args.model

# List of hardware architectures we run our experiment for
accelerator1 = "lab3.inputs.hardware.c_k"
accelerator2 = "lab3.inputs.hardware.ox_fx_fy"
accelerator3 = "lab3.inputs.hardware.ox_k"
mapping1 = "lab3.inputs.mapping.mapping_c_k"
mapping2 = "lab3.inputs.mapping.mapping_ox_fx_fy"
mapping3 = "lab3.inputs.mapping.mapping_ox_k"
hardwares = [accelerator1, accelerator2, accelerator3]
# List of mappings for each hardware (encodes the spatial dataflow)
mappings = [mapping1, mapping2, mapping3]

# hardwares = [accelerator_c_k, accelerator_ox_fx_fy, accelerator_ox_k]
# mappings = [mapping_c_k, mapping_ox_fx_fy, mapping_ox_k]
cmes = []
for (hardware, mapping) in zip(hardwares, mappings):
    mainstage = MainStage(
        [  # Initializes the MainStage as entry point
            WorkloadParserStage,  # Parses the ONNX Model into the workload
            AcceleratorParserStage,  # Parses the accelerator
            CompleteSaveStage,  # Saves all received CMEs information to a json
            WorkloadStage,  # Iterates through the different layers in the workload
            SpatialMappingGeneratorStage,  # Generates multiple spatial mappings (SM)
            MinimalLatencyStage,  # Reduces all CMEs, returning minimal latency one
            TemporalOrderingConversionStage,  # Converts defined temporal_ordering to temporal mapping
            CostModelStage,  # Evaluates generated SM and TM through cost model
        ],
        accelerator=hardware,  # required by AcceleratorParserStage
        workload=onnx_model_path,  # required by ONNXModelParserStage
        mapping=mapping,  # required by ONNXModelParserStage
        dump_filename_pattern=f"lab1/outputs/3{hardware}-?.json",  # output file save pattern, ? will be replaced
        loma_lpf_limit=6,  # required by LomaStage
        loma_show_progress_bar=True,  # shows a progress bar while iterating over temporal mappings
    )   
    mainstage.run() 
#     # Pickle filename to save list of cmes
#     pickle_filename = "lab3/outputs/list_of_cmes.pickle"
#     # Call the zigzag api, using a provided accelerator and mapping
#     energy, latency, results = get_hardware_performance_zigzag(
#         onnx_model_path,
#         hardware,
#         mapping,
#         opt="latency",
#         dump_filename_pattern=f"lab3/outputs/{hardware.name}.json",
#         pickle_filename=pickle_filename,
#     )
#     cmes.append(results[0][0])

# x_labels = [hardware.name for hardware in hardwares]
# bar_plot_cost_model_evaluations_total(
#     cmes,
#     labels=x_labels,
#     save_path="lab3/outputs/plot_total.png",
# )
