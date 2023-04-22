import os
import sys

sys.path.insert(0, os.getcwd())

from zigzag.api import get_hardware_performance_zigzag
from zigzag.visualization.results.plot_cme import bar_plot_cost_model_evaluations_total

from inputs.hardware.c_k import accelerator as accelerator_c_k
from inputs.hardware.ox_k import accelerator as accelerator_ox_k
from inputs.hardware.ox_fx_fy import accelerator as accelerator_ox_fx_fy
from inputs.mapping.mapping_c_k import mapping as mapping_c_k
from inputs.mapping.mapping_ox_k import mapping as mapping_ox_k
from inputs.mapping.mapping_ox_fx_fy import mapping as mapping_ox_fx_fy


# Path to the workload onnx model
# onnx_model_path = "zigzag/inputs/examples/workload/resnet18.onnx"
onnx_model_path = "lab1/resnet18_first_layer.onnx"

# List of hardware architectures we run our experiment for
hardwares = [accelerator_c_k, accelerator_ox_k, accelerator_ox_fx_fy]
# List of mappings for each hardware (encodes the spatial dataflow)
mappings = [mapping_c_k, mapping_ox_k, mapping_ox_fx_fy]

cmes = []
for (hardware, mapping) in zip(hardwares, mappings):
    # Pickle filename to save list of cmes
    pickle_filename = "lab3/outputs/list_of_cmes.pickle"
    # Call the zigzag api, using a provided accelerator and mapping
    energy, latency, results = get_hardware_performance_zigzag(
        onnx_model_path,
        hardware,
        mapping,
        opt="latency",
        dump_filename_pattern=f"lab3/outputs/{hardware.name}.json",
        pickle_filename=pickle_filename,
    )
    cmes.append(results[0][0])

x_labels = [hardware.name for hardware in hardwares]
bar_plot_cost_model_evaluations_total(
    cmes,
    labels=x_labels,
    save_path="lab3/outputs/plot_total.png",
)
