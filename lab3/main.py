import onnx
from zigzag.api import get_hardware_performance_zigzag
from zigzag.inputs.examples.mapping.tpu_like import mapping as mapping_tpu
from zigzag.inputs.examples.mapping.edge_tpu_like import mapping as mapping_edge_tpu
from zigzag.inputs.examples.mapping.tesla_npu_like import mapping as mapping_tesla_npu
from zigzag.inputs.examples.mapping.meta_prototype_like import (
    mapping as mapping_meta_prototype,
)

from zigzag.visualization.results.plot_cme import bar_plot_cost_model_evaluations_total

# Path to the workload onnx model
onnx_model_path = "lab1/resnet18_first_layer.onnx"
# OR directly load in the onnx model
workload_model = onnx.load(onnx_model_path, load_external_data=False)

# List of hardware architectures we run our experiment for
hardwares = ["tpu", "edge-tpu", "tesla-npu", "meta-prototype"]
# List of mappings for each hardware (encodes the spatial dataflow)
mappings = [mapping_tpu, mapping_edge_tpu, mapping_tesla_npu, mapping_meta_prototype]

cmes = []
for (hardware, mapping) in zip(hardwares, mappings):
    # Pickle filename to save list of cmes
    pickle_filename = "lab3/outputs/list_of_cmes.pickle"
    # Call the zigzag api, using a provided accelerator and mapping
    energy, latency, results = get_hardware_performance_zigzag(
        onnx_model_path,
        hardware,
        mapping=mapping,
        opt="latency",
        pickle_filename=pickle_filename,
    )
    print(
        f"Total onnx model (energy, latency) performance = ({energy:.3e}, {latency:.3e})."
    )
    cmes.append(results[0][0])


bar_plot_cost_model_evaluations_total(
    cmes,
    labels=hardwares,
    fig_title="Best latency mapping found",
    save_path="lab3/outputs/plot_total.png",
)
