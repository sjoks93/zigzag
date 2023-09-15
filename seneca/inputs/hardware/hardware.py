from zigzag.classes.hardware.architecture.accelerator import Accelerator
from zigzag.classes.hardware.architecture.operational_unit import Multiplier
from zigzag.classes.hardware.architecture.operational_array import MultiplierArray
from zigzag.classes.hardware.architecture.memory_instance import MemoryInstance
from zigzag.classes.hardware.architecture.memory_instance import MemoryInstance
from zigzag.classes.hardware.architecture.memory_hierarchy import MemoryHierarchy
from zigzag.classes.hardware.architecture.core import Core


def get_multiplier_array():
    """Multiplier array variables"""
    multiplier_input_precision = [16, 16]
    multiplier_energy = 1.2
    multiplier_area = 0
    dimensions = {"D0":4, "D1": 4, "D2": 4, "D3": 8}

    multiplier = Multiplier(
        multiplier_input_precision, multiplier_energy, multiplier_area
    )
    multiplier_array = MultiplierArray(multiplier, dimensions)

    return multiplier_array


def get_memory_hierarchy(multiplier_array):

    """Memory hierarchy variables"""
    """ size=#bit, bw=#bit"""
    # Defintion of register file for inputs and weights
    rf_m = MemoryInstance(
        name="rf_m",
        mem_type="rf",
        size= 64 * 16,
        r_bw= 16,
        w_bw = 16,
        r_port=2,
        w_port=1,
        rw_port=1,
        r_cost = 0.2,
        w_cost = 0.2,
        latency=1,
        auto_cost_extraction=False,
    )
    # Defintion of first SRAM for inputs and outputs
    sram_m = MemoryInstance(
        name="sram_m",
        mem_type="sram",
        size= 8*1024 * 32,
        r_bw = 16,
        w_bw = 16,
        rw_port=2,
        r_cost = 3.5,
        w_cost = 3.7,
        latency= 2,
        auto_cost_extraction=False,
    )
    # Defintion of first SRAM for weights
    mram = MemoryInstance(
        name="MRAM",
        mem_type="sram",
        size= 1024*265*144,  # 36.8 MB
        r_bw=32,
        w_bw = 32,
        rw_port=2,
        latency=25,
        r_cost = 32,
        w_cost = 32,
        auto_cost_extraction=False,
    )

    dram = MemoryInstance(
        name="dram",
        mem_type="dram",
        size=1073741824 * 8,
        r_bw= 32,
        w_bw= 32,
        rw_port = 2,
        latency= 135,
        r_cost = 112,
        w_cost = 112,
        auto_cost_extraction=False,
    )

    memory_hierarchy_graph = MemoryHierarchy(operational_array=multiplier_array)

    """
    fh: from high = wr_in_by_high = 
    fl: from low = wr_in_by_low 
    th: to high = rd_out_to_high = 
    tl: to low = rd_out_to_low = 
    """
    # Register file for weight
    memory_hierarchy_graph.add_memory(
        memory_instance=rf_m,
        operands=("I1", "I2", "O"),
        port_alloc=({"fh": "rw_port_1", "th": "rw_port_1", "fl": "w_port_1", "tl": "r_port_1", "tl": "r_port_2"},
                    {"fh": "rw_port_1", "th": "rw_port_1", "fl": "w_port_1", "tl": "r_port_1", "tl": "r_port_2"},
                    {"fh": "rw_port_1", "th": "rw_port_1", "fl": "w_port_1", "tl": "r_port_1", "tl": "r_port_2"},),
        served_dimensions=set(),
    )

    # First SRAM for weights
    memory_hierarchy_graph.add_memory(
        memory_instance=sram_m,
        operands=("I1", "I2", "O"),
        port_alloc=({"fh": "rw_port_1", "tl": "rw_port_1", "fl": "rw_port_2", "th": "rw_port_2"},
                    {"fh": "rw_port_1", "tl": "rw_port_1", "fl": "rw_port_2", "th": "rw_port_2"},
                    {"fh": "rw_port_1", "tl": "rw_port_1", "fl": "rw_port_2", "th": "rw_port_2"},),
        served_dimensions=set(),
    )

    # First SRAM for inputs and outputs
    memory_hierarchy_graph.add_memory(
        memory_instance=mram,
        operands=("I1", "I2", "O"),
        port_alloc=({"fh": "rw_port_1", "tl": "rw_port_1", "fl": "rw_port_2", "th": "rw_port_2"},
                    {"fh": "rw_port_1", "tl": "rw_port_1", "fl": "rw_port_2", "th": "rw_port_2"},
                    {"fh": "rw_port_1", "tl": "rw_port_1", "fl": "rw_port_2", "th": "rw_port_2"},),

        served_dimensions={(0, 1, 1, 1)},
    )

    memory_hierarchy_graph.add_memory(
        memory_instance=dram,
        operands=("I1", "I2", "O"),
        port_alloc=(
            {"fh": "rw_port_1", "tl": "rw_port_2", "fl": None, "th": None},
            {"fh": "rw_port_1", "tl": "rw_port_2", "fl": None, "th": None},
            {
                "fh": "rw_port_1",
                "tl": "rw_port_2",
                "fl": "rw_port_1",
                "th": "rw_port_2",
            },
        ),
        served_dimensions="all",
    )

    return memory_hierarchy_graph


# def get_dataflows():
#     return [
#         {"D1": ("C", 32), "D2": ("K", 32)},
#         {"D1": ("G", 32)},
#     ]


def get_core(id):
    operational_array = get_multiplier_array()
    memory_hierarchy = get_memory_hierarchy(operational_array)
    core = Core(id, operational_array, memory_hierarchy)
    return core


cores = {get_core(id=0)}
name = "accelerator-seneca"
accelerator = Accelerator(name, cores)
