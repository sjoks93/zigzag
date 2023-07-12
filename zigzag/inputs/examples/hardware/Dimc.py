import os
from zigzag.classes.hardware.architecture.memory_hierarchy import MemoryHierarchy
from zigzag.classes.hardware.architecture.memory_level import MemoryLevel
from zigzag.classes.hardware.architecture.operational_unit import Multiplier
from zigzag.classes.hardware.architecture.operational_array import MultiplierArray
from zigzag.classes.hardware.architecture.memory_instance import MemoryInstance
from zigzag.classes.hardware.architecture.accelerator import Accelerator
from zigzag.classes.hardware.architecture.core import Core


def memory_hierarchy_dut(imc_array):
    breakpoint()
    """Memory hierarchy variables"""
    """ size=#bit, bw=(read bw, write bw), cost=(read word energy, write work energy) """
    cell_group = MemoryInstance(
        name="cell_group",
        size=imc_array.hd_param["weight_precision"]*imc_array.hd_param["group_depth"],
        r_bw=imc_array.hd_param["weight_precision"],
        w_bw=imc_array.hd_param["weight_precision"],
        r_cost=0,
        w_cost=imc_array.hd_param["w_energy_per_bit"] * 8, # unit: pJ/weight
        area=0, # this area is already included in imc_array
        r_port=0,
        w_port=0,
        rw_port=1,
        latency=0,
    )
    reg_I1 = MemoryInstance(
        name="rf_1B",
        size=8,
        r_bw=8,
        w_bw=8,
        r_cost=0,
        w_cost=0.013608,
        area=3.684*8/1e6,
        r_port=1,
        w_port=1,
        rw_port=0,
        latency=1,
    )

    reg_O1 = MemoryInstance(
        name="rf_1B",
        size=8,
        r_bw=8,
        w_bw=8,
        r_cost=0,
        w_cost=0.013608,
        area=3.684*8/1e6,
        r_port=2,
        w_port=2,
        rw_port=0,
        latency=1,
    )

    ##################################### on-chip memory hierarchy building blocks #####################################

    sram_256KB_256_3r_3w = MemoryInstance(
        name="sram_256KB",
        size=256 * 1024 * 8,
        r_bw=256,
        w_bw=256,
        r_cost=41.098622399999996,
        w_cost=27.001188,
        area=0.778311504,
        r_port=3,
        w_port=3,
        rw_port=0,
        latency=1,
        min_r_granularity=256,
        min_w_granularity=256,
    )

    #######################################################################################################################

    dram_100MB_32_3r_3w = MemoryInstance(
        name="dram_100MB",
        size=100*1024*1024*8,
        r_bw=32,
        w_bw=32,
        r_cost=3.7*32,
        w_cost=3.7*32,
        area=0,
        r_port=3,
        w_port=3,
        rw_port=0,
        latency=1,
    )

    memory_hierarchy_graph = MemoryHierarchy(operational_array=imc_array)

    """
    fh: from high = wr_in_by_high 
    fl: from low = wr_in_by_low 
    th: to high = rd_out_to_high
    tl: to low = rd_out_to_low
    """
    memory_hierarchy_graph.add_memory(
        memory_instance=cell_group,
        operands=("I2",),
        port_alloc=({"fh": "rw_port_1", "tl": "rw_port_1", "fl": None, "th": None},),
        served_dimensions=set(),
    )
    memory_hierarchy_graph.add_memory(
        memory_instance=reg_I1,
        operands=("I1",),
        port_alloc=({"fh": "w_port_1", "tl": "r_port_1", "fl": None, "th": None},),
        served_dimensions={(1, 0, 0)},
    )
    memory_hierarchy_graph.add_memory(
        memory_instance=reg_O1,
        operands=("O",),
        port_alloc=(
            {"fh": "w_port_1", "tl": "r_port_1", "fl": "w_port_2", "th": "r_port_2"},),
        served_dimensions={(0, 1, 0)},
    )

    ##################################### on-chip highest memory hierarchy initialization #####################################

    memory_hierarchy_graph.add_memory(
        memory_instance=sram_256KB_256_3r_3w,
        operands=("I1","O",),
        port_alloc=(
            {"fh": "w_port_1", "tl": "r_port_1", "fl": None, "th": None},
            {"fh": "w_port_2", "tl": "r_port_2", "fl": "w_port_3", "th": "r_port_3"},
        ),
        served_dimensions="all",
    )

    ####################################################################################################################

    memory_hierarchy_graph.add_memory(
        memory_instance=dram_100MB_32_2r_2w,
        operands=("I1", "I2", "O"),
        port_alloc=(
            {"fh": "w_port_1", "tl": "r_port_1", "fl": None, "th": None},
            {"fh": "w_port_2", "tl": "r_port_2", "fl": None, "th": None},
            {"fh": "w_port_1", "tl": "r_port_1", "fl": "w_port_3", "th": "r_port_3"},
        ),
        served_dimensions="all",
    )

    from zigzag.visualization.graph.memory_hierarchy import (
        visualize_memory_hierarchy_graph,
    )

    # visualize_memory_hierarchy_graph(memory_hierarchy_graph)
    return memory_hierarchy_graph


def imc_array_dut():
    """Multiplier array variables"""
    tech_param_28nm = {
        "tech_node": 0.028,             # unit: um
        "vdd":      0.9,                # unit: V
        "nd2_cap":  0.7/1e3,            # unit: pF
        "xor2_cap": 0.7*1.5/1e3,        # unit: pF
        "dff_cap":  0.7*3/1e3,          # unit: pF
        "nd2_area": 0.614/1e6,          # unit: mm^2
        "xor2_area":0.614*2.4/1e6,      # unit: mm^2
        "dff_area": 0.614*6/1e6,        # unit: mm^2
        "nd2_dly":  0.0478,             # unit: ns
        "xor2_dly": 0.0478*1.5,         # unit: ns
        "dff_dly":  0.0478*3.4,         # unit: ns
    }
    hd_param = {
        "imc_type":             "DIMC", # "DIMC" or "AIMC". Or else: pure digital
        "input_precision":      8,      # activation precison
        "weight_precision":     8,      # weight precision
        "input_bit_per_cycle":  1,      # nb_bits of input/cycle
        "group_depth":          1,      # m factor
        "w_energy_per_bit":     0.0101119640625,    # unit: pJ/bit (energy of weight writing)
    }

    dimensions = {
        "D1": 4,    # Do
        "D2": 32,   # Di (fix meaning: the dimension where addition happens)
        "D3": 1,    # nb_macros (fix meaning: nb_arrays)
    }  # {"D1": ("K", 4), "D2": ("C", 32),}

    dimc_array = DimcArray(
        tech_param_28nm, hd_param, dimensions
    )

    return dimc_array

def cores_dut():
    imc_array1 = imc_array_dut()
    memory_hierarchy1 = memory_hierarchy_dut(imc_array1)

    core1 = Core(1, imc_array1, memory_hierarchy1)

    return {core1}


cores = cores_dut()
acc_name = os.path.basename(__file__)[:-3]
accelerator = Accelerator(acc_name, cores)
