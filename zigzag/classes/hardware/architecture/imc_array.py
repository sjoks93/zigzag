import numpy as np
import math
from zigzag.classes.hardware.architecture.dimension import Dimension
from typing import Dict
from zigzag.classes.hardware.architecture.get_cacti_cost import get_cacti_cost

###############################################################################################################
# README
# This file includes:
#   . class LogicUnit (defines the energy/area/delay cost of multipliers, adders, regs)
#   . class ImcArray (provides initialization function, function of calling CACTI for DimcArray and AimcArray)
#   . class DimcArray (defines the energy/area/delay cost of a DIMC array)
#   . class AimcArray (defines the energy/area/delay cost of an ADC, a DAC and an AIMC array)
# How to use this file?
#   . This file is internally called in ZigZag-IMC framework.
#   . It can also be run independently, for mainly debugging. An example is given at the end of the file.
###############################################################################################################


class LogicUnit:
    """cost (energy, area, delay) of 1b adder, 1b multiplier, 1b register is defined in this class"""
    def __init__(self, tech_param:dict):
        """
        Input example:
        tech_param_28nm = {
        "vdd":      0.9,            # unit: V
        "nd2_cap":  0.7/1e3,        # unit: pF
        "nd2_area": 0.614/1e6,      # unit: mm^2
        "nd2_dly":  0.0478,         # unit: ns
        "xor2_cap": 0.7*1.5/1e3,    # unit: pF
        "xor2_area":0.614*2.4/1e6,  # unit: mm^2
        "xor2_dly": 0.0478*1.5,     # unit: ns
        "dff_cap":  0.7*3/1e3,      # unit: pF
        "dff_area": 0.0614*6/1e6,   # unit: mm^2
        "dff_dly":  0.0478*3.4,     # unit: ns
    }
        """
        """check input firstly"""
        required_param = ["tech_node", "vdd", "nd2_cap", "nd2_area", "nd2_dly", "xor2_cap", "xor2_area", "xor2_dly", "dff_cap", "dff_area", "dff_dly"]
        for ii_a, a in enumerate(required_param):
            if a not in tech_param.keys():
                raise Exception(f"[LogicUnit] Incorrect input, required param [{a}] not found.")
            if not (isinstance(tech_param[a], int) or isinstance(tech_param[a], float)):
                raise Exception(f"[LogicUnit] Incorrect input, value [{tech_param[a]}] of param [{a}] is not a num.")
            if tech_param[a] <= 0:
                raise Exception(f"[LogicUnit] Incorrect input, value [{tech_param[a]}] of param [{a}] is not positive.")
        """initialization"""
        self.tech_param = tech_param

    def get_1b_adder_energy(self):
        """energy of 1b full adder"""
        """Assume a 1b adder has 3 ND2 gate and 2 XOR2 gate"""
        adder_cap =  3 * self.tech_param["nd2_cap"] + 2 *  self.tech_param["xor2_cap"]
        return adder_cap * (self.tech_param["vdd"]**2) # unit: pJ

    def get_1b_adder_energy_half_activated(self):
        """energy of 1b full adder when 1 input is 0"""
        adder_cap = 2 * self.tech_param["xor2_cap"]
        return adder_cap * (self.tech_param["vdd"] ** 2)  # unit: pJ

    def get_1b_multiplier_energy(self):
        """energy of 1b multiplier"""
        """1b mult includes 1 NOR gate, which is assumed as the same cost of ND2 gate"""
        """why 0.5: considering weight stays constant during multiplication"""
        return 0.5 * self.tech_param["nd2_cap"] * (self.tech_param["vdd"] ** 2) # unit: pJ

    def get_1b_reg_energy(self):
        """energy of 1b DFF"""
        return self.tech_param["dff_cap"] * (self.tech_param["vdd"] ** 2) # unit: pJ

    def get_1b_adder_area(self):
        """area of 1b full adder"""
        """Assume a 1b adder has 3 ND2 gate and 2 XOR2 gate"""
        adder_area = 3 * self.tech_param["nd2_area"] + 2 * self.tech_param["xor2_area"]
        return adder_area

    def get_1b_multiplier_area(self):
        """area of 1b multiplier"""
        """1b mult includes 1 NOR gate, which is assumed as the same cost of ND2 gate"""
        return self.tech_param["nd2_area"]

    def get_1b_reg_area(self):
        """area of 1b DFF"""
        return self.tech_param["dff_area"]

    def get_1b_adder_dly_in2sum(self):
        """delay of 1b adder: input to sum-out"""
        adder_dly = 2 * self.tech_param["xor2_dly"]
        return adder_dly

    def get_1b_adder_dly_in2cout(self):
        """delay of 1b adder: input to carry-out"""
        adder_dly = self.tech_param["xor2_dly"] + 2 * self.tech_param["nd2_dly"]
        return adder_dly

    def get_1b_adder_dly_cin2cout(self):
        """delay of 1b adder: carry-in to carry-out"""
        adder_dly = 2 * self.tech_param["nd2_dly"]
        return adder_dly

    def get_1b_multiplier_dly(self):
        """delay of 1b multiplier"""
        """1b mult includes 1 NOR gate, which is assumed as the same cost of ND2 gate"""
        return self.tech_param["nd2_dly"]

    def get_1b_reg_dly(self):
        """delay of 1b DFF"""
        """why 0? Compared to others, it's negligible"""
        return 0

class Imc:
    """definition of general initilization function for D/AIMC"""
    def __init__(self,tech_param:dict, hd_param:dict, dimensions:dict):
        """check input firstly"""
        required_hd_param = ["imc_type", "input_precision", "weight_precision", "input_bit_per_cycle", "group_depth", "w_energy_per_bit"]
        required_dimension_param = ["D1", "D2", "D3"]
        for ii_a, a in enumerate(required_hd_param):
            if a not in hd_param.keys():
                raise Exception(f"[ImcArray] Incorrect hd_param, required param [{a}] not found.")
            if a == "imc_type":
                if hd_param[a] not in ["AIMC", "DIMC"]:
                    raise Exception(f"[ImcArray] Incorrect imc_type in hd_param, [AIMC] or [DIMC] is expected.")
            else:
                if not (isinstance(hd_param[a], int) or isinstance(hd_param[a], float)):
                    raise Exception(f"[ImcArray] Incorrect hd_param, value [{hd_param[a]}] of param [{a}] is not a num.")
                if hd_param[a] <= 0:
                    raise Exception(f"[ImcArray] Incorrect hd_param, value [{hd_param[a]}] of param [{a}] is not positive.")
                if a == "input_bit_per_cycle" and hd_param[a] > hd_param["input_precision"]:
                    input_precision = hd_param["input_precision"]
                    raise Exception(f"[ImcArray] Incorrect hd_param, value [{hd_param[a]}] of param [{a}] is bigger than [input_precision] ({input_precision}).")
        for ii_a, a in enumerate(required_dimension_param):
            if a not in dimensions.keys():
                raise Exception(f"[ImcArray] Incorrect dimensions, required dimension [{a}] not found.")
            if not isinstance(dimensions[a], int) or dimensions[a] <= 0:
                raise Exception(f"[ImcArray] Incorrect dimensions, value [{dimensions[a]}] of param [{a}] is not a positive int.")
        if hd_param["imc_type"] == "AIMC":
            a = "adc_resolution"
            if a not in hd_param.keys():
                raise Exception(f"[ImcArray] Incorrect hd_param, required param [{a}] not found.")
            if not (isinstance(hd_param[a], int) or isinstance(hd_param[a], float)) or hd_param[a] <= 0:
                raise Exception(f"[ImcArray] Incorrect hd_param, value [{hd_param[a]}] of param [{a}] is not a positive num.")
        """initialization"""
        self.hd_param = hd_param
        self.dimensions = dimensions
        # tech_param will be checked in LogicUnit class
        self.logic_unit = LogicUnit(tech_param)
        self.tech_param = tech_param

    def get_cell_array_cost(self):
        """get the area, energy cost of a single cell array using CACTI"""
        array_rows = self.dimensions["D2"] * self.hd_param["group_depth"]
        array_cols = self.dimensions["D1"] * self.hd_param["weight_precision"]
        cell_array_size = array_rows * array_cols / 8 # unit: byte
        imc_bw = self.dimensions["D1"] * self.hd_param["weight_precision"] # unit: bit
        if __name__ == "__main__":
            cacti_path = "../../cacti/cacti_master"
        else:
            cacti_path = "zigzag/classes/cacti/cacti_master"
        access_time, area, r_cost, w_cost = get_cacti_cost(cacti_path=cacti_path, tech_node=self.tech_param["tech_node"], mem_type="sram", mem_size_in_byte=cell_array_size, bw=imc_bw) # unit: ns, mm^2, nJ/access, nJ/access
        return access_time, area, r_cost, w_cost

class DimcArray(Imc):
    """definition of a DIMC array"""
    """
    constraint:
        -- activation precision must be in the power of 2.
        -- input_bit_per_cycle must be in the power of 2.
        -- 
    assumption:
    """
    def __init__(self,tech_param:dict, hd_param:dict, dimensions:dict):
        super().__init__(tech_param, hd_param, dimensions)

    def __jsonrepr__(self):
        """
        JSON Representation of this class to save it to a json file.
        """
        # not implemented
        #return {"operational_unit": self.unit, "dimensions": self.dimensions}
        pass

    def get_area(self):
        """area of cell array (cells, mults, adders, adders_pv, accumulators. Not include input/output regs)"""
        self.hd_param["cell_array_area"] = self.get_cell_array_cost()[1]
        area_cells = self.hd_param["cell_array_area"] * self.dimensions["D3"]
        """area of multiplier array"""
        area_mults = self.logic_unit.get_1b_multiplier_area() * self.hd_param["input_bit_per_cycle"] * np.prod(list(self.dimensions.values())) * self.hd_param["weight_precision"]
        """area of adder trees (type: RCA)"""
        adder_input_pres = self.hd_param["weight_precision"]
        nb_inputs_of_adder = self.dimensions["D2"]
        adder_depth = math.log2(nb_inputs_of_adder)
        assert adder_depth%1==0, f"[DimcArray] The number of inputs [{nb_inputs_of_adder}] for the adder tree is not in the power of 2."
        adder_depth = int(adder_depth) # float -> int for simplicity
        adder_output_pres = adder_input_pres + adder_depth
        nb_of_1b_adder = nb_inputs_of_adder * (adder_input_pres+1) - (adder_input_pres+adder_depth+1) # nb of 1b adders in a single adder tree
        nb_of_1b_adder *= self.hd_param["input_bit_per_cycle"] # multiply with nb_of_adder_trees along D2
        nb_of_1b_adder *= self.dimensions["D1"] * self.dimensions["D3"] # multiply with nb_of_adder_trees
        area_adders = self.logic_unit.get_1b_adder_area() * nb_of_1b_adder
        """area of adders with place values when input_bit_per_cycle>1 (type: RCA)"""
        nb_inputs_of_adder_pv = self.hd_param["input_bit_per_cycle"]
        if nb_inputs_of_adder_pv == 1:
            nb_of_1b_adder_pv = 0
        else:
            adder_depth_pv = math.log2(nb_inputs_of_adder_pv)
            assert adder_depth_pv%1==0, f"[DimcArray] The value [{nb_inputs_of_adder_pv}] of [input_bit_per_cycle] is not in the power of 2."
            adder_depth_pv = int(adder_depth_pv) # float -> int for simplicity
            nb_of_1b_adder_pv = adder_output_pres * (nb_inputs_of_adder_pv-1) + nb_inputs_of_adder_pv * (adder_depth_pv-0.5) # nb of 1b adders in a single place-value adder tree
            nb_of_1b_adder_pv *= self.dimensions["D1"] * self.dimensions["D3"] # multiply with nb_of_adder_trees
        area_adders_pv = self.logic_unit.get_1b_adder_area() * nb_of_1b_adder_pv
        """area of accumulators (adder type: RCA)"""
        if self.hd_param["input_bit_per_cycle"] == self.hd_param["input_precision"]:
            area_accumulators = 0
        else:
            accumulator_output_pres = self.hd_param["input_precision"]+self.hd_param["weight_precision"]+math.log2(self.dimensions["D2"])
            nb_of_1b_adder_accumulator = accumulator_output_pres * self.dimensions["D1"] * self.dimensions["D3"]
            nb_of_1b_reg_accumulator = nb_of_1b_adder_accumulator # number of regs in an accumulator
            area_accumulators = self.logic_unit.get_1b_adder_area() * nb_of_1b_adder_accumulator + self.logic_unit.get_1b_reg_area() * nb_of_1b_reg_accumulator
        """total area of imc"""
        self.area = area_cells + area_mults + area_adders + area_adders_pv + area_accumulators
        self.area_breakdown = { # unit: same with in input hd file
            "cells":    area_cells,
            "mults":    area_mults,
            "adders":   area_adders,
            "adders_pv":area_adders_pv,
            "accumulators": area_accumulators
        }
        return self.area

    def get_delay(self):
        """delay of imc (worst path: mults -> adders -> adders_pv -> accumulators) """
        dly_mults = self.logic_unit.get_1b_multiplier_dly()
        """delay of adders (tree) (type: RCA)"""
        adder_input_pres = self.hd_param["weight_precision"]
        nb_inputs_of_adder = self.dimensions["D2"]
        adder_depth = math.log2(nb_inputs_of_adder)
        adder_depth = int(adder_depth)  # float -> int for simplicity
        adder_output_pres = adder_input_pres + adder_depth
        dly_adders = (adder_depth-1) * self.logic_unit.get_1b_adder_dly_in2sum() + self.logic_unit.get_1b_adder_dly_in2cout() + (adder_output_pres-1-1) * self.logic_unit.get_1b_adder_dly_cin2cout()
        """delay of adders_pv (type: RCA)"""
        nb_inputs_of_adder_pv = self.hd_param["input_bit_per_cycle"]
        if nb_inputs_of_adder_pv == 1:
            dly_adders_pv = 0
            accumulator_input_pres = adder_output_pres
        else:
            adder_depth_pv = math.log2(nb_inputs_of_adder_pv)
            adder_depth_pv = int(adder_depth_pv)  # float -> int for simplicity
            adder_pv_input_precision = adder_output_pres
            adder_pv_output_precision = nb_inputs_of_adder_pv + adder_output_pres  # output precision from adders_pv (depth + input_precision)
            accumulator_input_pres = adder_pv_output_precision
            dly_adders_pv = (adder_depth_pv - 1) * self.logic_unit.get_1b_adder_dly_in2sum() + self.logic_unit.get_1b_adder_dly_in2cout() + (adder_pv_output_precision - adder_pv_input_precision-1) * self.logic_unit.get_1b_adder_dly_cin2cout()
        """delay of accumulators (adder type: RCA)"""
        accumulator_output_pres = self.hd_param["input_precision"] + self.hd_param["weight_precision"] + math.log2(self.dimensions["D2"])
        accumulator_output_pres = int(accumulator_output_pres) # float -> int for simplicity
        if accumulator_output_pres == accumulator_input_pres: # no accumulator
            dly_accumulators = 0
        else:
            dly_accumulators = self.logic_unit.get_1b_adder_dly_in2cout() + (accumulator_output_pres - accumulator_input_pres) * self.logic_unit.get_1b_adder_dly_cin2cout()
        """total delay of imc"""
        self.delay = dly_mults + dly_adders + dly_adders_pv + dly_accumulators
        self.delay_breakdown = {
            "mults":    dly_mults,
            "adders":   dly_adders,
            "adders_pv":dly_adders_pv,
            "accumulators": dly_accumulators
        }
        return self.delay

    def get_peak_performance(self, level="macro"):
        pass

    def get_energy(self, layer, mapping):
        pass

class AimcArray(Imc):
    def __init__(self,tech_param:dict, hd_param:dict, dimensions:dict):
        super().__init__(tech_param, hd_param, dimensions)

    def __jsonrepr__(self):
        """
        JSON Representation of this class to save it to a json file.
        """
        # not implemented
        # return {"operational_unit": self.unit, "dimensions": self.dimensions}
        pass

    def get_adc_cost(self):
        """single ADC cost calculation"""
        """area (mm^2)"""
        if self.hd_param["adc_resolution"] == 1:
            adc_area = 0
        else: # formula extracted and validated against 3 AIMC papers on 28nm
            k1 = -0.0369
            k2 = 1.206
            adc_area = 10**(k1*self.hd_param["adc_resolution"]+k2) * 2**self.hd_param["adc_resolution"] * (10**-6) # unit: mm^2
        """delay (ns)"""
        k3 = 0.00653 # ns
        k4 = 0.640 # ns
        adc_delay = self.hd_param["adc_resolution"] * (k3*self.dimensions["D2"] + k4) # unit: ns
        """energy (fJ)"""
        k5 = 100 # fF
        k6 = 0.001 # fF
        adc_energy = (k5 * self.hd_param["adc_resolution"] + k6 * 4**self.hd_param["adc_resolution"]) * self.tech_param["vdd"]**2 # unit: fJ
        return adc_area, adc_delay, adc_energy

    def get_dac_cost(self):
        """single DAC cost calculation"""
        """area (mm^2)"""
        dac_area = 0 # neglected
        """delay (ns)"""
        dac_delay = 0 # neglected
        """energy (fJ)"""
        if self.hd_param["input_bit_per_cycle"] == 1:
            dac_energy = 0
        else:
            k0 = 50 # fF
            dac_energy = k0 * self.hd_param["input_bit_per_cycle"] * self.tech_param["vdd"]**2 # unit: fJ
        return dac_area, dac_delay, dac_energy

    def get_area(self):
        """area of cell array (cells, mults, adders, adders_pv, accumulators. Not include input/output regs)"""
        self.hd_param["cell_array_area"] = self.get_cell_array_cost()[1]
        area_cells = self.hd_param["cell_array_area"] * self.dimensions["D3"]
        """area of multiplier array"""
        area_mults = self.logic_unit.get_1b_multiplier_area() * np.prod(list(self.dimensions.values())) * self.hd_param["weight_precision"]
        """area of ADCs"""
        area_adcs = self.get_adc_cost()[0] * self.hd_param["weight_precision"] * self.dimensions["D1"] * self.dimensions["D3"] # one ADC along D2
        """area of DACs"""
        area_dacs = self.get_dac_cost()[0] * self.dimensions["D2"] * self.dimensions["D3"] # one ADC along D1
        """area of adders with place values after ADC conversion (type: RCA)"""
        nb_inputs_of_adder_pv = self.hd_param["weight_precision"]
        if nb_inputs_of_adder_pv == 1:
            nb_of_1b_adder_pv = 0
        else:
            adder_depth_pv = math.log2(nb_inputs_of_adder_pv)
            assert adder_depth_pv % 1 == 0, f"[AimcArray] The value [{nb_inputs_of_adder_pv}] of [weight_precision] is not in the power of 2."
            adder_depth_pv = int(adder_depth_pv)  # float -> int for simplicity
            adder_input_precision = self.hd_param["adc_resolution"]
            nb_of_1b_adder_pv = adder_input_precision * (nb_inputs_of_adder_pv - 1) + nb_inputs_of_adder_pv * (adder_depth_pv - 0.5)  # nb of 1b adders in a single place-value adder tree
            nb_of_1b_adder_pv *= self.dimensions["D1"] * self.dimensions["D3"]  # multiply with nb_of_adder_trees
        area_adders_pv = self.logic_unit.get_1b_adder_area() * nb_of_1b_adder_pv
        """area of accumulators (adder type: RCA)"""
        if self.hd_param["input_bit_per_cycle"] == self.hd_param["input_precision"]:
            area_accumulators = 0
        else:
            accumulator_output_pres = self.hd_param["weight_precision"] + self.hd_param["adc_resolution"] + self.hd_param["input_precision"] # output precision from adders_pv + required shifted bits
            nb_of_1b_adder_accumulator = accumulator_output_pres * self.dimensions["D1"] * self.dimensions["D3"]
            nb_of_1b_reg_accumulator = nb_of_1b_adder_accumulator # number of regs in an accumulator
            area_accumulators = self.logic_unit.get_1b_adder_area() * nb_of_1b_adder_accumulator + self.logic_unit.get_1b_reg_area() * nb_of_1b_reg_accumulator
        """total area of imc"""
        self.area = area_cells + area_mults + area_adcs + area_dacs + area_adders_pv + area_accumulators
        self.area_breakdown = { # unit: same with in input hd file
            "cells":    area_cells,
            "mults":    area_mults,
            "adcs":     area_adcs,
            "dacs":     area_dacs,
            "adders_pv":area_adders_pv,
            "accumulators": area_accumulators
        }
        return self.area

    def get_delay(self):
        """delay of imc (worst path: dacs -> mults -> adcs -> adders -> accumulators) """
        dly_dacs = self.get_dac_cost()[1]
        dly_mults = self.logic_unit.get_1b_multiplier_dly()
        dly_adcs = self.get_adc_cost()[1]
        """delay of adders_pv (adder type: RCA, worst path: in-to-sum -> in-to-sum -> ... -> in-to-cout -> cin-to-cout -> ... -> cin-to-cout)"""
        nb_inputs_of_adder_pv = self.hd_param["weight_precision"]
        if nb_inputs_of_adder_pv == 1:
            dly_adders_pv = 0
        else:
            adder_depth_pv = math.log2(nb_inputs_of_adder_pv)
            adder_depth_pv = int(adder_depth_pv)  # float -> int for simplicity
            adder_pv_output_precision = nb_inputs_of_adder_pv + self.hd_param["adc_resolution"] # output precision from adders_pv
            dly_adders_pv = (adder_depth_pv-1) * self.logic_unit.get_1b_adder_dly_in2sum() + self.logic_unit.get_1b_adder_dly_in2cout() + (adder_pv_output_precision-1) * self.logic_unit.get_1b_adder_dly_cin2cout()
        """delay of accumulators (adder type: RCA)"""
        if self.hd_param["input_bit_per_cycle"] == self.hd_param["input_precision"]:
            dly_accumulators = 0
        else:
            accumulator_input_pres = adder_pv_output_precision
            accumulator_output_pres = self.hd_param["weight_precision"] + self.hd_param["adc_resolution"] + self.hd_param["input_precision"]  # output precision from adders_pv + required shifted bits
            dly_accumulators = self.logic_unit.get_1b_adder_dly_in2cout() + (accumulator_output_pres-accumulator_input_pres) * self.logic_unit.get_1b_adder_dly_cin2cout()
        """total delay of imc"""
        self.delay = dly_dacs + dly_mults + dly_adcs + dly_adders_pv + dly_accumulators
        self.delay_breakdown = {
            "dacs":     dly_dacs,
            "mults":    dly_mults,
            "adcs":     dly_adcs,
            "adders_pv":dly_adders_pv,
            "accumulators": dly_accumulators
        }
        return self.delay

class ImcArray:
    def __init__(self, tech_param: Dict[str, float], hd_param: dict, dimensions: Dict[str, int]):
        """
        This class defines the general IMC array (including AIMC and DIMC)
        :param tech_param: definition of technology-related parameters
        :param hd_param: hardware architecture parameters except dimensions
        :param dimensions: dimensions definition
        """
        if hd_param["imc_type"] == "DIMC":
            self.unit = DimcArray(tech_param, hd_param, dimensions)
        elif hd_param["imc_type"] == "AIMC":
            self.unit = AimcArray(tech_param, hd_param, dimensions)
        self.total_area = self.unit.get_area()
        base_dims = [
            Dimension(idx, name, size)
            for idx, (name, size) in enumerate(dimensions.items())
        ]
        self.dimensions = base_dims
        self.dimension_sizes = [dim.size for dim in base_dims]
        self.nb_dimensions = len(base_dims)

if __name__ == "__main__":
#
##### IMC hardware dimension illustration (keypoint: adders' accumulation happens on D2)
#
#       |<------------------------ D1 ----------------------------->| (nb_of_columns/macro = D1 * weight_precision)
#    -  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    \
#    ^  +                                                           +  +  D3 (nb_of_macros)
#    |  +         ^     +++++++                                     +   +  \
#    |  +         |     +  W  +                                     +   +
#    |  +   group_depth +++++++                                     +   +
#    |  +         |     +  W  +                                     +   +
#    |  +         v     +++++++                                     +   +
#    |  +                  |                                        +   +
#    |  +                  v                                        +   +
#    |  +               multipliers -\                              +   +
#    |  +        .                    \                             +   +
#       +        .                     - adders (DIMC)              +   +
#   D2  +        .                    / OR adcs (AIMC)              +   +
#       +               multipliers -/       |                      +   +
#    |  +                  ^                 |                      +   +
#    |  +                  |                 |                      +   +
#    |  +         ^     +++++++              v                      +   +
#    |  +         |     +  W  +          adders_pv (place value)    +   +
#    |  +   group_depth +++++++              |                      +   +
#    |  +         |     +  W  +              v                      +   +
#    |  +         v     +++++++         accumulators                +   +
#    |  +                                    |                      +   +
#    v  +                                    |                      +   +
#    -  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   +
#          +                                 |                        + +
#           +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   (nb_of_rows/macro = D2 * group_depth)    |
#                                            v
#                                        outputs
#

    tech_param_28nm = {
        "tech_node":0.028,              # unit: um
        "vdd":      0.9,                # unit: V
        "nd2_cap":  0.7/1e3,            # unit: pF
        "xor2_cap": 0.7*1.5/1e3,        # unit: pF
        "dff_cap":  0.7*3/1e3,          # unit: pF
        "nd2_area": 0.614/1e6,          # unit: mm^2
        "xor2_area":0.614*2.4/1e6,      # unit: mm^2
        "dff_area": 0.614*6/1e6,        # unit: mm^2
        "nd2_dly":  0.0478,             # unit: ns
        "xor2_dly": 0.0478*2.4,         # unit: ns
        "dff_dly":  0.0478*3.4,         # unit: ns
    }
    dimensions = {
        "D1": 4,    # Do
        "D2": 32,   # Di (fix meaning: the dimension where addition happens)
        "D3": 1,    # nb_macros (fix meaning: nb_arrays)
    }  # {"D1": ("K", 4), "D2": ("C", 32),}

    """hd_param example for DIMC"""
    hd_param = {
        "imc_type":             "DIMC", # "DIMC" or "AIMC". Or else: pure digital
        "input_precision":      8,      # activation precison
        "weight_precision":     8,      # weight precision
        "input_bit_per_cycle":  1,      # nb_bits of input/cycle
        "group_depth":          1,      # m factor
        "w_energy_per_bit":     0.0101119640625,    # unit: pJ/bit
    }
    dimc = DimcArray(tech_param_28nm, hd_param, dimensions)
    #print(dimc.get_delay())
    #print(dimc.delay_breakdown)
    #exit()

    """hd_param example for AIMC"""
    hd_param_aimc = {
        "imc_type":             "AIMC", # "DIMC" or "AIMC". Or else: pure digital
        "input_precision":      8,      # activation precison
        "weight_precision":     8,      # weight precision
        "input_bit_per_cycle":  2,      # nb_bits of input/cycle
        "group_depth":          1,      # m factor
        "adc_resolution":       8,     # adc resolution
        "w_energy_per_bit":     0.0101119640625,    # unit: pJ/bit
    }
    hd_param_aimc["adc_resolution"] = hd_param_aimc["input_bit_per_cycle"] + 0.5*math.log2(dimensions["D2"])
    aimc = AimcArray(tech_param_28nm, hd_param_aimc, dimensions)
    print(aimc.get_delay())
    print(aimc.delay_breakdown)
    exit()
