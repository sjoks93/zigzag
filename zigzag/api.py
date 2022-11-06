from zigzag.classes.stages import *


def get_hardware_performance(onnx_model, accelerator, mapping=None):
    
    # Initialize the logger
    import logging as _logging
    _logging_level = _logging.INFO
    # _logging_format = '%(asctime)s - %(name)s.%(funcName)s +%(lineno)s - %(levelname)s - %(message)s'
    _logging_format = '%(asctime)s - %(funcName)s +%(lineno)s - %(levelname)s - %(message)s'
    _logging.basicConfig(level=_logging_level,
                        format=_logging_format)

    mainstage = MainStage([  # Initialize the MainStage as entry point
        ONNXModelParserStage,  # Parse the ONNX Model into the workload
        AcceleratorParserStage,  # Parse the accelerator module/passthrough given accelerator
        SimpleSaveStage,  # Save the summed CME to a json
        SumStage,  # Sum up the received best CME across all layers of he workload
        WorkloadStage,  # Iterate through the different layers in the workload
        MinimalLatencyStage,  # Reduce all CMEs, returning minimal latency one
        SpatialMappingGeneratorStage,  # Generate multiple spatial mappings (SM)
        MinimalLatencyStage,  # Reduce all CMEs, returning minimal latency one
        LomaStage,  # Generate multiple temporal mappings (TM)
        CostModelStage  # Evaluate generated SM and TM through cost model
    ],
        accelerator=accelerator,  # required by AcceleratorParserStage
        onnx_model=onnx_model,  # required by ONNXModelParserStage
        mapping_path=mapping,  # required by ONNXModelParserStage
        dump_filename_pattern="outputs/{datetime}.json",  # output file save pattern
        loma_lpf_limit=6  # required by LomaStage
    )

    # Launch the MainStage
    answers = mainstage.run()
    # Sanity check on the results
    assert len(answers) == 1, "Mainstage returned more than one CME."
    # Get CME from answer
    cme = answers[0][0]

    return cme.energy_total, cme.latency_total2