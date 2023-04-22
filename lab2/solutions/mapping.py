mapping = {
    "/conv1/Conv": {  # first ResNet18 layer name in onnx model
        "spatial_mapping": {"D1": ("K", 32), "D2": ("C", 32)},
        "core_allocation": 1,
        "memory_operand_links": {
            "O": "O",
            "W": "I2",
            "I": "I1",
        },
    },
    "default": {  # used if layer.name is not found in dict
        "spatial_mapping": {"D1": ("K", 32), "D2": ("C", 32)},
        "core_allocation": 1,  # not important for now
        "memory_operand_links": {
            "O": "O",
            "W": "I2",
            "I": "I1",
        },  # not important for now
    },
}
