mapping = {
    "default": {  # first ResNet18 layer name in onnx model
        "spatial_mapping": {"D1": ("K", 32), "D2": ("C", 32)},
        "temporal_ordering": [
            # Innermost loop
            ("OX", 112),
            ("OY", 112),
            ("FX", 7),
            ("FY", 7),
            ("K", 2),
            # Outermost loop
        ],
        "core_allocation": 1,
        "memory_operand_links": {
            "O": "O",
            "W": "I2",
            "I": "I1",
        },
    },
}