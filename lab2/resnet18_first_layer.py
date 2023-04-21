workload = {
    0: {  # conv1, stride 2
        "operator_type": "Conv",
        "equation": "O[b][k][oy][ox]+=W[k][c][fy][fx]*I[b][c][iy][ix]",
        "dimension_relations": ["ix=2*ox+1*fx", "iy=2*oy+1*fy"],
        "loop_dim_size": {
            "B": 1,
            "K": 128,
            "C": 3,
            "OY": 112,
            "OX": 112,
            "FY": 7,
            "FX": 7,
        },
        "operand_precision": {"O": 16, "O_final": 8, "W": 8, "I": 8},
        "operand_source": {"W": [], "I": []},
        "constant_operands": ["I", "W"],
        "padding": {"IY": (3, 2), "IX": (3, 2)},
    },
}
