workload = {
    0: {  # conv1, stride 2
        "operator_type": "Conv",
        "equation": "O[b][k][oy][ox]+=W[k][c][fy][fx]*I[b][c][iy][ix]",
        "dimension_relations": ["ix=ox+1*fx", "iy=oy+1*fy"],
        "loop_dim_size": {
            "B": 1,
            "K": 32,
            "C": 32,
            "OY": 112,
            "OX": 112,
            "FY": 7,
            "FX": 7,
            #"C": "16"
        },
        "operand_precision": {"V": 16, "O": 8, "W": 8, "I": 8},
        #"operand_precision": {"O_final": 8, "O": 16, "W": 8, "I": 8},
        "operand_source": {"W": [], "I": []},
        "operand_state": {"V": []},
        "constant_operands": ["I", "W"],
        "state": True,
    },
}