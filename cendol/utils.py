def replace_R_with_dummy(smiles: str):
    return re.sub(r"\(\[R([0-9]+)]\)", r"([\1*:\1])", smiles)

def replace_dummy_with_R(smiles: str):
    return re.sub(r"\[\*:([0-9]+)\]", r"([R\1])", smiles)