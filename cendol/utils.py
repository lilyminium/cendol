from typing import List
import re

def replace_R_with_dummy(smiles: str):
    smiles = re.sub(r"\(\[R([0-9]+)]\)", r"([\1*:\1])", smiles)
    smiles = re.sub(r"\[R\]", r"[*]", smiles)
    return smiles

def replace_dummy_with_R(smiles: str):
    return re.sub(r"\[\*:([0-9]+)\]", r"([R\1])", smiles)

def get_r_group_numbers_from_smiles(smiles: str) -> List[int]:
    return list(map(int, re.findall(r"\(\[R([0-9]+)]\)", smiles)))

def replace_dummy_with_wildcard(smiles: str):
    return re.sub(r"\[[0-9]*\*:([0-9]+)\]", r"[*:\1]", smiles)