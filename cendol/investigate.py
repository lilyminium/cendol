import re
from typing import List

from openff.toolkit.topology import Molecule

try:
    from .oefuncs import sdf_to_offmols, fragment_into_substituent_smiles
except ImportError:
    from .rdfuncs import sdf_to_offmols, fragment_into_substituent_smiles


def get_breakable_bonds(offmol: Molecule,
                        get_bonds_only: bool = True,
                        n_neighbor_bonds: int = 1) -> set:
    # single bonds -- break the middle one
    ATOM = "[!$(*#*)&!$(*=*)&A&!D1:{i}]"
    CHAIN = [ATOM.format(i=i + 1) for i in range(2 * (n_neighbor_bonds + 1))]
    SINGLE = "-;!@".join(CHAIN)
    try:
        matches = offmol.chemical_environment_matches(SINGLE)
    except Exception:  # stereochemistry error without custom type
        return set()
    
    unique_bonds = set()
    unique_matches = set()
    seen = set()

    for group in matches:
        # central atoms
        bond = (group[n_neighbor_bonds], group[n_neighbor_bonds + 1])
        if bond[0] > bond[1]:
            group = group[::-1]
        # avoid cutting every bond in a long chain
        if not any(x in seen for x in group):
            unique_bonds.add(bond)
            unique_matches.add(group)
            seen |= set(group)

    if get_bonds_only:
        return unique_bonds
    return unique_matches


def replace_dummy_with_R(smiles: str, number_r_groups: bool = True) -> str:
    PATTERN = r"\[\*:([0-9]+)\]"
    if number_r_groups:
        return re.sub(PATTERN, r"([R\1])", smiles)
    return re.sub(PATTERN, r"([R])", smiles)


def get_scaffolds(offmol: Molecule,
                  n_neighbor_bonds: int = 1,
                  replace_with_r: bool = True) -> List[str]:
    bonds = get_breakable_bonds(offmol, get_bonds_only=True,
                                n_neighbor_bonds=n_neighbor_bonds)
    smiles = fragment_into_substituent_smiles(offmol, bonds)
    if replace_with_r:
        smiles = [replace_dummy_with_R(s) for s in smiles]
    return smiles
