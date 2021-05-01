import re

from openff.toolkit.topology import Molecule

try:
    from .oefuncs import sdf_to_offmols, fragment_into_substituent_smiles
except ImportError:
    from .rdfuncs import sdf_to_offmols, fragment_into_substituent_smiles


def get_breakable_bonds(offmol: Molecule, get_bonds_only: bool=True, n_neighbor_bonds: int=1):
    # single bonds -- break the middle one
    ATOM = "[!$(*#*)&!$(*=*)&A&!D1:{i}]"
    # CARBONYL = "[$(C=O)&A&!D1:{i}]"
    SINGLE = "-;!@".join([ATOM.format(i=i+1) for i in range(2*(n_neighbor_bonds+1))])
    # AMIDE_BONDS = [ATOM] * n_neighbor_bonds + [CARBONYL] + [ATOM] * (n_neighbor_bonds + 1)
    # AMIDE_BONDS = [x.format(i=i+1) for i, x in enumerate(AMIDE_BONDS)]
    # AMIDE = "-;!@".join(AMIDE_BONDS)
    try:
        matches = offmol.chemical_environment_matches(SINGLE)
    except Exception:  # stereochemistry error
        return set()
    # matches += offmol.chemical_environment_matches(AMIDE)
    unique_bonds = set()
    unique_matches = set()
    seen = set()
    
    for group in matches:
        i, j = n_neighbor_bonds, n_neighbor_bonds + 1
        if group[i] > group[j]:
            group = group[::-1]
        bond = (group[i], group[j])
        # avoid cutting every bond in a long chain

        if not any(x in seen for x in group):
            unique_bonds.add(bond)
            unique_matches.add(group)
            seen |= set(group)
    
    if get_bonds_only:
        return unique_bonds
    return unique_matches


def replace_dummy_with_R(smiles: str, number_r_groups: bool=True):
    PATTERN = r"\[\*:([0-9]+)\]"
    if number_r_groups:
        return re.sub(PATTERN, r"([R\1])", smiles)
    return re.sub(PATTERN, r"([R])", smiles)


def get_scaffolds(offmol: Molecule,
                  n_neighbor_bonds: int=1,
                  replace_with_r: bool=True):
    bonds = get_breakable_bonds(offmol, get_bonds_only=True,
                                n_neighbor_bonds=n_neighbor_bonds)
    smiles = fragment_into_substituent_smiles(offmol, bonds)
    if replace_with_r:
        smiles = [replace_dummy_with_R(s) for s in smiles]
    return smiles