from typing import Set, Tuple

from openff.toolkit.topology import Molecule

try:
    from oefuncs import cleave_substituents_from_bond
except ImportError:
    from rdfuncs import cleave_substituents_from_bond

def get_heavy_degree(offmol: Molecule, atom_index: int) -> int:
    """Get heavy degree of atom in molecule"""
    atom = offmol.atoms[atom_index]
    return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)


def get_rotatable_bonds(offmol: Molecule) -> List[Tuple[int, int]]:
    ROT_BOND = "[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]"
    matches = offmol.chemical_environment_matches(ROT_BOND)
    unique_matches = {tuple(sorted(m)) for m in matches}
    non_terminal = {m for m in unique_matches
                    if all(get_heavy_degree(offmol, i) > 1 for i in m)}
    return sorted(non_terminal)



def get_dihedral_bonds(offmol: Molecule) -> List[Tuple[int, int]]:
    rotatable = get_rotatable_bonds(offmol)
    dihedral_bonds = []
    # remove if one end is symmetric
    for bond in rotatable:
        for i in bond:
            subs = cleave_substituents_from_bond(offmol, bond, i)
            if len(subs) == 1:
                break
        else:
            dihedral_bonds.append(bond)
    
    return dihedral_bonds

