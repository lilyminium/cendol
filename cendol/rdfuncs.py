import re
from typing import List, Tuple, Set

from constructure.scaffolds import Scaffold
from rdkit import Chem
import numpy as np

from openff.toolkit.topology import Molecule

def fragment_smiles(fragmenter, smiles: str, combine: bool=False, **kwargs) -> List[str]:
    """
    Fragment SMILES string.

    Parameters
    ----------
    fragmenter:
        fragmenter class
    smiles: str
        SMILES string
    combine: bool (optional)
        Whether to combine resulting fragments so that original
        residue is whole. This requires the original residue to
        be labelled with atom map numbers.
    **kwargs:
        passed to fragmenter.fragment()

    Returns
    -------
    List of isomeric SMILES fragments
    """
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=False)
    rdmol = Chem.MolFromSmiles(smiles)
    central = set()
    for atom in rdmol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            central.add(atom.GetIdx())
    fragger = fragmenter(mol)
    fragger.fragment(**kwargs)
    fragments = []
    for (a1, a2), oemol in fragger.fragments.items():
        if a1 in central or a2 in central:
            fragments.append(oemol)
    if not combine:
        mols = [x.to_rdkit() for x in fragments]
        return [Chem.MolToSmiles(x, isomericSmiles=True, allHsExplicit=True) for x in mols]
    # match all fragments to original molecule
    fragment_indices = set()
    for oemol in fragments:
        rdfrag = Chem.MolFromSmiles(oemol.to_smiles())
        matches = rdmol.GetSubstructMatches(rdfrag)
        for m in matches:
            mset = set(m)
            if mset.intersection(central):
                fragment_indices |= mset
    atoms_to_keep = [i for i in range(rdmol.GetNumAtoms()) if i not in fragment_indices]
    rwmol = Chem.RWMol(rdmol)
    for idx in atoms_to_keep[::-1]:
        rwmol.RemoveAtom(idx)
    rwmol.UpdatePropertyCache()
    Chem.SanitizeMol(rwmol)
    fragments = [Chem.MolToSmiles(rwmol, isomericSmiles=True, allHsExplicit=True)]
    return fragments


def get_neighbor_atom_indices(mol, atom_indices: List[int]=[]) -> List[int]:
    """
    Get indices of neighboring atoms

    Parameters
    ----------
    mol: RDKit.Chem.rdchem.Mol
        RDKit molecule
    atom_indices: list (optional)
        Atom indices to get the neighbors for

    Returns
    -------
    Indices of neighboring atoms
    """
    neighbor_indices = []
    for i in atom_indices:
        atom = mol.GetAtomWithIdx(i)
        for nbr in atom.GetNeighbors():
            neighbor_indices.append(nbr.GetIdx())
    return neighbor_indices


def count_h_neighbors(mol, atom_indices: List[int]=[]) -> int:
    """
    Count number of neighboring atoms that are hydrogen

    Parameters
    ----------
    mol: RDKit.Chem.rdchem.Mol
        RDKit molecule
    atom_indices: list (optional)
        Atom indices to get the neighbors for

    Returns
    -------
    Number of neighboring hydrogens
    """
    neighbor_indices = get_neighbor_atom_indices(mol, atom_indices)
    return sum(1 for i in neighbor_indices if mol.GetAtomWithIdx(i).GetSymbol() == "H")


def label_scaffold(smiles: str, scaffold: Scaffold) -> str:
    """
    Label region of central scaffold pattern in the ``smiles`` pattern
    with atom map numbers

    Parameters
    ----------
    smiles: str
        SMILES string
    scaffold: Scaffold
        Scaffold with pattern to label in ``smiles``

    Returns
    -------
    Labelled SMILES string
    """
    sub = re.sub(r"\(([\\/]*)\[R([0-9]+)]\)", "", scaffold.smiles)
    pattern = Chem.MolFromSmiles(sub)
    mol = Chem.MolFromSmiles(smiles)
    newmol = Chem.AddHs(mol)
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return smiles
    n_hs = [count_h_neighbors(newmol, indices) for indices in matches]
    central_i = np.argmin(n_hs)
    #     n_orig = mol.GetNumAtoms()
    for i, index in enumerate(matches[central_i], 1):
        try:
            mol.GetAtomWithIdx(index).SetAtomMapNum(i)
        except RuntimeError:
            pass
    return Chem.MolToSmiles(mol)


def cleave_substituents_from_bond(offmol: Molecule,
                                  bond_tuple: Tuple[int, int],
                                  atom_index: int) -> List[str]:
    rdmol = Chem.RWMol(offmol.to_rdkit())
    rdatom = rdmol.GetAtomWithIdx(atom_index)
    sub_atoms = set()
    for bond in rdatom.GetBonds():
        other_idx = bond.GetOtherAtomIdx(atom_index)
        bond_atoms = tuple(sorted([atom_index, other_idx]))
        if bond_atoms != bond_tuple:
            sub_atoms.add(other_idx)
    
    for other_index in sub_atoms:
        rdmol.RemoveBond(atom_index, other_index)
    rdmol.UpdatePropertyCache()

    atoms_to_keep = set()
    fragments = Chem.GetMolFrags(rdmol, sanitizeFrags=False)
    for frag in fragments:
        if any(i in frag for i in sub_atoms):
            atoms_to_keep.add(tuple(sorted(frag)))
    
    smiles = []
    for atoms in atoms_to_keep:
        root_atom = [i for i in sub_atoms if i in atoms][0]
        smiles.append(Chem.MolFragmentToSmiles(rdmol, list(atoms)),
                      rootedAtAtom=root_atom, canonical=True,
                      isomericSmiles=True, allHsExplicit=True)
    return smiles