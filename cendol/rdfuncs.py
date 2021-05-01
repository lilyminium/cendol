import re
from typing import List, Tuple, Set

from constructure.scaffolds import Scaffold
from rdkit import Chem
import numpy as np

from openff.toolkit.topology import Molecule

def fragment_molecule(fragmenter, rdmol, combine: bool=False, **kwargs):
    molecule = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    central = set()
    for atom in rdmol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            central.add(atom.GetIdx())
    fragger = fragmenter()
    results = fragger.fragment(molecule)
    fragments = []
    for (a1, a2), oemol in results.fragments_by_bond.items():
        if a1 in central or a2 in central:
            fragments.append(oemol)
    if not combine:
        if not fragments:
            return [rdmol]
        mols = [x.to_rdkit() for x in fragments]
        return mols
    # match all fragments to original molecule
    fragment_indices = set()
    for oemol in fragments:
        rdfrag = Chem.MolFromSmiles(oemol.smiles)
        rdfrag = Chem.AddHs(rdfrag)
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
    if not rwmol.GetNumAtoms():
        return [rdmol]
    Chem.SanitizeMol(rwmol)
    return [rwmol]


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
    fragmented = fragment_molecule(Chem.MolFromSmiles(smiles), combine=combine, **kwargs)
    return [Chem.MolToSmiles(x, isomericSmiles=True, allHsExplicit=True) for x in fragmented]


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


def rdmol_from_pdb_and_smiles(pdbfile, smiles):
    pdb = Chem.MolFromPDBFile(pdbfile, removeHs=False)
    return map_smiles_onto_pdb_rdmol(pdb, smiles)
    
    
def map_smiles_onto_pdb_rdmol(rdmol, smiles):
    pdb = rdmol
    smi = Chem.MolFromSmiles(smiles)
    pdb_no_h_to_h_index = {}
    for atom in pdb.GetAtoms():
        if atom.GetSymbol() != "H":
            pdb_no_h_to_h_index[len(pdb_no_h_to_h_index)] = atom.GetIdx()
            
    smi_no_h_to_h_index = {}
    for atom in smi.GetAtoms():
        if atom.GetSymbol() != "H":
            smi_no_h_to_h_index[len(smi_no_h_to_h_index)] = atom.GetIdx()
            
    pdb2 = Chem.RemoveHs(pdb)
    smi2 = Chem.RemoveHs(smi)
    # PDB reads in all bonds as single, I think, and all atoms as no charge
    for bond in smi2.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in smi2.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(False)
    Chem.SanitizeMol(smi2)
    smi2 = Chem.AddHs(smi2)
    smi2 = Chem.RemoveHs(smi2)
    
    match = pdb2.GetSubstructMatch(smi2)
    if not len(match):
        raise ValueError("PDB file and SMILES do not match")
    
    for smi2_index, pdb2_index in enumerate(match):
        smi_atom = smi.GetAtomWithIdx(smi_no_h_to_h_index[smi2_index])
        pdb_atom = pdb.GetAtomWithIdx(pdb_no_h_to_h_index[pdb2_index])
        pdb_atom.SetFormalCharge(smi_atom.GetFormalCharge())
    
    for bond in smi2.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        smi_i = smi_no_h_to_h_index[i]
        smi_j = smi_no_h_to_h_index[j]
        pdb_i = pdb_no_h_to_h_index[match[i]]
        pdb_j = pdb_no_h_to_h_index[match[j]]
        smi_bond = smi.GetBondBetweenAtoms(smi_i, smi_j)
        pdb_bond = pdb.GetBondBetweenAtoms(pdb_i, pdb_j)
        pdb_bond.SetBondType(smi_bond.GetBondType())

    Chem.SanitizeMol(pdb)
    return pdb
    
