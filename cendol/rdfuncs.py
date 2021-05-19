import re
from typing import List, Tuple, Set

from constructure.scaffolds import Scaffold
from constructure.constructors import RDKitConstructor as Constructor
from rdkit import Chem
import numpy as np

from openff.toolkit.topology import Molecule

from . import utils

def label_atoms(smiles, scaffold, substituents,
                label_central_atoms=False,
                label_monomers=False,
                monomer_smiles=[]):
    dummy_substituents = {}
    for i, sub in substituents.items():
        dummy_substituents[i] = clip_r_from_smiles(sub, remove_hs=False, label=False)
    rdmol = Chem.MolFromSmiles(smiles)
    rdmol = Chem.AddHs(rdmol)
    rdmol = label_scaffold_atoms(rdmol, scaffold.smiles, dummy_substituents)
    if label_monomers:
        rdmol = label_monomer_atoms(rdmol, monomer_smiles=monomer_smiles)
    # if not label_central_atoms:
    #     for atom in rdmol.GetAtoms():
    #         if atom.GetAtomMapNum() > 0:
    #             atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(rdmol, isomericSmiles=True, allHsExplicit=True)


def label_monomer_atoms(rdmol, monomer_smiles=[]):
    clipped = [clip_r_from_smiles(smi, label=False, remove_hs=False) for smi in monomer_smiles]
    clipped_mols = [Chem.MolFromSmarts(smi) for smi in clipped]

    for n_mol, monomer in enumerate(clipped_mols, 1):
        matches = rdmol.GetSubstructMatches(monomer)
        for match in matches:
            anums = [rdmol.GetAtomWithIdx(i).GetAtomMapNum() for i in match]
            if not any(i != 0 for i in anums):
                for atom_index in match:
                    atom = rdmol.GetAtomWithIdx(atom_index)
                    atom.SetAtomMapNum(-n_mol)
    return rdmol

    

def is_correct_match(rdmol, match_indices, index_to_r, sub_mols):
    rdcopy = Chem.RWMol(rdmol)
    for atom in rdcopy.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    for index in sorted(match_indices)[::-1]:
        if index not in index_to_r:
            rdcopy.RemoveAtom(index)
    frag_indices = []
    fragments = Chem.GetMolFrags(rdcopy, asMols=True, sanitizeFrags=False,
                                 fragsMolAtomMapping=frag_indices)
    for frag, indices in zip(fragments, frag_indices):
        for index_ in indices:
            index = rdcopy.GetAtomWithIdx(index_).GetAtomMapNum()
            if index in index_to_r:
                r_num = index_to_r[index]
                break
        frag = Chem.RemoveHs(frag)
        Chem.SanitizeMol(frag)
        sub = sub_mols[r_num]
        # print(r_num)
        sub = Chem.RemoveAllHs(sub_mols[r_num])
        match = frag.GetSubstructMatch(sub)
        # match = frag.GetSubstructMatch(sub_mols[r_num])
        if not match:
            # print("not match")
            # print(Chem.MolToSmiles(frag, isomericSmiles=True))
            # print(Chem.MolToSmiles(sub, isomericSmiles=True))
            return False
    return True

def label_scaffold_atoms(rdmol, r_group_smiles, dummy_substituents):
    dummy_smiles = utils.replace_R_with_dummy(r_group_smiles)
    wild_smiles = utils.replace_dummy_with_wildcard(dummy_smiles)
    pattern = Chem.MolFromSmarts(wild_smiles)
    matches = rdmol.GetSubstructMatches(pattern)
    if not matches:
        raise ValueError(f"Could not label scaffold with r_group {r_group_smiles}, "
                         f"{Chem.MolToSmiles(rdmol, isomericSmiles=True, allHsExplicit=True)}, "
                         f"{dummy_substituents}")

    # which match is it? Test by cutting bonds and comparing
    # to substituents
    r_pattern_indices = {}
    for atom in pattern.GetAtoms():
        anum = atom.GetAtomMapNum()
        if anum != 0:
            r_pattern_indices[anum] = atom.GetIdx()
    
    sub_mols = {}
    for i, smi in dummy_substituents.items():
        mol = Chem.MolFromSmarts(smi)
        sub_mols[i] = mol
    
    for match in matches:
        index_to_r = {match[i]: a for a, i in r_pattern_indices.items()}
        indices = [x for x in match if x not in index_to_r]
        if is_correct_match(rdmol, match, index_to_r, sub_mols):
            break
    for i, index in enumerate(indices, 1):
        rdmol.GetAtomWithIdx(index).SetAtomMapNum(i)
        
    return rdmol
    
def set_atom_numbers(smiles, number=0):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(int(number))
    return Chem.MolToSmiles(mol, isomericSmiles=True, allHsExplicit=True)

def fragment_into_dummy_smiles(offmol, cleave_bonds=[]):
    rdmol = Chem.RWMol(offmol.to_rdkit())
    dummy = Chem.Atom("*")
    r_linkages = {}
    counter = 1
    for bond in cleave_bonds:
        bond_type = rdmol.GetBondBetweenAtoms(*bond).GetBondType()
        rdmol.RemoveBond(*bond)
        r_linkages[counter] = [counter + 1]
        for atom_index in bond:
            dummy_copy = Chem.Atom(dummy)
            dummy_copy.SetAtomMapNum(counter)
            new_atom_index = rdmol.AddAtom(dummy_copy)
            rdmol.AddBond(atom_index, new_atom_index, bond_type)
            counter += 1
    mols = Chem.GetMolFrags(rdmol, asMols=True)
    smiles = [Chem.MolToSmiles(m, isomericSmiles=True, allHsExplicit=True)
              for m in mols]
    return smiles, r_linkages


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
    

def clip_r_from_smiles(smiles: str, remove_hs: bool = True, label: bool = True) -> str:
    subbed = utils.replace_R_with_dummy(smiles)
    rdmol = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(subbed)))
    to_remove = []
    # loop over this twice just in case
    # removing atoms screws up the iteration
    for atom in rdmol.GetAtoms():
        if atom.HasProp("dummyLabel"):
            to_remove.append(atom.GetIdx())
    for idx in sorted(to_remove)[::-1]:
        rdmol.RemoveAtom(idx)
    rdmol.UpdatePropertyCache()
    if remove_hs:
        rdmol = Chem.RemoveHs(rdmol, sanitize=False)
    if label:
        for atom in rdmol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
    return Chem.MolToSmiles(rdmol, isomericSmiles=True, allHsExplicit=True)


def get_subrdmol(rdmol, indices: list=[],
                 label_indices: list=[],
                 sanitize: bool=False):
    """Create new sub-molecule from selected atom indices
    
    Parameters
    ----------
    rdmol: rdkit.Chem.rdchem.Mol
        Input molecule
    indices: iterable of ints
        atom indices to include from input molecule, indexed from 0
    sanitize: bool
        whether to sanitize the molecule (recommend: no)
        
    Returns
    -------
    rdkit.Chem.rdchem.Mol: subset of molecule
    """
    submol = Chem.RWMol(rdmol)
    for i, lix in enumerate(label_indices, 1):
        at = submol.GetAtomWithIdx(int(lix))
        at.SetAtomMapNum(i)
    indices = sorted(set(indices) | set(label_indices))

    ix = sorted([at.GetIdx() for at in rdmol.GetAtoms()
                 if at.GetIdx() not in indices])
    for i in ix[::-1]:
        submol.RemoveAtom(int(i))
    if sanitize:
        Chem.SanitizeMol(submol)
    return submol

def get_sub_smarts(mol, indices: list=[], label_indices: list=[]):
    """Get SMARTS of selected atoms in molecule
    
    Parameters
    ----------
    mol: openff.toolkit.topology.molecule.Molecule
        Input molecule
    indices: iterable of ints
        atom indices to include in SMARTS, indexed from 0
    label_indices: iterable of ints
        atom indices to label, indexed from 0. The atoms
        will be labelled in the order specified. Labels
        begin from 1.
        
    Returns
    -------
    str: SMARTS string
    """
    rdmol = mol.to_rdkit()
    submol =  get_subrdmol(rdmol, indices=indices,
                           label_indices=label_indices)
    return Chem.MolToSmarts(submol, isomericSmiles=False)

def get_sub_smiles(mol, indices: list=[], label_indices: list=[]):
    """Get SMARTS of selected atoms in molecule
    
    Parameters
    ----------
    mol: openff.toolkit.topology.molecule.Molecule
        Input molecule
    indices: iterable of ints
        atom indices to include in SMARTS, indexed from 0
    label_indices: iterable of ints
        atom indices to label, indexed from 0. The atoms
        will be labelled in the order specified. Labels
        begin from 1.
        
    Returns
    -------
    str: SMARTS string
    """
    rdmol = mol.to_rdkit()
    submol =  get_subrdmol(rdmol, indices=indices,
                           label_indices=label_indices)
    return Chem.MolToSmiles(submol, allHsExplicit=True)