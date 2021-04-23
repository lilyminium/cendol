import re

from rdkit import Chem
import numpy as np


def get_neighbor_atom_indices(mol, atom_indices=[]):
    neighbor_indices = []
    for i in atom_indices:
        atom = mol.GetAtomWithIdx(i)
        for nbr in atom.GetNeighbors():
            neighbor_indices.append(nbr.GetIdx())
    return neighbor_indices


def count_h_neighbors(mol, atom_indices=[]):
    neighbor_indices = get_neighbor_atom_indices(mol, atom_indices)
    return sum(1 for i in neighbor_indices if mol.GetAtomWithIdx(i).GetSymbol() == "H")


def label_scaffold(smiles, scaffold):
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
