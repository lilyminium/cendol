from collections import defaultdict
from typing import List, Dict, Set
import itertools

from constructure.scaffolds import Scaffold

try:
    from .oefuncs import label_scaffold, Constructor, label_atoms
except ImportError:
    from .rdfuncs import label_scaffold, Constructor, label_atoms

from . import utils


def combine_substituents(constructor, scaffold, substituents,
                         monomer_smiles=[],
                         label_central_atoms=False,
                         label_monomers=False):
    r_subs = []
    for i, subs in substituents.items():
        r_subs.append([(i, sub) for sub in subs])
    combinations = map(dict, itertools.product(*r_subs))
    smiles = []
    for comb in combinations:
        smi = constructor.attach_substituents(scaffold, comb)
        if label_central_atoms or label_monomers:
            smi = label_atoms(smi, scaffold, comb,
                              label_central_atoms=label_central_atoms,
                              label_monomers=label_monomers,
                              monomer_smiles=monomer_smiles)
        smiles.append(smi)
    return smiles


def get_r_substituents(
    r_linkages: Dict[int, Set[int]] = {},
    substituents: Dict[int, List[str]] = {},
) -> Dict[int, List[str]]:
    """
    Create a dictionary of substituents to add to the scaffold,
    where each key corresponds to a labelled R group (e.g. 1 corresponds to R1)

    Parameters
    ----------
    r_linkages: dict (optional)
        Linkages between labelled R-groups. e.g. {1: [2, 3]} means that
        monomers will join together between R1-R2 and R1-R3.
    substituents: dict (optional)
        dictionary of capped substituents items, where
        each key corresponds to the labelled R-group that was not capped,
        (e.g. 1 corresponds to R1)
    
    Returns
    -------
    Dictionary of substituents to add to the scaffold
    """
    # translate links -> substituents
    link_subs = defaultdict(list)
    for i, partners in r_linkages.items():
        for j in partners:
            link_subs[i].extend(substituents[j])
    return link_subs


def cap_r_groups(constructor,
                 scaffold: Scaffold,
                 substituents: Dict[int, List[str]] = {},
                 r: int = 1,
                 r_groups: List[int] = [],
                 monomer_smiles=[],
                 label_central_atoms=False,
                 label_monomers=False,
                ) -> List[str]:
    """
    Cap each scaffold with substituents for each R-group except ``r``.

    Parameters
    ----------
    constructor: Constructor
        Constructor object from constructure
    scaffold: Scaffold
        Scaffold object from constructure with R groups defined
    substituents: dict (optional)
        A dictionary of substituents to add to the scaffold, where
        each key corresponds to a labelled R group (e.g. 1 corresponds to R1)
        and each value a list of SMILES patterns corresponding to possible
        substituents to attach. If groups are not available for a labelled
        R-group, it will be capped with a hydrogen.
    r: int (optional)
        R-group to not cap
    r_groups: list (optional)
        Available labelled R-groups

    Returns
    -------
    List of SMILES strings with [R] in place of uncapped R-group
    """
    if not r_groups:
        r_groups = constructor.get_replaceable_r_groups(scaffold)
    cap_r_groups = [i for i in r_groups if i != r]
    r_subs = {i: substituents.get(i, ["[R][H]"]) for i in cap_r_groups}
    smiles = combine_substituents(constructor, scaffold, r_subs,
                                  monomer_smiles=monomer_smiles,
                                  label_central_atoms=label_central_atoms,
                                  label_monomers=label_monomers)
    smiles = [x.replace(f"[{r}*]", "[R]") for x in smiles]
    smiles = [x.replace(f"[R{r}]", "[R]") for x in smiles]
    return smiles


def join_monomers(monomer_smiles: List[str] = [],
                  r_linkages: Dict[int, Set[int]] = {},
                  n_neighbor_monomers: int = 1,
                  label_central_atoms: bool = True,
                  label_monomers: bool = False) -> List[str]:
    """Join monomers together at specified ``r_linkages``

    Parameters
    ----------
    monomer_smiles: list (optional)
        List of SMILES with labelled R-groups describing
        scaffolds for constructure. e.g. "C([R1])([R2])([R3])"
    r_linkages: dict (optional)
        Linkages between labelled R-groups. e.g. {1: [2, 3]} means that
        monomers will join together between R1-R2 and R1-R3.
    n_neighbor_monomers: int (optional)
        How many neighbors to add onto each labelled R-group
    label_central_atoms: bool (optional)
        Whether to label the atoms of the original scaffold with
        atom map numbers

    Returns
    -------
    List of list of SMILES
    """
    constructor = Constructor()

    if n_neighbor_monomers < 0:
        raise ValueError("Cannot have negative numbers of neighbors")
    
    if not r_linkages:
        r_groups = set()
        for smi in monomer_smiles:
            rs = utils.get_r_group_numbers_from_smiles(smi)
            r_groups |= set(rs)
        r_linkages = {i: [i] for i in r_groups}

    # unify all links (i.e. A=>B means A<=>B)
    links = defaultdict(set)
    for i, partners in r_linkages.items():
        links[i] |= set(partners)
        for j in partners:
            links[j].add(i)

    # each monomer_smiles is a scaffold
    scaffolds = []
    scaffold_r_numbers = []
    for i, smi in enumerate(monomer_smiles, 1):
        # skip r groups for now
        scaffold = Scaffold(smiles=smi, r_groups={})
        scaffolds.append(scaffold)
        r_groups = constructor.get_replaceable_r_groups(scaffold)
        scaffold_r_numbers.append(r_groups)

    caps = {}

    # create capped substituents
    while n_neighbor_monomers:
        capped = defaultdict(list)
        for scaffold, r_groups in zip(scaffolds, scaffold_r_numbers):
            for r in r_groups:
                smiles = cap_r_groups(constructor, scaffold, r=r,
                                      substituents=caps, r_groups=r_groups)
                capped[r].extend(smiles)

        caps = get_r_substituents(substituents=capped, r_linkages=links)
        n_neighbor_monomers -= 1

    # now combine with each central scaffold
    combinations = []
    for scaffold, r_groups in zip(scaffolds, scaffold_r_numbers):
        smiles = cap_r_groups(constructor, scaffold, r=-1,
                              substituents=caps, r_groups=r_groups,
                              monomer_smiles=monomer_smiles,
                              label_central_atoms=label_central_atoms,
                              label_monomers=label_monomers)
        combinations.append(smiles)

    return combinations
