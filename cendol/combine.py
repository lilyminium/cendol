from collections import defaultdict
from typing import List, Dict, Set

from constructure.scaffolds import Scaffold
from constructure.constructors import Constructor

try:
    from .oefuncs import label_scaffold
except ImportError:
    from .rdfuncs import label_scaffold


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
        (e.g. 1 corresopnds to R1)
    
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


def cap_r_groups(
    constructor: Constructor,
    scaffold: Scaffold,
    substituents: Dict[int, List[str]] = {},
    r: int = 1,
    r_groups: List[int] = [],
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
    smiles = constructor.enumerate_combinations(scaffold, r_subs, validate=False)
    smiles = [x.replace(f"[{r}*]", "[R]") for x in smiles]
    smiles = [x.replace(f"[R{r}]", "[R]") for x in smiles]
    return smiles


def join_monomers(constructor: Constructor,
                  monomer_smiles: List[str] = [],
                  r_linkages: Dict[int, Set[int]] = {},
                  n_neighbor_monomers: int = 1,
                  label_central_atoms: bool = True,) -> List[str]:
    """Join monomers together at specified ``r_linkages``

    Parameters
    ----------
    constructor: Constructor
        Constructor object from constructure
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
    List of SMILES
    """
    if n_neighbor_monomers < 0:
        raise ValueError("Cannot have negative numbers of neighbors")

    # unify all links
    links = defaultdict(set)
    for i, partners in r_linkages.items():
        links[i] |= set(partners)
        for j in partners:
            links[j].add(i)

    # each monomer_smiles is a scaffold
    scaffolds = [Scaffold(smiles=m, r_groups={})  # skip r groups for now
                 for m in monomer_smiles]
    scaffold_r_groups = []
    for scaffold in scaffolds:
        r_groups = constructor.get_replaceable_r_groups(scaffold)
        scaffold_r_groups.append(r_groups)

    caps = {}

    # create capped substituents
    while n_neighbor_monomers:
        capped = defaultdict(list)
        for scaffold, r_groups in zip(scaffolds, scaffold_r_groups):
            for r in r_groups:
                smiles = cap_r_groups(constructor, scaffold, r=r,
                                      substituents=caps, r_groups=r_groups)
                capped[r].extend(smiles)

        caps = get_r_substituents(substituents=capped, r_linkages=links)
        n_neighbor_monomers -= 1

    # now combine with each central scaffold
    combinations = []
    for scaffold, r_groups in zip(scaffolds, scaffold_r_groups):
        smiles = cap_r_groups(constructor, scaffold, r=-1,
                              substituents=caps, r_groups=r_groups)
        if label_central_atoms:
            smiles = [label_scaffold(x, scaffold) for x in smiles]
        combinations.extend(smiles)

    return combinations
