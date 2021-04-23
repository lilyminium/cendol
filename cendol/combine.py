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
    cap_r_groups = [i for i in r_groups if i != r]
    r_subs = {i: substituents.get(i, ["[R][H]"]) for i in cap_r_groups}
    smiles = constructor.enumerate_combinations(scaffold, r_subs, validate=False)
    smiles = [x.replace(f"[{r}*]", "[R]") for x in smiles]
    smiles = [x.replace(f"[R{r}]", "[R]") for x in smiles]
    return smiles


def join_monomers(
    constructor: Constructor,
    monomer_smiles: List[str] = [],
    r_linkages: Dict[int, Set[int]] = {},
    n_neighbor_monomers: int = 1,
    label_central_atoms: bool = True,
):
    if n_neighbor_monomers < 0:
        raise ValueError("Cannot have negative numbers of neighbors")

    # unify all links
    links = defaultdict(set)
    for i, partners in r_linkages.items():
        links[i] |= set(partners)
        for j in partners:
            links[j].add(i)

    # each monomer_smiles is a scaffold
    scaffolds = [
        Scaffold(smiles=m, r_groups={}) for m in monomer_smiles  # skip r groups for now
    ]
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
                smiles = cap_r_groups(
                    constructor, scaffold, substituents=caps, r=r, r_groups=r_groups
                )
                capped[r].extend(smiles)

        caps = get_r_substituents(substituents=capped, r_linkages=links)
        n_neighbor_monomers -= 1

    # now combine with each central scaffold
    combinations = []
    for scaffold, r_groups in zip(scaffolds, scaffold_r_groups):
        smiles = cap_r_groups(
            constructor, scaffold, substituents=caps, r=-1, r_groups=r_groups
        )
        if label_central_atoms:
            smiles = [label_scaffold(x, scaffold) for x in smiles]
        combinations.extend(smiles)

    return combinations
