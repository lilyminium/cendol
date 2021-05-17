
try:
    from .oefuncs import fragment_smiles, fragment_into_dummy_smiles
except ImportError:
    from .rdfuncs import fragment_smiles, fragment_into_dummy_smiles


def get_breakable_bonds(offmol, get_bonds_only=True, n_neighbor_bonds=1,
                        ignore_neighbors=False):
    # single bonds -- break the middle one
    ATOM = "[!$(*#*)&!$(*=*)&A&!D1:{i}]"
    SINGLE = "-;!@".join([ATOM.format(i=i+1) for i in range(2*(n_neighbor_bonds+1))])
    matches = offmol.chemical_environment_matches(SINGLE)
    unique_bonds = set()
    unique_matches = set()  # avoid long chains
    seen = set()
    
    for group in matches:
        i, j = n_neighbor_bonds, n_neighbor_bonds + 1
        if group[i] > group[j]:
            group = group[::-1]
        bond = (group[i], group[j])
        if bond not in unique_bonds:
            if not ignore_neighbors or not any(x in seen for x in group):
                unique_bonds.add(bond)
                unique_matches.add(group)
                seen |= set(group)
    
    if get_bonds_only:
        return unique_bonds
    return unique_matches


def fragment_into_substituent_smiles(offmol, cleave_bonds=[]):
    smiles = fragment_into_dummy_smiles(offmol, cleave_bonds)


def fragment_capped_monomers(fragmenter, smiles, combine=False):
    all_smiles = set()
    for smi in smiles:
        fragmented = fragment_smiles(fragmenter, smi, combine=combine)
        all_smiles |= set(fragmented)
    return all_smiles