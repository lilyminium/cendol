
try:
    from .oefuncs import fragment_smiles
except ImportError:
    from .rdfuncs import fragment_smiles


def fragment_capped_monomers(fragmenter, smiles, combine=False):
    all_smiles = set()
    for smi in smiles:
        fragmented = fragment_smiles(fragmenter, smi, combine=combine)
        all_smiles |= set(fragmented)
    return all_smiles