import importlib


from cendol import join_monomers
from constructure.constructors import Constructor, OpenEyeConstructor, RDKitConstructor

try:
    importlib.import_module("openeye")
    CONSTRUCTORS = [OpenEyeConstructor, RDKitConstructor]
except ImportError:
    CONSTRUCTORS = [RDKitConstructor]


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
@pytest.mark.parametrize("n_neighbor_monomers, n_combinations, last", [
    (1, 110, 'C=C(O)C/[CH:1]=[C:2](\\[CH2:3]/C=C(\\C)O)[OH:4]'),
    (2, 1802, 'C=C(O)C/C=C(/O)C/C=C(/O)C/[CH:1]=[C:2](\\[CH2:3]/C=C(\\C)O)[OH:4]'),
])
def test_join_multiple_rs(constructor, n_neighbor_monomers, n_combinations,
                          last):
    monomer1_smiles = [
        "C(/[R2])",  # start cap
        "C(=C)(/[R1])",  # end cap
        "C(=C(/C([R2]))N)(/[R1])",  # res 1
        "C(=C(/C)C([R2]))(/[R1])",  # res 2
        "C(=C(/C([R2]))CCC)(/[R1])",  # res 3
        "C(=C(/C([R2]))O)(/[R1])",  # res 4
    ]

    link1_smiles = {1: [2]}

    smiles = join_monomers(constructor, monomer1_smiles, link1_smiles,
                           n_neighbor_monomers=n_neighbor_monomers)
    assert len(smiles) == n_combinations
    assert smiles[-1] == last


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
def test_unique_rs(constructor):
    monomer2_smiles = [
        "C([R1])",  # start cap
        "C(=C)([R10])",  # end cap
        "C(=C(C([R3]))N)([R2])",  # res 1
        "C(=C(C)C([R5]))([R4])",  # res 2
        "C(=C(C([R7]))CCC)([R6])",  # res 3
        "C(=C(C([R9]))O)([R8])",  # res 4
    ]

    link2_smiles = {1: [2], 3: [4], 5: [6], 7: [8], 9: [10]}
    smiles2 = join_monomers(monomer2_smiles, link2_smiles)
    assert len(smiles2) == 8
    assert smiles2[-1] == 'C=C[CH2:3][C:2](=[CH:1]CC(=C)CCC)[OH:4]'


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
def test_branching_monomers(constructor):
    monomer3_smiles = [
        "C([R1])",  # cap
        "C(COCN([R3])([R4]))([R2])", # res
    ]
    link3_smiles = {2: [1], 3: [1, 2], 4: [1, 2],}

    smiles3 = join_monomers(monomer3_smiles, link3_smiles)
    assert len(smiles3) == 8
    assert smiles3[-1] == 'CCOCN[CH2:1][CH2:2][O:3][CH2:4][N:5](CCOCN)CCOCN'