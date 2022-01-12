def test_descriptor():
    from rdkit.Chem import Descriptors

    assert len(Descriptors._descList) == 208


def test_3d_descriptors():
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors3D

    m2 = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(m2, randomSeed=1)
    assert round(Descriptors3D.NPR1(m2), 10) == 0.2553516286
