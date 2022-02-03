def test_descriptor():
    from rdkit.Chem import Descriptors

    assert len(Descriptors._descList) == 208


def test_3d_descriptors():
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors3D

    m2 = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(m2, randomSeed=1)
    assert round(Descriptors3D.NPR1(m2), 10) == 0.2553516286


def test_data_dir_and_chemical_features():
    """Checks if data directory is correctly set
    and if ChemicalFeatures work
    """
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig
    from rdkit import Chem
    import os

    fdefName = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    m = Chem.MolFromSmiles("OCc1ccccc1CN")
    feats = factory.GetFeaturesForMol(m)
    assert len(feats) == 8
