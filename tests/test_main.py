def test_descriptor():
    from rdkit.Chem import Descriptors
    # Was 209 but changed to 211 in Release_2023_09_1
    # Is 210 from Release_2023_09_3
    assert len(Descriptors._descList) == 210


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
    import os

    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures

    fdefName = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    m = Chem.MolFromSmiles("OCc1ccccc1CN")
    feats = factory.GetFeaturesForMol(m)
    assert len(feats) == 8


def test_rdkit_chem_draw_import():
    # This segfaults if the compiled cairo version from centos is used
    from rdkit.Chem.Draw import ReactionToImage  # noqa: F401
