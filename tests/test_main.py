import pytest

def test_descriptor():
    from rdkit.Chem import Descriptors
    # Was 209 but changed to 211 in Release_2023_09_1
    # Is 210 from Release_2023_09_3
    # Is 217 from Release_2024_09_4
    assert len(Descriptors._descList) == 217


def test_3d_descriptors():
    # from https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/UnitTestDescriptors.py
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors3D

    mol = Chem.MolFromSmiles('CCCO')
    
    # test function returns expected outputs
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    descs = Descriptors3D.CalcMolDescriptors3D(mol)
    assert 'InertialShapeFactor' in descs
    assert 20.9582649071385 == pytest.approx(descs['PMI1'], 1e-4)


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


def test_chemdraw_parser_roundtrip():
    # Exercises the new (2026.03.4) expat-based CDXML parser: write a molecule
    # to a ChemDraw block and read it back, checking the structure survives.
    from rdkit import Chem
    from rdkit.Chem import rdChemDraw

    mol = Chem.MolFromSmiles("NCc1ccccc1")
    block = rdChemDraw.MolToChemDrawBlock(mol)
    parsed = rdChemDraw.MolsFromChemDrawBlock(block)

    assert len(parsed) == 1
    assert Chem.MolToSmiles(Chem.RemoveHs(parsed[0])) == Chem.MolToSmiles(mol)
