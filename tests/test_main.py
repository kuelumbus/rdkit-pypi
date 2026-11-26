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
    from rdkit.Chem import Descriptors3D

    # Fixed coordinates
    mol = Chem.MolFromSmiles('CCCO |(1.44534,-0.585581,0.158885;0.667797,0.646552,-0.278384;'
                             '-0.741018,0.544094,0.296045;-1.37212,-0.605065,-0.176546)|')
    descs = Descriptors3D.CalcMolDescriptors3D(mol)
    assert 'InertialShapeFactor' in descs
    assert descs['PMI1'] == pytest.approx(20.9583, abs=1e-4)


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
    """
    Exercises the new (2026.03.4) expat-based CDXML parser: write a molecule
    to a ChemDraw block and read it back, checking the structure survives.
    """
    from rdkit import Chem
    from rdkit.Chem import rdChemDraw

    mol = Chem.MolFromSmiles("NCc1ccccc1")
    block = rdChemDraw.MolToChemDrawBlock(mol)
    parsed = rdChemDraw.MolsFromChemDrawBlock(block)

    assert len(parsed) == 1
    assert Chem.MolToSmiles(Chem.RemoveHs(parsed[0])) == Chem.MolToSmiles(mol)


def test_all_compiled_extensions_import():
    """
    Import every compiled extension shipped in the wheel
    """
    import importlib
    import pathlib

    import rdkit

    root = pathlib.Path(rdkit.__file__).parent
    modules = []
    for path in sorted(root.rglob("*")):
        if path.suffix not in (".so", ".pyd"):
            continue
        rel = path.relative_to(root)
        # the bundled shared libraries are not importable modules
        if rel.parts[0] in ("rdkit.libs", ".dylibs"):
            continue
        modules.append("rdkit." + ".".join(rel.parts[:-1] + (rel.name.split(".")[0],)))

    failed = []
    for name in modules:
        try:
            importlib.import_module(name)
        except BaseException as exc:  # noqa: BLE001 - report, don't mask
            failed.append(f"{name}: {type(exc).__name__}: {exc}")

    assert not failed, "extension modules failed to import:\n" + "\n".join(failed)


def test_inchi():
    # RDK_BUILD_INCHI_SUPPORT
    from rdkit import Chem
    from rdkit.Chem import inchi

    assert inchi.INCHI_AVAILABLE
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    assert inchi.MolToInchi(mol).startswith("InChI=1S/C9H8O4")
    assert len(inchi.MolToInchiKey(mol).split("-")) == 3


def test_avalon():
    # RDK_BUILD_AVALON_SUPPORT
    from rdkit import Chem
    from rdkit.Avalon import pyAvalonTools

    assert pyAvalonTools.GetCanonSmiles("c1ccccc1O", True) == "Oc1ccccc1"
    fp = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    assert fp.GetNumOnBits() > 0


def test_freesasa():
    # RDK_BUILD_FREESASA_SUPPORT
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFreeSASA

    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    radii = rdFreeSASA.classifyAtoms(mol)
    assert len(radii) == mol.GetNumAtoms()
    assert rdFreeSASA.CalcSASA(mol, radii) > 0


def test_yaehmop():
    # RDK_BUILD_YAEHMOP_SUPPORT
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdEHTTools

    mol = Chem.AddHs(Chem.MolFromSmiles("CO"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    ok, res = rdEHTTools.RunMol(mol)
    assert ok
    assert res.totalEnergy < 0
    assert len(res.GetAtomicCharges()) == mol.GetNumAtoms()


def test_xyz2mol():
    # RDK_BUILD_XYZ2MOL_SUPPORT
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    xyz = (
        "3\n\n"
        "O      0.000000    0.000000    0.000000\n"
        "H      0.758602    0.000000    0.504284\n"
        "H     -0.758602    0.000000    0.504284\n"
    )
    mol = Chem.MolFromXYZBlock(xyz)
    rdDetermineBonds.DetermineBonds(mol, charge=0)
    assert Chem.MolToSmiles(mol) == "[H]O[H]"


def test_cairo_png_rendering():
    # RDK_BUILD_CAIRO_SUPPORT
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, Chem.MolFromSmiles("c1ccccc1O"))
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    assert png.startswith(b"\x89PNG")
    assert len(png) > 1000


def test_pickle_roundtrip():
    import pickle

    from rdkit import Chem

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    assert Chem.MolToSmiles(pickle.loads(pickle.dumps(mol))) == Chem.MolToSmiles(mol)


def test_numpy_interop():
    # The wheels are built against NumPy 2.x
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)
    arr = np.array(fp)
    assert arr.shape == (2048,)
    assert arr.sum() == fp.GetNumOnBits()


def test_molblock_roundtrip():
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    parsed = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
    assert Chem.MolToSmiles(parsed) == Chem.MolToSmiles(mol)
