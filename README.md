# RDKit Python platform wheels for pip

Use [cibuildwheel](https://github.com/joerick/cibuildwheel) and Github Actions to build [RDKit](https://github.com/rdkit/rdkit) wheels for Linux and Mac OS. Wheels are available on PyPi using `pip`.

Available for

| OS | Version | Python |
| ----------- | ----------- | ----------- |
| Linux  | glibc >= 2.17 (e.g., Ubuntu 16.04+, CentOS 6+, ...) | 3.6, 3.7, 3.8, 3.9 |
| Mac OS | >= 10.9 (Mavericks)  | 3.6, 3.7, 3.8, 3.9 |

## Install RDKit 

### PIP

```bash
pip install rdkit-pypi
python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

### Poetry
```bash
poetry add rdkit-pypi
poetry run python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

## Testing: Build wheels locally (Linux only)

cibuildwheel uses `patchelf` (`apt install patchelf`) 

```bash
git clone https://github.com/kuelumbus/rdkit_platform_wheels.git
cd rdkit_platform_wheels

python3.8 -m pip install cibuildwheel

CIBW_BUILD_VERBOSITY=1 CIBW_MANYLINUX_X86_64_IMAGE=manylinux2014 CIBW_BEFORE_BUILD_LINUX="bash pre_linux.sh" cibuildwheel --platform linux --output-dir wheelhouse
```
