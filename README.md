# 🔥 RDKit Python Wheels

This repository holds the code to build [RDKit](https://github.com/rdkit/rdkit) platform wheels for Linux, macOS, and Windows. These wheels contain the necessary dynamic libraries (`*.so`, `*.dylib`, and `*.dll`) and can be installed without compilation. The wheels are available at the [PyPi](https://pypi.org/project/rdkit-pypi/) repository and can be installed using pip (`pip install rdkit-pypi`).

Please open an issue if you find something missing or not working as expected.

[![PyPI version shields.io](https://img.shields.io/pypi/v/rdkit-pypi.svg?style=for-the-badge&logo=PyPI&logoColor=blue)](https://pypi.python.org/pypi/rdkit-pypi/)
[![PyPI download month](https://img.shields.io/pypi/dm/rdkit-pypi.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit-pypi/)
[![PyPI download day](https://img.shields.io/pypi/dd/rdkit-pypi.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit-pypi/)
[![GitHub stars](https://img.shields.io/github/stars/kuelumbus/rdkit_platform_wheels.svg?style=social&label=Star&maxAge=2592000)](https://github.com/kuelumbus/rdkit_platform_wheels)

## Available Wheels

| OS      | Arch  | Bit | Conditions                                          | Python              |
| ------- | ----- | --- | --------------------------------------------------- | ------------------- |
| Linux   | intel | 64  | glibc >= 2.17 (e.g., Ubuntu 16.04+, CentOS 6+, ...) | 3.7, 3.8, 3.9, 3.10 |
| macOS   | intel | 64  | >= macOS-11                                         | 3.7, 3.8, 3.9, 3.10 |
| Windows | intel | 64  |                                                     | 3.7, 3.8, 3.9, 3.10 |

## Install RDKit

### PIP

```bash
python -m pip install rdkit-pypi
python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

### [Poetry](https://python-poetry.org/)

```bash
poetry add rdkit-pypi
poetry run python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

## Build wheels locally (Linux only)

cibuildwheel uses `patchelf` (`apt install patchelf`)

```bash
git clone https://github.com/kuelumbus/rdkit_platform_wheels.git
cd rdkit_platform_wheels

python3 -m pip install cibuildwheel
cibuildwheel --platform linux --output-dir wheelhouse
```
