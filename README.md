# 🔥 RDKit Python Wheels

This repository holds the code to build [RDKit](https://github.com/rdkit/rdkit) platform wheels for Linux, macOS, and Windows. The wheels contain the platform-specific dynamic libraries (`*.so`, `*.dylib`, and `*.dll`) and are available at [PyPi](https://pypi.org/project/rdkit-pypi/). You can install RDKit using pip (`pip install rdkit-pypi`).


Please open an issue if you find something missing or not working as expected.

[![PyPI version shields.io](https://img.shields.io/pypi/v/rdkit-pypi.svg?style=for-the-badge&logo=PyPI&logoColor=blue)](https://pypi.python.org/pypi/rdkit-pypi/)
[![PyPI download month](https://img.shields.io/pypi/dm/rdkit-pypi.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit-pypi/)
[![PyPI download day](https://img.shields.io/pypi/dd/rdkit-pypi.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit-pypi/)
[![GitHub stars](https://img.shields.io/github/stars/kuelumbus/rdkit-pypi.svg?style=social&label=Star&maxAge=2592000)](https://github.com/kuelumbus/rdkit-pypi)

## Available Builds

| OS      | Arch    | Bit | Conditions                                          | 3.7 | 3.8 | 3.9 | 3.10 | CI             |
| ------- | ------- | --- | --------------------------------------------------- | --- | --- | --- | ---- | -------------- |
| Linux   | intel   | 64  | glibc >= 2.17 (e.g., Ubuntu 16.04+, CentOS 6+, ...) | ✔️  | ✔️  | ✔️  | ✔️   | Github Actions |
| Linux   | aarch64 | 64  | glibc >= 2.17 (e.g., Raspberry Pi, ...)             | ✔️  | ✔️  | ✔️  | ✔️   | Circle CI      |
| macOS   | intel   | 64  | >= macOS-11                                         | ✔️  | ✔️  | ✔️  | ✔️   | Github Actions |
| macOS   | armv8   | 64  | >= macOS-11 (M1 hardware)                           |     | ✔️  | ✔️  | ✔️   | Github Actions |
| Windows | intel   | 64  |                                                     | ✔️  | ✔️  | ✔️  | ✔️   | Github Actions |

## Install RDKit from PyPi

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

## Local builds on Linux

`cibuildwheel` requires `patchelf` (`apt install patchelf`)

```bash
python3 -m pip install cibuildwheel

git clone https://github.com/kuelumbus/rdkit-pypi.git
cd rdkit-pypi

CIBW_BUILD=cp37-manylinux_x86_64 python3 -m cibuildwheel --platform linux --output-dir wheelhouse --config-file pyproject.toml
```

Replace `cp37-manylinux_x86_64` with `cp38-manylinux_x86_64`, `cp39-manylinux_x86_64`, and `cp310-manylinux_x86_64` to build for other python versions.
