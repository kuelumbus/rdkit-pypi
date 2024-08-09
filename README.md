# 🔥 RDKit Python Wheels

This repository holds the code to build [RDKit](https://github.com/rdkit/rdkit) platform wheels for Linux, macOS, and Windows on Github Action and Circle CI. The wheels contain the compiled platform-specific dynamic libraries (`*.so`, `*.dylib`, and `*.dll`) and are available at [PyPI](https://pypi.org/project/rdkit/). RDKit can easily be installed using

```sh
pip install rdkit
```

**NOTE:** Older versions of RDKit might be available at the [`rdkit-pypi`](https://pypi.org/project/rdkit-pypi/) PyPI repository (`pip install rdkit-pypi`). `rdkit-pypi` is the old name of this project at PyPI. Future RDKit versions will be available at the `rdkit` PyPI repository. Please update your dependencies, i.e., change `rdkit-pypi` to `rdkit`.

Please open an issue if you find something missing or not working as expected.


[![PyPI version shields.io](https://img.shields.io/pypi/v/rdkit.svg?style=for-the-badge&logo=PyPI&logoColor=blue)](https://pypi.python.org/pypi/rdkit/)
[![PyPI download day](https://img.shields.io/pypi/dm/rdkit.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit/)
[![PyPI download month](https://img.shields.io/pypi/dw/rdkit.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit/)
[![PyPI download day](https://img.shields.io/pypi/dd/rdkit.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit/)
![GitHub Repo stars](https://img.shields.io/github/stars/kuelumbus/rdkit-pypi?style=for-the-badge&logo=github)
## Available Builds

| OS      | Arch    | Bit | Conditions                                          | 3.8            | 3.9 | 3.10 | 3.11 | 3.12 | CI                |
| ------- | ------- | --- | --------------------------------------------------- | -------------- | --- | ---- | ---- | ---- | ----------------- |
| Linux   | intel   | 64  | glibc >= 2.28 (e.g., Ubuntu 18.04+, CentOS 6+, ...) | last: 2024.3.5 | ✔️   | ✔️    | ✔️    | ✔️    | Github Actions    |
| Linux   | aarch64 | 64  | glibc >= 2.28 (e.g., Raspberry Pi, ...)             | last: 2024.3.5 | ✔️   | ✔️    | ✔️    | ✔️    | Circle CI         |
| macOS   | intel   | 64  | >= macOS 10.13                                      | last: 2024.3.5 | ✔️   | ✔️    | ✔️    | ✔️    | Github Actions    |
| macOS   | armv8   | 64  | >= macOS 11, M1 hardware                            | last: 2024.3.5 | ✔️   | ✔️    | ✔️    | ✔️    | Github Actions |
| Windows | intel   | 64  |                                                     | last: 2024.3.5 | ✔️   | ✔️    | ✔️    | ✔️    | Github Actions    |

## Installation

### PIP

```bash
python -m pip install rdkit
python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

### [Poetry](https://python-poetry.org/)

```bash
poetry add rdkit
poetry run python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

## Local builds on Linux

`cibuildwheel` requires `patchelf` (`apt install patchelf`)

```bash
python3 -m pip install cibuildwheel

git clone https://github.com/kuelumbus/rdkit-pypi.git
cd rdkit-pypi

CIBW_BUILD=cp39-manylinux_x86_64 python3 -m cibuildwheel --platform linux --output-dir wheelhouse --config-file pyproject.toml
```

Replace `cp39-manylinux_x86_64` with `cp310-manylinux_x86_64`, `cp311-manylinux_x86_64`, or `cp312-manylinux_x86_64` to build for different Python