# RDKit Python platform wheels

Use [cibuildwheel](https://github.com/joerick/cibuildwheel) and Github Actions to build [RDKit](https://github.com/rdkit/rdkit) wheels for Linux and MacOS. Wheels are available on PyPi.

Versions:

- Linux: 3.6 <= Python <= 3.9 and glibc >= 2.17 (e.g., Ubuntu 16.04+, CentOS 6+, ...)
- MacOS 10.9: 3.6 <= Python <= 3.9 

## Install RDKit using pip

```bash
pip install rdkit-pypi
python -c "from rdkit import Chem"
```

## Install RDKit using poetry
```bash
poetry add rdkit-pypi
```

## Build wheels locally (Linux only)

cibuildwheel uses `patchelf` (`apt install patchelf`) 

```bash
git clone https://github.com/kuelumbus/rdkit_platform_wheels.git
cd rdkit_platform_wheels

python3.8 -m pip install cibuildwheel

CIBW_BUILD_VERBOSITY=1 CIBW_MANYLINUX_X86_64_IMAGE=manylinux2014 CIBW_BEFORE_BUILD_LINUX="bash pre_linux.sh" cibuildwheel --platform linux --output-dir wheelhouse
```