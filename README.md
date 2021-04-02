# RDKit Python platform wheels

Dockerfile and script to build [RDKit](https://github.com/rdkit/rdkit) Python wheels for Linux using [manylinux](https://github.com/pypa/manylinux). 

Versions:

- Linux: Python >= 3.6 and glibc >= 2.17 (e.g., Ubuntu 16.04+, CentOS 6+, ...)


## Install RDKit using pip

```bash
pip install rdkit-pypi
python -c "from rdkit import Chem"
```

## Build wheels locally (works only for Linux) 

cibuildwheel uses `patchelf` (`apt install patchelf`) 

```bash
git clone https://github.com/kuelumbus/rdkit_platform_wheels.git
cd rdkit_platform_wheels

python3.8 -m pip install cibuildwheel

CIBW_BUILD_VERBOSITY=1 CIBW_MANYLINUX_X86_64_IMAGE=manylinux2014 CIBW_BEFORE_BUILD_LINUX="bash pre_linux.sh" cibuildwheel --platform linux --output-dir wheelhouse
```