# RDKit Python platform wheels

Dockerfile and script to build [RDKit](https://github.com/rdkit/rdkit) Python wheels for Linux using [manylinux](https://github.com/pypa/manylinux). 

Versions:

- Linux: Python >= 3.6 and glibc >= 2.17 (e.g., Ubuntu 16.04+, CentOS 6+, ...)


## Install RDKit using pip

```bash
pip install rdkit-pypi
python -c "from rdkit import Chem"
```

## Build wheels 

Clone the repository
```bash
git clone https://github.com/kuelumbus/rdkit_platform_wheels.git
cd rdkit_platform_wheels
```

- RDKIT_RELEASE (tag from [RDKit](https://github.com/rdkit/rdkit))
- RDKIT_VERSION (format must be yyyy.m.d)

### `Dockerfile.manylinux2014_x86_64` 
Wheels for glibc >= 2.17 (Ubuntu 16.04 +)

```bash
docker build -f "Dockerfile.manylinux2014_x86_64" -t rdkitdocker:latest .
docker run --rm -e PLAT=manylinux2014_x86_64 -e RDKIT_RELEASE=Release_2021_03_1 -e RDKIT_VERSION=2021.03.1 -v `pwd`:/io rdkitdocker bash /io/wheeling.sh
```

### `Dockerfile.manylinux_2_24_x86_64`
Wheels for glibc >= 2.24  (Ubuntu 18.04 +)

```bash
docker build -f "Dockerfile.manylinux_2_24_x86_64" -t rdkitdocker:latest .
docker run --rm -e PLAT=manylinux_2_24_x86_64 -e RDKIT_RELEASE=Release_2020_09_5 -e RDKIT_VERSION=2020.9.5 -v `pwd`:/io rdkitdocker bash /io/wheeling.sh
```

Afterward, you should see a `wheelhouse` directory with the wheels.

