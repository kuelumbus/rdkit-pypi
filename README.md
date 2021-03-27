# RDKit Python platform wheels

Dockerfile and script to build [RDKit](https://github.com/rdkit/rdkit) python wheels for linux using [manylinux](https://github.com/pypa/manylinux). The script is currently building platform wheels for python 3.7 and 3.8. 

## Install RDKit python package


### Python 3.7 
```bash
pip install https://github.com/kuelumbus/rdkit_platform_wheels/releases/download/2020.9.5/rdkit-2020.9.5-cp37-cp37m-manylinux2014_x86_64.whl
or
poetry add  https://github.com/kuelumbus/rdkit_platform_wheels/releases/download/2020.9.5/rdkit-2020.9.5-cp37-cp37m-manylinux2014_x86_64.whl
```

### Python 3.8
```bash
pip install https://github.com/kuelumbus/rdkit_platform_wheels/releases/download/2020.9.5/rdkit-2020.9.5-cp38-cp38-manylinux2014_x86_64.whl
or
poetry add https://github.com/kuelumbus/rdkit_platform_wheels/releases/download/2020.9.5/rdkit-2020.9.5-cp38-cp38-manylinux2014_x86_64.whl
```

## Build wheels

Clone the repository
```bash
git clone https://github.com/kuelumbus/rdkit_platform_wheels.git
cd rdkit_platform_wheels
```

### `Dockerfile.manylinux2014_x86_64` 
Wheels for glibc >= 2.17 (Ubuntu 16.04 +)

```bash
docker build -f "Dockerfile.manylinux2014_x86_64" -t rdkitdocker:latest .
docker run --rm -e PLAT=manylinux2014_x86_64 -v `pwd`:/io rdkitdocker bash /io/wheeling.sh
```

### `Dockerfile.manylinux_2_24_x86_64`
Wheels for glibc >= 2.24  (Ubuntu 18.04 +)

```bash
docker build -f "Dockerfile.manylinux_2_24_x86_64" -t rdkitdocker:latest .
docker run --rm -e PLAT=manylinux_2_24_x86_64 -v `pwd`:/io rdkitdocker bash /io/wheeling.sh
```


Afterward, you should see a `wheelhouse` directory with the wheels.

