# 🔥 RDKit Python Wheels

This repository holds the code to build [RDKit](https://github.com/rdkit/rdkit) platform wheels for Linux, macOS, and Windows on Github Action and Circle CI. The wheels contain the compiled platform-specific dynamic libraries (`*.so`, `*.dylib`, and `*.dll`) and are available at [PyPI](https://pypi.org/project/rdkit/). RDKit can easily be installed using

```sh
pip install rdkit
```

Please open an issue if you find something missing or not working as expected.


[![PyPI version shields.io](https://img.shields.io/pypi/v/rdkit.svg?style=for-the-badge&logo=PyPI&logoColor=blue)](https://pypi.python.org/pypi/rdkit/)
[![PyPI download day](https://img.shields.io/pypi/dm/rdkit.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit/)
[![PyPI download month](https://img.shields.io/pypi/dw/rdkit.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit/)
[![PyPI download day](https://img.shields.io/pypi/dd/rdkit.svg?style=for-the-badge&logo=PyPI)](https://pypi.python.org/pypi/rdkit/)
![GitHub Repo stars](https://img.shields.io/github/stars/kuelumbus/rdkit-pypi?style=for-the-badge&logo=github)
## Available Builds

| OS      | Arch    | Bit | Conditions                                          | 3.8            | 3.9 | 3.10 | 3.11 | 3.12 | 3.13 | 3.14 | CI             |
| ------- | ------- | --- | --------------------------------------------------- | -------------- | --- | ---- | ---- | ---- | ---- | ---- | -------------- |
| Linux   | intel   | 64  | glibc >= 2.28 (e.g., Ubuntu 18.04+, CentOS 6+, ...) | last: 2024.3.5 | last: 2025.9.2    | ✔️    | ✔️    | ✔️    | ✔️    | ✔️    | Github Actions |
| Linux   | aarch64 | 64  | glibc >= 2.28 (e.g., Raspberry Pi, ...)             | last: 2024.3.5 | last: 2025.9.2   | ✔️    | ✔️    | ✔️    | ✔️    | ✔️    | Circle CI      |
| macOS   | intel   | 64  | >= macOS 10.15                                      | last: 2024.3.5 | last: 2025.9.2   | last: 2025.9.3    | last: 2025.9.3    | last: 2025.9.3    | last: 2025.9.3    | last: 2025.9.3    | Github Actions |
| macOS   | armv8   | 64  | >= macOS 11, M-Hardware                             | last: 2024.3.5 | last: 2025.9.2   | ✔️    | ✔️    | ✔️    | ✔️    | ✔️    | Github Actions | 
| Windows | intel   | 64  |                                                     | last: 2024.3.5 | last: 2025.9.2   | ✔️    | ✔️    | ✔️    | ✔️    | ✔️    | Github Actions |

## Installation

### PIP

```bash
python -m pip install rdkit
python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

### [uv](https://docs.astral.sh/uv/getting-started/installation/)

```bash
uv add rdkit
uv run python -c "from rdkit import Chem; print(Chem.MolToMolBlock(Chem.MolFromSmiles('C1CCC1')))"
```

## Local builds on Linux

`cibuildwheel` requires `patchelf` (`apt install patchelf`)

```bash
python3 -m pip install cibuildwheel

git clone https://github.com/kuelumbus/rdkit-pypi.git
cd rdkit-pypi

CIBW_BUILD=cp313-manylinux_x86_64 python3 -m cibuildwheel --platform linux --output-dir wheelhouse --config-file pyproject.toml
```

Replace `*` in `cp*-manylinux_x86_64` with `310`, `311`, `312`, `313`, or `314` to build for different Python versions.

## Local builds on macOS

### Prerequisites

Install the required dependencies via Homebrew:

```bash
brew install lapack eigen gnu-sed
```

### Building with pip wheel

For local macOS builds, use `pip wheel` directly (avoids cibuildwheel's requirement for system Python):

```bash
git clone https://github.com/kuelumbus/rdkit-pypi.git
cd rdkit-pypi

# For macOS ARM64 (Apple Silicon M1/M2/M3)
RDKIT_OSMORDRED_VERSION=2025-9-4-v2 \
CIBW_BUILD=cp311-macosx_arm64 \
pip wheel . --no-deps -w wheelhouse -v

# For macOS x86_64 (Intel)
RDKIT_OSMORDRED_VERSION=2025-9-4-v2 \
CIBW_BUILD=cp311-macosx_x86_64 \
pip wheel . --no-deps -w wheelhouse -v
```

### macOS SDK Compatibility (RDKIT_CATCH2_LEGACY)

If you encounter a compilation error like:

```
error: no member named 'uncaught_exception' in namespace 'std'; did you mean 'uncaught_exceptions'?
```

This is a C++ compatibility issue between different macOS SDK versions:

- **Modern macOS (SDK 26+)**: Uses `std::uncaught_exceptions()` (C++17/20 standard)
- **Older macOS (SDK < 26)**: Uses `std::uncaught_exception()` (deprecated, removed in C++20)

Use the `RDKIT_CATCH2_LEGACY` environment variable to select the appropriate mode:

```bash
# For older macOS (SDK < 26) - DEFAULT
# Legacy mode is enabled by default for maximum compatibility
RDKIT_OSMORDRED_VERSION=2025-9-4-v2 \
CIBW_BUILD=cp311-macosx_arm64 \
pip wheel . --no-deps -w wheelhouse -v

# For modern macOS (SDK 26+) - disable legacy mode
RDKIT_OSMORDRED_VERSION=2025-9-4-v2 \
RDKIT_CATCH2_LEGACY=0 \
CIBW_BUILD=cp311-macosx_arm64 \
pip wheel . --no-deps -w wheelhouse -v
```

The build output will show which mode is being used:
- `Catch2 Legacy: True` = legacy mode (default, for older macOS)
- `Catch2 Legacy: False` = modern mode for macOS SDK 26+

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `RDKIT_OSMORDRED_VERSION` | Version string in format `YYYY-M-P-vN` (e.g., `2025-9-4-v2`) | `2025-9-5-v3` |
| `CIBW_BUILD` | Target platform (e.g., `cp311-macosx_arm64`, `cp311-manylinux_x86_64`) | Required |
| `RDKIT_CATCH2_LEGACY` | Set to `1` for older macOS SDK compatibility | `0` |