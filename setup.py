from setuptools import setup, find_packages
from cmake_build_extension import BuildExtension, CMakeExtension
import os

from pathlib import Path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="rdkit-pypi",
    version=f"{os.environ['RDKIT_VERSION']}",
    description="A collection of chemoinformatics and machine-learning software written in C++ and Python",
    url="https://github.com/kuelumbus/rdkit_platform_wheels",
    project_urls={
        "RDKit": "http://rdkit.org/",
        "RDKit Github": "https://github.com/rdkit/rdkit",
    },
    license="BSD-3-Clause",
    packages=find_packages("."),
    install_requires=["numpy == 1.20", "pillow"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
    # Include all libs and files in /rdkit/rdkit/* 
    package_data={'': ['/rdkit/lib/*.so' ] + [str(i) for i in Path('/rdkit/rdkit').absolute().rglob('*')] },
    ext_modules=[
        CMakeExtension(name="RdkitCompile",
                       install_prefix="rdkit",
                       cmake_configure_options=[
                            f'-DPYTHON_EXECUTABLE={os.environ["PYBIN"]}/python',
                            f'-DPYTHON_INCLUDE_DIR={os.environ["PYINC"]}',
                            f"-DRDK_BUILD_INCHI_SUPPORT=ON",
                            f"-DRDK_BUILD_AVALON_SUPPORT=ON",
                            f"-DRDK_BUILD_PYTHON_WRAPPERS=ON",
                            f"-DBOOST_ROOT={os.environ['BOOST_ROOT']}",
                       ]),
    ],
    cmdclass=dict(build_ext=BuildExtension),
)