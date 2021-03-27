from setuptools import setup, find_packages
from cmake_build_extension import BuildExtension, CMakeExtension
import os

from pathlib import Path

setup(
    name="rdkit",
    version="2020.9.5",
    description="A collection of chemoinformatics and machine-learning software written in C++ and Python",
    url="http://rdkit.org/",
    license="BSD-3-Clause",
    packages=find_packages("."),
    install_requires=["numpy == 1.20", "pillow"],
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