import imp
from shutil import copyfile
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

import os
from subprocess import check_call
import sys

from pathlib import Path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()



class RDKit(Extension):

    def __init__(self, name, **kwargs):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])
        self.__dict__.update(kwargs)


class BuildRDKit(build_ext_orig):

    def run(self):
        for ext in self.extensions:
            # Build boost
            self.build_boost(ext)
            # Then RDKit
            self.build_rdkit(ext)
            # Copy files so that a wheels package can be created
            self.create_package(ext)
            
        # invoke the creation of the wheels package
        super().run()

    def create_package(self, ext):
        from distutils.file_util import copy_file
        from shutil import copytree, copy
        wheel_path = Path(self.get_ext_fullpath(ext.name)).absolute()
        wheel_path.mkdir(parents=True, exist_ok=True)
        rdkit_root = Path(self.build_temp).absolute() / 'rdkit_install/' / 'lib'

        # copy rkdit
        rdkit_pyfiles = list(rdkit_root.glob('python*'))[0] / 'site-packages' 
        copytree(str(rdkit_pyfiles), wheel_path)

        # copy *.so
        libs = wheel_path / 'rdkit' / 'libs'
        libs.mkdir(parents=True, exist_ok=True)
        [copy_file(i, libs ) for i in Path(rdkit_root).glob('*.so')]
    
    def build_boost(self, ext):
        # Build boost libraries
        cwd = Path().absolute()

        build_temp = Path(self.build_temp).absolute() / 'boost'
        build_temp.mkdir(parents=True, exist_ok=True)

        boost_root = Path(self.build_temp).absolute() / 'boost_install'
        boost_root.mkdir(parents=True, exist_ok=True)

        # Download Boost
        os.chdir(str(build_temp))
        cmds = [
            f'wget {ext.boost_download_url}',
            f'tar -xzf {Path(ext.boost_download_url).name}',]

        [check_call(c.split()) for c in cmds]

        # Compile Boost
        os.chdir(Path(ext.boost_download_url).with_suffix('').with_suffix('').name)
        cmds = [
            'ln -fs  /opt/python/cp36-cp36m/include/python3.6m /opt/python/cp36-cp36m/include/python3.6',
            'ln -fs  /opt/python/cp37-cp37m/include/python3.7m /opt/python/cp37-cp37m/include/python3.7',
            f'./bootstrap.sh --with-libraries=python,serialization,iostreams,system,regex --with-python={sys.executable} --with-python-root={Path(sys.executable).parent}/..',
            f'./b2 install --prefix={boost_root} -j 20',
        ]
        [check_call(c.split()) for c in cmds]
        os.chdir(str(cwd))

    def build_rdkit(self, ext):
        # Build RDKit
        cwd = Path().absolute()
        build_temp = Path(self.build_temp).absolute() / 'rdkit'
        build_temp.mkdir(parents=True, exist_ok=True) 

        boost_root = Path(self.build_temp).absolute() / 'boost_install'

        rdkit_root = Path(self.build_temp).absolute() / 'rdkit_install'
        rdkit_root.mkdir(parents=True, exist_ok=True)

        # Clone from GIT
        os.chdir(str(build_temp))
        cmds = [
            f'git clone -b {ext.rdkit_tag} https://github.com/rdkit/rdkit'
            ]
        [check_call(c.split()) for c in cmds]
        os.chdir(str('rdkit'))

        # Get the Python include path (is there a better way?)
        pyinc_path = Path(sys.executable).parent / ".." / "include"
        pyinc_path = list(pyinc_path.glob('python*'))[0]

        # Invoke cmake and compile RDKit
        options = [
                    f'-DPYTHON_EXECUTABLE={sys.executable}',
                    f'-DPYTHON_INCLUDE_DIR={pyinc_path}',
                    f"-DRDK_BUILD_INCHI_SUPPORT=ON",
                    f"-DRDK_BUILD_AVALON_SUPPORT=ON",
                    f"-DRDK_BUILD_PYTHON_WRAPPERS=ON",
                    # f"-DRDK_BUILD_CAIRO_SUPPORT=ON",
                    f"-DBoost_NO_SYSTEM_PATHS=ON",
                    f"-DBOOST_ROOT={boost_root}",
                    f"-DRDK_INSTALL_INTREE=off",
                    f"-DCMAKE_INSTALL_PREFIX={rdkit_root}"
                ]

        cmds = [
            f"cmake -S . -B build {' '.join(options)} ",
            "cmake --build build --verbose -j 20 --config Release",
            "make -C build install -j 20"
        ]    
        [check_call(c.split()) for c in cmds]
        os.chdir(str(cwd))

            

setup(
    name="rdkit-pypi",
    version=f"2021.03.01",
    description="A collection of chemoinformatics and machine-learning software written in C++ and Python",
    url="https://github.com/kuelumbus/rdkit_platform_wheels",
    project_urls={
        "RDKit": "http://rdkit.org/",
        "RDKit Github": "https://github.com/rdkit/rdkit",
    },
    license="BSD-3-Clause",
    install_requires=["numpy == 1.20", "pillow"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[
        RDKit(
            'rdkit',
            boost_download_url='https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz',
            rdkit_tag='Release_2021_03_1'
            ),        
    ],
    cmdclass=dict(build_ext=BuildRDKit),

    
)