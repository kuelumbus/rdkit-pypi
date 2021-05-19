from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
from sysconfig import get_paths
import os
from subprocess import check_call, call, run, PIPE
import sys

from pathlib import Path

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()



<<<<<<< HEAD
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


    def get_ext_filename(self, ext_name):
        ext_path = ext_name.split('.')
        return os.path.join(*ext_path)
                
    def create_package(self, ext):
        from distutils.file_util import copy_file
        from shutil import copytree, rmtree

        # copy RDKit package
        rdkit_root = Path(self.build_temp).absolute() / 'rdkit_install/' / 'lib'
        rdkit_pyfiles = list(rdkit_root.glob('python*'))[0] / 'site-packages' / 'rdkit' 

        wheel_path = Path(self.get_ext_fullpath(ext.name)).absolute()
        # remove if exists
        if wheel_path.exists():
            rmtree(str(wheel_path))
        copytree(str(rdkit_pyfiles), str(wheel_path))

        # Copy *.so files to /usr/local/lib
        # auditwheel finds the libs at /usr/local/lib
        libs_rdkit_linux = Path(rdkit_root).glob('*.so*')
        libs_rdkit_macos = Path(rdkit_root).glob('*dylib')
        libs_rdkit = list(libs_rdkit_linux) + list(libs_rdkit_macos)

        libs_boost = Path(self.build_temp).absolute() / 'boost_install' / 'lib'

        libs_boost_linux = libs_boost.glob('*.so*')
        libs_boost_mac = libs_boost.glob('*dylib')
        libs_boost = list(libs_boost_linux) + list(libs_boost_mac)

        [copy_file(i, '/usr/local/lib' ) for i in libs_rdkit]
        [copy_file(i, '/usr/local/lib' ) for i in libs_boost]
    
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
            f'wget {ext.boost_download_url} --no-check-certificate',
            f'tar -xzf {Path(ext.boost_download_url).name}',]

        [check_call(c.split()) for c in cmds]

        # Compile Boost
        os.chdir(Path(ext.boost_download_url).with_suffix('').with_suffix('').name)

        cmds = [
            # for linux
            'ln -fs /opt/python/cp36-cp36m/include/python3.6m /opt/python/cp36-cp36m/include/python3.6',
            'ln -fs /opt/python/cp37-cp37m/include/python3.7m /opt/python/cp37-cp37m/include/python3.7',
            # same for MacOS
            'ln -fs /Library/Frameworks/Python.framework/Versions/3.6/include/python3.6m /Library/Frameworks/Python.framework/Versions/3.6/include/python3.6',
            'ln -fs /Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m /Library/Frameworks/Python.framework/Versions/3.7/include/python3.7',
            ]
        
        # Ok to fail
        [call(c.split()) for c in cmds]


        cmds = [
            f'./bootstrap.sh --with-libraries=python,serialization,iostreams,system,regex --with-python={sys.executable} --with-python-root={Path(sys.executable).parent}/..',
            f'./b2 install --prefix={boost_root} -j 20',
        ]
        [check_call(c.split()) for c in cmds]

        os.chdir(str(cwd))

        libs_boost = Path(self.build_temp).absolute() / 'boost_install' / 'lib'
        libs_boost_mac = libs_boost.glob('*dylib')

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

        # Invoke cmake and compile RDKit
        options = [
                    f'-DPYTHON_EXECUTABLE={sys.executable}',
                    f'-DPYTHON_INCLUDE_DIR={get_paths()["include"]}',

                    # RDKIT build flags
                    f"-DRDK_BUILD_INCHI_SUPPORT=ON",
                    f"-DRDK_BUILD_AVALON_SUPPORT=ON",
                    f"-DRDK_BUILD_PYTHON_WRAPPERS=ON",
                    f"-DRDK_INSTALL_INTREE=OFF",
                    f"-DRDK_BUILD_CAIRO_SUPPORT=ON",

                    f"-DBOOST_ROOT={boost_root}",
                    f"-DBoost_NO_SYSTEM_PATHS=ON",

                    f"-DCMAKE_INSTALL_PREFIX={rdkit_root}",
                    f"-DCMAKE_C_FLAGS=-Wno-implicit-function-declaration",
                    f"-DCMAKE_CXX_FLAGS=-Wno-implicit-function-declaration",
                ]

        cmds = [
            f"cmake -S . -B build {' '.join(options)} ",
            f"cmake --build build --verbose -j 5 --config Release",
            f"make -C build install -j 5"
        ]    
        [check_call(c.split()) for c in cmds]
        os.chdir(str(cwd))


=======
>>>>>>> 6722637f123d34574cb136099da4771e15150735
setup(
    name="rdkit-pypi",
    version=f"2021.3.1.5",
    description="A collection of chemoinformatics and machine-learning software written in C++ and Python",
    url="https://github.com/kuelumbus/rdkit_platform_wheels",
    project_urls={
        "RDKit": "http://rdkit.org/",
        "RDKit Github": "https://github.com/rdkit/rdkit",
    },
    license="BSD-3-Clause",
<<<<<<< HEAD
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    ext_modules=[
        RDKit(
            'rdkit',
            boost_download_url='https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.gz',
            rdkit_tag='Release_2021_03_1'
            ),        
=======
    packages=find_packages("."),
    # install_requires=["numpy"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
    # Include all libs and files in /rdkit/rdkit/* 
    package_data={'': [str(i) for i in Path('/rdkit/rdkit').absolute().rglob('*')] },
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
                            f"-DPYTHON_NUMPY_INCLUDE_PATH={os.environ['NUMPY_INC']}",
                            f"-DRDK_BUILD_CAIRO_SUPPORT=ON",
                       ]),
>>>>>>> 6722637f123d34574cb136099da4771e15150735
    ],
    cmdclass=dict(build_ext=BuildRDKit),

    
)