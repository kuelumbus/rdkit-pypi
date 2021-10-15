from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
from sysconfig import get_paths
import os
from subprocess import check_call, call, run, PIPE
import sys
from sys import platform
import shutil

from pathlib import Path

# get vcpkg path on Github
vcpkg_path = Path('C:/vcpkg')


def towin(pt: Path):
    """Returns a windows path from a Path object"""
    return str(pt).replace('\\', '/')

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
            # Build boot
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

        # copy RDKit
        if platform == 'win32':
            rdkit_root = Path(self.build_temp).absolute() / 'rdkit_install/' / 'Lib'
            rdkit_pyfiles = rdkit_root / 'site-packages' / 'rdkit' 
        else:
            rdkit_root = Path(self.build_temp).absolute() / 'rdkit_install/' / 'lib'
            rdkit_pyfiles = list(rdkit_root.glob('python*'))[0] / 'site-packages' / 'rdkit' 

        # rdkit needs some files from the Data directory to run correctly 
        rdkit_data_path = Path(self.build_temp).absolute() / 'rdkit' / 'Data'

        # replace line with _share=... with _share = os.path.dirname(__file__) in RDPaths.py
        rdpaths = rdkit_pyfiles / 'RDPaths.py'

        # For linux
        if platform == "linux" or platform == "linux2" or platform == 'win32':
            # linux
            call(["sed", "-i", "/_share =/c\_share = os.path.dirname(__file__)", str(rdpaths)])
        elif platform == "darwin":
            # OS X
            call(["gsed", "-i", "/_share =/c\_share = os.path.dirname(__file__)", str(rdpaths)])

        wheel_path = Path(self.get_ext_fullpath(ext.name)).absolute()
        # remove if exists
        if wheel_path.exists():
            rmtree(str(wheel_path))

        # copy rdkit files 
        copytree(str(rdkit_pyfiles), str(wheel_path))

        # copy the Data directory to the wheel path
        copytree(str(rdkit_data_path), str(wheel_path / 'Data'))

        # Copy *.so files to /usr/local/lib
        # auditwheel finds the libs at /usr/local/lib
        libs_rdkit_linux = Path(rdkit_root).glob('*.so*')
        libs_rdkit_macos = Path(rdkit_root).glob('*dylib')
        libs_rdkit = list(libs_rdkit_linux) + list(libs_rdkit_macos)

        libs_boost = Path(self.build_temp).absolute() / 'boost_install' / 'lib'
        libs_boost_linux = libs_boost.glob('*.so*')
        libs_boost_mac = libs_boost.glob('*dylib')
        libs_boost_tmp = list(libs_boost_linux) + list(libs_boost_mac)

        if platform != 'win32':
            [copy_file(i, '/usr/local/lib' ) for i in libs_rdkit]
            [copy_file(i, '/usr/local/lib' ) for i in libs_boost_tmp]
        else:
            libs_rdkit_win = Path(rdkit_root).glob('*.dll')
            libs_boost_win = libs_boost.glob('*.dll')
            
            libs_vcpkg = list((vcpkg_path / 'installed' / 'x64-windows' / 'bin').glob('*.dll'))

            [copy_file(i, 'C://libs' ) for i in libs_rdkit_win]
            [copy_file(i, 'C://libs' ) for i in libs_boost_win]
            [copy_file(i, 'C://libs' ) for i in libs_vcpkg]
    
    def build_boost(self, ext):
        """Build the Boost libraries"""
        
        cwd = Path().absolute()
        boost_build_path = Path(self.build_temp).absolute() / 'boost'
        boost_build_path.mkdir(parents=True, exist_ok=True)

        boost_install_path = Path(self.build_temp).absolute() / 'boost_install'
        boost_install_path.mkdir(parents=True, exist_ok=True)

        # Download and unpack Boost
        os.chdir(str(boost_build_path))

        cmds = [
            f'wget {ext.boost_download_url} --no-check-certificate -q',
            f'tar -xzf {Path(ext.boost_download_url).name}',
            ]

        [check_call(c.split()) for c in cmds]

        # Compile Boost
        os.chdir(Path(ext.boost_download_url).with_suffix('').with_suffix('').name)

        # This fixes a bug in the boost configuration. Boost expects python include paths without "m"
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
         
   
        # Change commands for windows
        if sys.platform == 'win32':

            cmds = [
                f'bootstrap.bat',
            ]
            [check_call(c.split()) for c in cmds]

            # Compile for many python versions at the same time?
            python_inc = Path(get_paths()["include"])
            python_libs = Path(get_paths()["data"]) / 'libs'
            python_exe = Path(get_paths()["data"]) / 'python.exe'
            
            zlib_include = vcpkg_path / 'packages/zlib_x64-windows/include'
            zlib_lib = vcpkg_path / 'packages/zlib_x64-windows/lib'

            bzip2_include = vcpkg_path / 'packages/bzip2_x64-windows/include'
            bzip2_lib = vcpkg_path / 'packages/bzip2_x64-windows/lib'
            with open('project-config.jam', 'a') as fl:
                print(f'using python : {sys.version_info[0]}.{sys.version_info[1]} : {towin(python_exe)} : {towin(python_inc)} : {towin(python_libs)} ;', file=fl)
                print(f' ', file=fl)
                print(f'using zlib : 2 : <include>{towin(zlib_include)} <search>{towin(zlib_lib)} ;', file=fl)
                print(f' ', file=fl)
                print(f'using bzip2 : 1 : <include>{towin(bzip2_include)} <search>{towin(bzip2_lib)} ;', file=fl)
                print(f' ', file=fl)
            
            cmds = [                
                f'./b2 address-model=64 architecture=x86 link=static link=shared threading=single threading=multi ' \
                f'variant=release -d0 --abbreviate-paths ' \
                f'--with-python --with-serialization --with-iostreams --with-system --with-regex --with-program_options ' \
                f'--prefix={boost_install_path} -j 20 install',
            ]
            [check_call(c.split()) for c in cmds]

        else:
            cmds = [
            f'./bootstrap.sh --with-libraries=python,serialization,iostreams,system,regex --with-python={sys.executable} --with-python-root={Path(sys.executable).parent}/..',
            f'./b2 install --prefix={boost_install_path} -j 20',
            ]
            [check_call(c.split()) for c in cmds]

            
        check_call(['ls', towin(boost_install_path / 'lib')])

        os.chdir(str(cwd))

    def build_rdkit(self, ext):
        """ Build RDKit """

        cwd = Path().absolute()
        rdkit_build_path = Path(self.build_temp).absolute() 
        rdkit_build_path.mkdir(parents=True, exist_ok=True) 

        boost_install_path = Path(self.build_temp).absolute() / 'boost_install'

        rdkit_install_path = Path(self.build_temp).absolute() / 'rdkit_install'
        rdkit_install_path.mkdir(parents=True, exist_ok=True)

        # Clone RDKit from git at rdkit_tag
        os.chdir(str(rdkit_build_path))
        cmds = [
            f'git clone -b {ext.rdkit_tag} https://github.com/rdkit/rdkit'
            ]
        [check_call(c.split()) for c in cmds]
        
        os.chdir(str('rdkit'))
        
        # all includes are here
        vcpkg_install_path = vcpkg_path / 'installed' / 'x64-windows' 
        vcpkg_include_path = vcpkg_path / 'installed' / 'x64-windows' / 'include'
        vcpkg_lib_path = vcpkg_path / 'installed' / 'x64-windows' / 'lib'

        # Invoke cmake and compile RDKit
        options = [ 
            # Defines the paths to many include and libaray paths for windows
            # Does not work for some reason??
            # f"-DCMAKE_TOOLCHAIN_FILE=C:\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake" if sys.platform == 'win32' else "",
            # f"-DVCPKG_TARGET_TRIPLET=x64-windows-static" if sys.platform == 'win32' else "",

            f'-DPYTHON_EXECUTABLE={sys.executable}',
            f'-DPYTHON_INCLUDE_DIR={get_paths()["include"]}',

            # RDKIT build flags
            f"-DRDK_BUILD_INCHI_SUPPORT=ON",
            f"-DRDK_BUILD_AVALON_SUPPORT=ON",
            f"-DRDK_BUILD_PYTHON_WRAPPERS=ON",
            f"-DRDK_INSTALL_INTREE=OFF",

            # Boost              
            f"-DBOOST_ROOT={boost_install_path}",
            f"-DBoost_NO_SYSTEM_PATHS=ON",
            f"-DBoost_DEBUG=OFF",        

            # Does not work (this is fixed in future rdkit versions I believe)
            f"-DRDK_INSTALL_STATIC_LIBS=OFF" if sys.platform == 'win32' else "",

            # ##### for windows 
            # cairo
            f"-DRDK_BUILD_CAIRO_SUPPORT=ON",
            f"-DCAIRO_INCLUDE_DIR={towin(vcpkg_include_path)}" if sys.platform == 'win32' else "",
            f"-DCAIRO_LIBRARY_DIR={towin(vcpkg_lib_path)}" if sys.platform == 'win32' else "",
            
            # zlib
            f"-DZLIB_ROOT={towin(vcpkg_install_path)}" if sys.platform == 'win32' else "",

            # freetype
            f"-DFREETYPE_INCLUDE_DIRS={towin(vcpkg_include_path)}" if sys.platform == 'win32' else "",
            f"-DFREETYPE_LIBRARY={towin(vcpkg_lib_path / 'freetype.lib')}" if sys.platform == 'win32' else "",

            # eigen3
            f"-DEIGEN3_INCLUDE_DIR={towin(vcpkg_include_path)}" if sys.platform == 'win32' else "",

            # instruct to build x64
            "-Ax64" if sys.platform == 'win32' else "",

            # Mac needs these flags to compile 
            f"-DCMAKE_C_FLAGS=-Wno-implicit-function-declaration" if sys.platform == 'darwin' else "", 
            f"-DCMAKE_CXX_FLAGS=-Wno-implicit-function-declaration" if sys.platform == 'darwin' else "", 

            # build stuff
            f"-DCMAKE_INSTALL_PREFIX={rdkit_install_path}",
            f"-DCMAKE_BUILD_TYPE=Release",
        ]
        
        cmds = [
            f"cmake -S . -B build {' '.join(options)} ",
            f"cmake --build build -j 10 --config Release",
            f"cmake --install build"
        ]    
        [check_call(c.split()) for c in cmds]
        os.chdir(str(cwd))



setup(
    name="rdkit-pypi",
    version=f"2021.3.5.1",
    description="A collection of chemoinformatics and machine-learning software written in C++ and Python",
    author='Christopher Kuenneth',
    author_email='chris@kuenneth.dev',
    url="https://github.com/kuelumbus/rdkit_platform_wheels",
    project_urls={
        "RDKit": "http://rdkit.org/",
        "RDKit on Github": "https://github.com/rdkit/rdkit",
    },
    license="BSD-3-Clause",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=[
          'numpy>=1.19',
      ],
    ext_modules=[
        RDKit(
            'rdkit',
            # 1.73 does now compile on win for some reason
            boost_download_url='https://boostorg.jfrog.io/artifactory/main/release/1.69.0/source/boost_1_69_0.tar.gz' if sys.platform != 'win32' else 'https://boostorg.jfrog.io/artifactory/main/release/1.69.0/source/boost_1_69_0.tar.gz',
            rdkit_tag='Release_2021_03_5'
            ),        
    ],
    cmdclass=dict(build_ext=BuildRDKit),
)