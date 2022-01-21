from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
from sysconfig import get_paths
import os
from subprocess import check_call, call, run, PIPE
import sys
from sys import platform
from distutils.file_util import copy_file
from shutil import copytree, rmtree
from pathlib import Path
from textwrap import dedent

# get vcpkg path on Github
vcpkg_path = Path("C:/vcpkg")


def towin(pt: Path):
    """Returns a windows path from a Path object"""
    return str(pt).replace("\\", "/")


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
            self.build_rdkit(ext)

        # invoke the creation of the wheels package
        super().run()

    def get_ext_filename(self, ext_name):
        ext_path = ext_name.split(".")
        return os.path.join(*ext_path)

    def conan_install(self, conan_toolchain_path):
        """Run the Conan build"""
        boost_version = "1.77.0"

        # The modified conanfile.py from boost does not link libpython*.so
        # When building a platform wheel, we don't want to link libpython*.so.
        mod_conan_path = "conan_boost_mod"

        # Export the modified boost version
        check_call(
            [
                f"conan",
                f"export",
                f"{mod_conan_path}/all/",
                f"{boost_version}@chris/mod_boost",
            ]
        )

        # needed for windows builds
        win = """cairo/1.17.4
            freetype/2.11.1
            eigen/3.4.0
        """
        if platform != "win32":
            win = ""

        conanfile = f"""\
            [requires]
            boost/{boost_version}@chris/mod_boost
            {win}

            [generators]
            cmake_paths
            virtualrunenv

            [options]
            boost:shared=True
            boost:without_python=False
            boost:without_python_lib=True
            boost:python_executable={sys.executable}
        """

        Path("conanfile.txt").write_text(dedent(conanfile))

        # run conan install
        cmd = [
            f"conan",
            f"install",
            f"conanfile.txt",
            # build all missing
            f"--build=missing",
            f"-if",
            f"{conan_toolchain_path}",
        ]

        # but force build b2 on linux
        if "linux" in sys.platform:
            cmd += [f"--build=b2"]
        # if "win" in sys.platform:
        #     cmd += ["-s", "compiler.cppstd=20"]

        check_call(cmd)

    def build_rdkit(self, ext):
        """Build RDKit"""

        cwd = Path().absolute()

        conan_toolchain_path = cwd / "conan"
        conan_toolchain_path.mkdir(parents=True, exist_ok=True)

        # Install boost via conan
        self.conan_install(conan_toolchain_path)

        build_path = Path(self.build_temp).absolute()
        build_path.mkdir(parents=True, exist_ok=True)
        os.chdir(str(build_path))

        rdkit_install_path = Path(self.build_temp).absolute() / "rdkit_install"
        rdkit_install_path.mkdir(parents=True, exist_ok=True)

        # Clone RDKit from git at rdkit_tag
        cmds = [f"git clone https://github.com/rdkit/rdkit"]
        [check_call(c.split()) for c in cmds]

        os.chdir(str("rdkit"))

        # all includes are here
        vcpkg_install_path = vcpkg_path / "installed" / "x64-windows"
        vcpkg_include_path = vcpkg_path / "installed" / "x64-windows" / "include"
        vcpkg_lib_path = vcpkg_path / "installed" / "x64-windows" / "lib"

        # Invoke cmake and compile RDKit
        options = [
            # Defines the paths to many include and library paths for windows
            # Does not work for some reason??
            f"-DCMAKE_TOOLCHAIN_FILE={conan_toolchain_path / 'conan_paths.cmake'}",
            # f"-DVCPKG_TARGET_TRIPLET=x64-windows-static" if sys.platform == 'win32' else "",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f'-DPYTHON_INCLUDE_DIR={get_paths()["include"]}',
            # RDKIT build flags
            f"-DRDK_BUILD_INCHI_SUPPORT=ON",
            f"-DRDK_BUILD_AVALON_SUPPORT=ON",
            f"-DRDK_BUILD_PYTHON_WRAPPERS=ON",
            f"-DRDK_BUILD_YAEHMOP_SUPPORT=ON",
            f"-DRDK_INSTALL_INTREE=OFF",
            # Boost
            # f"-DBOOST_ROOT={boost_install_path}",
            f"-DBoost_NO_SYSTEM_PATHS=ON",
            f"-DBoost_DEBUG=OFF",
            # Does not work (this is fixed in future rdkit versions I believe)
            f"-DRDK_INSTALL_STATIC_LIBS=OFF" if sys.platform == "win32" else "",
            # ##### for windows
            # cairo
            f"-DRDK_BUILD_CAIRO_SUPPORT=ON",
            # f"-DCAIRO_INCLUDE_DIR={towin(vcpkg_include_path)}"
            # if sys.platform == "win32"
            # else "",
            # f"-DCAIRO_LIBRARY_DIR={towin(vcpkg_lib_path)}"
            # if sys.platform == "win32"
            # else "",
            # zlib
            # f"-DZLIB_ROOT={towin(vcpkg_install_path)}"
            # if sys.platform == "win32"
            # else "",
            # freetype
            # f"-DFREETYPE_INCLUDE_DIRS={towin(vcpkg_include_path)}"
            # if sys.platform == "win32"
            # else "",
            # f"-DFREETYPE_LIBRARY={towin(vcpkg_lib_path / 'freetype.lib')}"
            # if sys.platform == "win32"
            # else "",
            # eigen3
            # f"-DEIGEN3_INCLUDE_DIR={towin(vcpkg_include_path)}"
            # if sys.platform == "win32"
            # else "",
            # instruct to build x64
            "-Ax64" if sys.platform == "win32" else "",
            # Mac needs these flags to compile
            f"-DCMAKE_C_FLAGS=-Wno-implicit-function-declaration"
            if sys.platform == "darwin"
            else "",
            f"-DCMAKE_CXX_FLAGS=-Wno-implicit-function-declaration"
            if sys.platform == "darwin"
            else "",
            # build stuff
            f"-DCMAKE_INSTALL_PREFIX={rdkit_install_path}",
            f"-DCMAKE_BUILD_TYPE=Release",
            f"-GNinja" if sys.platform != "win32" else "",
        ]

        cmds = [
            f"cmake -S . -B build {' '.join(options)} ",
            f"cmake --build build --config Release",
            f"cmake --install build",
        ]
        [check_call(c.split()) for c in cmds]

        os.chdir(str(cwd))

        # ### copy files
        # rdkit_install_path = Path(self.build_temp).absolute() / "rdkit_install"
        py_name = "python" + ".".join(map(str, sys.version_info[:2]))
        rdkit_files = rdkit_install_path / "lib" / py_name / "site-packages" / "rdkit"

        if sys.platform == "win32":
            rdkit_files = rdkit_install_path / "Lib" / "site-packages" / "rdkit"

        rdpaths_file = rdkit_files / "RDPaths.py"
        rdkit_data_path = rdkit_install_path / "share" / "RDKit" / "Data"

        # copy rdkit files here
        wheel_path = Path(self.get_ext_fullpath(ext.name)).absolute()
        if wheel_path.exists():
            rmtree(str(wheel_path))

        # Modify RDPaths.py
        sed = "gsed" if platform == "darwin" else "sed"
        call(
            [
                sed,
                "-i",
                "/_share =/c\_share = os.path.dirname(__file__)",
                str(rdpaths_file),
            ]
        )

        # Normal files
        copytree(str(rdkit_files), str(wheel_path))
        # copy the Data directory to the wheel path
        copytree(str(rdkit_data_path), str(wheel_path / "Data"))

        # Copy so and dylib files to /usr/local/lib
        rdkit_root = rdkit_install_path / "lib"

        libs_rdkit_linux = Path(rdkit_root).glob("*.so*")
        libs_rdkit_macos = Path(rdkit_root).glob("*dylib")
        libs_rdkit = list(libs_rdkit_linux) + list(libs_rdkit_macos)

        if platform != "win32":
            [copy_file(i, "/usr/local/lib") for i in libs_rdkit]
        else:
            # windows is Lib or libs?
            rdkit_root = rdkit_install_path / "Lib"
            libs_rdkit_win = Path(rdkit_root).glob("*.dll")

            libs_vcpkg = list(
                (vcpkg_path / "installed" / "x64-windows" / "bin").glob("*.dll")
            )
            [copy_file(i, "C://libs") for i in libs_rdkit_win]
            [copy_file(i, "C://libs") for i in libs_vcpkg]


setup(
    name="rdkit-pypi",
    version=f"2021.9.4",
    description="A collection of chemoinformatics and machine-learning software written in C++ and Python",
    author="Christopher Kuenneth",
    author_email="chris@kuenneth.dev",
    url="https://github.com/kuelumbus/rdkit-pypi",
    project_urls={
        "RDKit": "http://rdkit.org/",
        "RDKit on Github": "https://github.com/rdkit/rdkit",
    },
    license="BSD-3-Clause",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19",
        "Pillow",
    ],
    ext_modules=[
        RDKit(
            "rdkit",
            # 1.73 does now compile on win for some reason
            boost_download_urls={
                "win": "https://boostorg.jfrog.io/artifactory/main/release/1.69.0/source/boost_1_69_0.tar.gz",
                "mac": "https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz",
                "linux": "https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz",
            },
            rdkit_tag="Release_2021_09_4",
        ),
    ],
    cmdclass=dict(build_ext=BuildRDKit),
)
