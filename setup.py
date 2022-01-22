from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
from sysconfig import get_paths
import os
from subprocess import check_call, call, run, PIPE
import sys
from distutils.file_util import copy_file
from shutil import copytree, rmtree
from pathlib import Path
from textwrap import dedent


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
        boost_version = "1.69.0"

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
        # cairo/1.16.0
        win = """eigen/3.4.0
            freetype/2.11.1
        """
        if sys.platform != "win32":
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

        check_call(cmd)

    def build_rdkit(self, ext):
        """Build RDKit
        (1) Use Conan to install boost and other libs (for windows only)
        (2) Build RDKit
        (3) Copy the shared library to correct paths before start building the wheel
        """

        cwd = Path().absolute()

        # Install boost via conan
        conan_toolchain_path = cwd / "conan"
        conan_toolchain_path.mkdir(parents=True, exist_ok=True)
        self.conan_install(conan_toolchain_path)

        # Build path for everything
        build_path = Path(self.build_temp).absolute()
        build_path.mkdir(parents=True, exist_ok=True)
        os.chdir(str(build_path))

        # RDKit install path
        rdkit_install_path = build_path / "rdkit_install"
        rdkit_install_path.mkdir(parents=True, exist_ok=True)

        # Clone RDKit from git at rdkit_tag
        cmds = [f"git clone https://github.com/rdkit/rdkit"]
        [check_call(c.split()) for c in cmds]

        # For windows
        if sys.platform == "win32":
            check_call(
                [
                    "wget",
                    "https://github.com/preshing/cairo-windows/releases/download/1.15.12/cairo-windows-1.15.12.zip",
                    "--no-check-certificate",
                    "-q",
                ]
            )
            check_call(["tar", "-xzf", "cairo-windows-1.15.12.zip"])
            cairo_lib = build_path / "cairo-windows-1.15.12" / "lib" / "x64"
            cairo_inc = build_path / "cairo-windows-1.15.12" / "include"

        # Build path of rdkit
        os.chdir(str("rdkit"))

        # vcpkg_path = Path("C:/vcpkg")
        # vcpkg_include_path = vcpkg_path / "installed" / "x64-windows" / "include"
        # vcpkg_lib_path = vcpkg_path / "installed" / "x64-windows" / "lib"

        def path_to_win(pt: Path):
            return str(pt).replace("\\", "/")

        # CMake options
        options = [
            f"-DCMAKE_TOOLCHAIN_FILE={conan_toolchain_path / 'conan_paths.cmake'}",
            # Select correct python interpreter
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f'-DPYTHON_INCLUDE_DIR={get_paths()["include"]}',
            # RDKit build flags
            f"-DRDK_BUILD_INCHI_SUPPORT=ON",
            f"-DRDK_BUILD_AVALON_SUPPORT=ON",
            f"-DRDK_BUILD_PYTHON_WRAPPERS=ON",
            f"-DRDK_BUILD_YAEHMOP_SUPPORT=ON",
            f"-DRDK_INSTALL_INTREE=OFF",
            f"-DRDK_BUILD_CAIRO_SUPPORT=ON",
            # Disable system libs for finding boost
            f"-DBoost_NO_SYSTEM_PATHS=ON",
            # build stuff
            f"-DCMAKE_INSTALL_PREFIX={rdkit_install_path}",
            f"-DCMAKE_BUILD_TYPE=Release",
            f"-GNinja" if sys.platform != "win32" else "",
        ]

        if sys.platform == "win32":
            # DRDK_INSTALL_STATIC_LIBS should be fixed in newer RDKit builds
            options += ["-Ax64", "-DRDK_INSTALL_STATIC_LIBS=OFF"]

            # link cairo on windows
            options += [
                f"-DCAIRO_INCLUDE_DIR={path_to_win(cairo_inc)}",
                f"-DCAIRO_LIBRARY_DIR={path_to_win(cairo_lib)}",
            ]

        if sys.platform == "darwin":
            options += [
                "-DCMAKE_C_FLAGS=-Wno-implicit-function-declaration",
                "-DCMAKE_CXX_FLAGS=-Wno-implicit-function-declaration",
            ]

        cmds = [
            f"cmake -S . -B build {' '.join(options)}",
            f"cmake --build build",
            f"cmake --install build",
        ]
        [check_call(c.split()) for c in cmds]

        os.chdir(str(cwd))

        # copy libs to make them finable by wheels creators
        py_name = "python" + ".".join(map(str, sys.version_info[:2]))
        rdkit_files = rdkit_install_path / "lib" / py_name / "site-packages" / "rdkit"

        if sys.platform == "win32":
            rdkit_files = rdkit_install_path / "Lib" / "site-packages" / "rdkit"

        # Modify RDPaths.py
        sed = "gsed" if sys.platform == "darwin" else "sed"
        call(
            [
                sed,
                "-i",
                "/_share =/c\_share = os.path.dirname(__file__)",
                f"{rdkit_files / 'RDPaths.py'}",
            ]
        )

        # Data dir needs to be present in wheel
        rdkit_data_path = rdkit_install_path / "share" / "RDKit" / "Data"

        # copy rdkit files here, make sure it's empty
        wheel_path = Path(self.get_ext_fullpath(ext.name)).absolute()
        if wheel_path.exists():
            rmtree(str(wheel_path))

        # Copy python files
        copytree(str(rdkit_files), str(wheel_path))
        # Copy the data directory
        copytree(str(rdkit_data_path), str(wheel_path / "Data"))

        rdkit_root = rdkit_install_path / "lib"

        if "linux" in sys.platform:
            to_path = Path("/usr/local/lib")
            [copy_file(i, str(to_path)) for i in rdkit_root.glob("*.so*")]
        elif "win32" in sys.platform:
            # windows is Lib or libs?
            to_path = Path("C://libs")
            to_path.mkdir(parents=True, exist_ok=True)
            [copy_file(i, str(to_path)) for i in rdkit_root.glob("*.dll")]

        elif "darwin" in sys.platform:
            to_path = Path("/usr/local/lib")
            [copy_file(i, str(to_path)) for i in rdkit_root.glob("*dylib")]


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
