import os
import shlex
import sys
from distutils.file_util import copy_file
from pathlib import Path
from shutil import copytree, rmtree, ignore_patterns
from subprocess import call, check_call
import sysconfig
from textwrap import dedent

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext as build_ext_orig

# RDKit version to build (tag from github repository)
rdkit_tag = "Release_2024_09_1"

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


class RDKit(Extension):
    def __init__(self, name, **kwargs):
        super().__init__(name, sources=[])
        self.__dict__.update(kwargs)


class BuildRDKit(build_ext_orig):
    def run(self):
        for ext in self.extensions:
            self.build_rdkit(ext)
        super().run()

    def get_ext_filename(self, ext_name):
        ext_path = ext_name.split(".")
        return os.path.join(*ext_path)

    def conan_install(self, boost_version, conan_toolchain_path):
        """Run the Conan"""

        # This modified conanfile.py for boost does not link libpython*.so
        # When building a platform wheel, we don't want to link libpython*.so.
        mod_conan_path = "conan_boost_mod"

        # Export the modified boost version
        check_call(
            [
                "conan",
                "export",
                f"{mod_conan_path}/all/",
                "--version",
                f"{boost_version}",
                "--user",
                "chris",
                "--channel",
                "mod_boost",
            ]
        )

        # conan2 needs a default profile
        check_call(
            [
                "conan",
                "profile",
                "detect",
            ]
        )

        without_python_lib = "boost/*:without_python_lib=False"
        boost_version_string = f"boost/{boost_version}@chris/mod_boost"
        without_stacktrace = "False"

        if sys.platform != "win32":
            # if no windows builds, compile boost without python lib.a/.so/.dylib
            without_python_lib = "boost/*:without_python_lib=True"

        if "macosx_arm64" in os.environ["CIBW_BUILD"]:
            # does not work on macos arm64 for some reason
            without_stacktrace = "True"

        macos_libs = ""
        if "macosx" in os.environ["CIBW_BUILD"]:
            ## install these libraries to meet the development target
            macos_libs = """
pixman/0.43.4
cairo/1.18.0
libpng/1.6.43
fontconfig/2.15.0
freetype/2.13.2
"""

        conanfile = f"""\
            [requires]
            {boost_version_string}
            {macos_libs}

            [generators]
            CMakeDeps
            CMakeToolchain
            VirtualRunEnv

            [options]
            boost/*:shared=True
            boost/*:without_python=False
            {without_python_lib}
            boost/*:python_executable={sys.executable}
            boost/*:without_stacktrace={without_stacktrace}
        """
        # boost:debug_level=1

        Path("conanfile.txt").write_text(dedent(conanfile))

        # run conan install
        cmd = [
            "conan",
            "install",
            "conanfile.txt",
            # build all missing
            "--build=missing",
            "--output-folder",
            f"{conan_toolchain_path}",
        ]

        if sys.platform == "win32":
            cmd += ["-pr:b", "default"]

        # but force build b2 on linux
        if "linux" in sys.platform:
            cmd += ["--build=b2", "-pr:b", "default"]

        check_call(cmd)

    def build_rdkit(self, ext):
        """Build RDKit

        Steps:
        (1) Use Conan to install boost and other libraries
        (2) Build RDKit
        (3) Copy RDKit and additional files to the wheel path
        (4) Copy the libraries to system paths
        """

        cwd = Path().absolute()

        # Install boost and other libraries using Conan
        conan_toolchain_path = cwd / "conan"
        conan_toolchain_path.mkdir(parents=True, exist_ok=True)
        boost_version = "1.85.0"

        boost_lib_version = "_".join(boost_version.split(".")[:2])
        self.conan_install(boost_version, conan_toolchain_path)

        # Build RDkit
        # Define paths
        build_path = Path(self.build_temp).absolute()
        build_path.mkdir(parents=True, exist_ok=True)
        os.chdir(str(build_path))

        rdkit_install_path = build_path / "rdkit_install"
        rdkit_install_path.mkdir(parents=True, exist_ok=True)

        # Clone RDKit from git at rdkit_tag
        check_call(
            ["git", "clone", "-b", f"{ext.rdkit_tag}", "https://github.com/rdkit/rdkit"]
        )

        # Location of license file
        license_file = build_path / "rdkit" / "license.txt"

        # Start build process
        os.chdir(str("rdkit"))

       
        import fileinput

        def replace_all(file, search_exp, replace_exp):
            with fileinput.input(file, inplace=True) as f:
                for line in f:
                    if search_exp in line:
                        line = line.replace(search_exp, replace_exp)
                    print(line, end="")

        # Fix a bug in conan or rdkit: target name for numpy is boost::numpy{pyversion} with small 'b'
        # and not Boost::numpy{pyversion}
        # Line 312 in 2024_03_06 in CMakeLists.txt
        # NEW: target_link_libraries(rdkit_py_base INTERFACE "Boost::python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}" "Boost::numpy${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")
        replace_all(
            "CMakeLists.txt",
            'target_link_libraries(rdkit_py_base INTERFACE "Boost::python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}" "Boost::numpy${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")',
            'target_link_libraries(rdkit_py_base INTERFACE "boost::python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}" "boost::numpy${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")',
        )

        if "macosx" in os.environ["CIBW_BUILD"]:
            # Replace Cairo with cairo because conan uses lower case target names
            # only on MacOS cairo is installed using conan
            replace_all(
                "Code/GraphMol/MolDraw2D/CMakeLists.txt",
                'target_link_libraries(MolDraw2D PUBLIC Cairo::Cairo)',
                'target_link_libraries(MolDraw2D PUBLIC cairo::cairo)',
            )
            replace_all(
                "Code/GraphMol/MolDraw2D/CMakeLists.txt",
                'target_link_libraries(MolDraw2D_static PUBLIC Cairo::Cairo)',
                'target_link_libraries(MolDraw2D_static PUBLIC cairo::cairo)',
            )

        # introduced in 2024_09_01 for compiling pubchem shape.
        replace_all(
            "External/pubchem_shape/Wrap/CMakeLists.txt",
            'find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)',
            'find_package(Python3 COMPONENTS Interpreter Development NumPy)',
        )
        


        print("---- Conf vars", file=sys.stderr)
        print(sysconfig.get_paths(), file=sys.stderr)
        print(sysconfig.get_config_vars(), file=sys.stderr)
        print("---- Conf vars", file=sys.stderr)
        


        print("---- Conf vars", file=sys.stderr)
        print(sysconfig.get_paths(), file=sys.stderr)
        print(sysconfig.get_config_vars(), file=sys.stderr)
        print("---- Conf vars", file=sys.stderr)
        

        # Define CMake options
        options = [
            # f"-DCMAKE_FIND_DEBUG_MODE=ON", # Enable debug mode
            f"-DCMAKE_TOOLCHAIN_FILE={conan_toolchain_path / 'conan_toolchain.cmake'}",
            # For the toolchain file this needs to be set
            f"-DCMAKE_POLICY_DEFAULT_CMP0091=NEW",
            # Boost_VERSION_STRING is set but Boost_LIB_VERSION is not set by conan.
            # Boost_LIB_VERSION is required by RDKit => Set manually
            f"-DBoost_LIB_VERSION={boost_lib_version}",
            # Select correct python 3 version
            f"-DPython3_ROOT_DIR={Path(sys.prefix)}",
            # RDKit build flags
            "-DRDK_BUILD_INCHI_SUPPORT=ON",
            "-DRDK_BUILD_AVALON_SUPPORT=ON",
            "-DRDK_BUILD_PYTHON_WRAPPERS=ON",
            "-DRDK_BUILD_YAEHMOP_SUPPORT=ON",
            "-DRDK_BUILD_XYZ2MOL_SUPPORT=ON",
            "-DRDK_INSTALL_INTREE=OFF",
            "-DRDK_BUILD_CAIRO_SUPPORT=ON",
            "-DRDK_BUILD_FREESASA_SUPPORT=ON",
            # Disable system libs for finding boost
            "-DBoost_NO_SYSTEM_PATHS=ON",
            # build stuff
            f"-DCMAKE_INSTALL_PREFIX={rdkit_install_path}",
            "-DCMAKE_BUILD_TYPE=Release",
            # Speed up builds
            "-DRDK_BUILD_CPP_TESTS=OFF",
            # Fix InChi download
            "-DINCHI_URL=https://rdkit.org/downloads/INCHI-1-SRC.zip",
        ]

        # Modifications for Windows
        vcpkg_path = cwd
        vcpkg_inc = vcpkg_path / "vcpkg_installed" / "x64-windows" / "include"
        vcpkg_lib = vcpkg_path / "vcpkg_installed" / "x64-windows" / "lib"

        if sys.platform == "win32":
            def to_win_path(pt: Path):
                return str(pt).replace("\\", "/")
        
            options += [
                # DRDK_INSTALL_STATIC_LIBS should be fixed in newer RDKit builds. Remove?
                "-DRDK_INSTALL_STATIC_LIBS=OFF",
                "-DRDK_INSTALL_DLLS_MSVC=ON",
            ]

            # Link cairo and freetype
            options += [
                f"-DCAIRO_INCLUDE_DIR={to_win_path(vcpkg_inc)}",
                f"-DCAIRO_LIBRARY_DIR={to_win_path(vcpkg_lib)}",
                f"-DFREETYPE_INCLUDE_DIRS={to_win_path(vcpkg_inc)}",
                f"-DFREETYPE_LIBRARY={to_win_path(vcpkg_lib / 'freetype.lib')}",
            ]

        # Modifications for MacOS all
        if sys.platform == "darwin":
            options += [
                "-DCMAKE_C_FLAGS=-Wno-implicit-function-declaration",
                # CATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS because MacOS does not fully support C++17.
                '-DCMAKE_CXX_FLAGS="-Wno-implicit-function-declaration -DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS"',
            ]

        # Modification for MacOS x86_64
        if "macosx_x86_64" in os.environ["CIBW_BUILD"]:
            options += [
                # macOS < 10.13 has a incomplete C++17 implementation
                # See https://github.com/kuelumbus/rdkit-pypi/pull/85 for a discussion
                f"-DCMAKE_OSX_DEPLOYMENT_TARGET={os.environ.get('MACOSX_DEPLOYMENT_TARGET', '10.13')}",
            ]

        # Modifications for MacOS arm64 (M1 hardware)
        if "macosx_arm64" in os.environ["CIBW_BUILD"]:
            options += [
                "-DRDK_OPTIMIZE_POPCNT=OFF",
                # Otherwise, cmake tries to link the system freetype
                # "-DFREETYPE_LIBRARY=/opt/homebrew/lib/libfreetype.dylib",
                # "-DFREETYPE_INCLUDE_DIRS=/opt/homebrew/include",
                # Arm64 build start with development target 11.0
                f"-DCMAKE_OSX_DEPLOYMENT_TARGET={os.environ.get('MACOSX_DEPLOYMENT_TARGET', '11.0')}",
                f"-DCMAKE_OSX_ARCHITECTURES=arm64",
                f"-DCMAKE_VERBOSE_MAKEFILE=ON" # Increase verbosity
            ]

        if "linux" in sys.platform:
            # Use ninja for linux builds
            cmds = [
                f"cmake -S . -B build -G Ninja --debug-find-pkg=Python3 {' '.join(options)} ",
                "cmake --build build --config Release",
                "cmake --install build",
            ]
        elif sys.platform == "win32":
            cmds = [
                f"cmake -S . -B build --debug-find-pkg=Python3 {' '.join(options)} ",
                "cmake --build build --config Release -v",
                "cmake --install build",
            ]
        else:
            cmds = [
                f"cmake -S . -B build --debug-find-pkg=Python3 {' '.join(options)} ",
                "cmake --build build --config Release",
                "cmake --install build",
            ]

        # Define the rdkit_files path
        py_name = "python" + ".".join(map(str, sys.version_info[:2]))

        path_site_packages = rdkit_install_path / "lib" / py_name / "site-packages"
        if sys.platform == "win32":
            path_site_packages = rdkit_install_path / "Lib" / "site-packages"

        print("!!! --- CMAKE build command and variables for RDKit", file=sys.stderr)
        print(cmds, file=sys.stderr)
        variables = {}
        print(variables, file=sys.stderr)

        # Run CMake and install RDKit
        [
            check_call(
                shlex.split(c, posix="win32" not in sys.platform),
                env=dict(os.environ, **variables),
            )
            for c in cmds
        ]

        # --- Copy libs to system path
        # While repairing the wheels, the built libs need to be copied to the platform wheels
        # Also, the libs needs to be accessible for building the stubs
        rdkit_lib_path = rdkit_install_path / "lib"
        boost_lib_path = conan_toolchain_path / "boost" / "lib"
        boost_lib_path_bin_windows_only = conan_toolchain_path / "boost" / "bin"

        cmds = []
        if "linux" in sys.platform:
            # Libs end with .so
            to_path = Path("/usr/local/lib")
            [copy_file(i, str(to_path)) for i in rdkit_lib_path.rglob("*.so*")]
            [copy_file(i, str(to_path)) for i in boost_lib_path.rglob("*.so*")]
            cmds.append("ldconfig")

        elif "win32" in sys.platform:
            # Libs end with .dll
            # windows paths are case insensitive
            # C://libs is specified as search dir for repairing the wheel
            dll_paths = ["C://Windows//SYSWOW64", "C://Windows//System32", "C://libs"]
            for pt in dll_paths:
                to_path = Path(pt)
                to_path.mkdir(parents=True, exist_ok=True)
                [copy_file(i, str(to_path)) for i in rdkit_lib_path.rglob("*.dll")]
                [copy_file(i, str(to_path)) for i in rdkit_lib_path.rglob("*.pyd")]
                [copy_file(i, str(to_path)) for i in rdkit_lib_path.rglob("*.lib")]

                [copy_file(i, str(to_path)) for i in boost_lib_path.rglob("*.lib")]
                [
                    copy_file(i, str(to_path))
                    for i in boost_lib_path_bin_windows_only.rglob("*.dll")
                ]

                variables["PATH"] = os.environ["PATH"] + os.pathsep + str(to_path)

            # VCPKG libs
            variables["PATH"] = os.environ["PATH"] + os.pathsep + str(vcpkg_lib)

        elif "darwin" in sys.platform:
            # Github actions
            to_path = Path("/Users/runner/work/lib")
            if "CIRRUS_CI" in os.environ:
                # on cirrus CI
                to_path = Path("/Users/admin/lib")

            # Make sure path exists?
            to_path.mkdir(parents=True, exist_ok=True)

            # Add path to DYLD_LIBRARY_PATH for generating stubs
            variables["DYLD_LIBRARY_PATH"] = str(to_path)

            # copy all boost and rdkit libs to one path
            [copy_file(i, str(to_path)) for i in rdkit_lib_path.rglob("*dylib")]
            [copy_file(i, str(to_path)) for i in boost_lib_path.rglob("*dylib")]

        # Build the RDKit stubs
        cmds += [
            f"cmake --build build --config Release --target stubs -v",
        ]

        # rdkit-stubs require the site-package path to be in sys.path / PYTHONPATH
        variables["PYTHONPATH"] = (
            os.environ.get("PYTHONPATH", "") + os.pathsep + str(path_site_packages)
        )

        print(
            "!!! --- CMAKE build command and variables for building stubs",
            file=sys.stderr,
        )
        print(cmds, file=sys.stderr)
        print(variables, file=sys.stderr)

        [
            check_call(
                shlex.split(c, posix="win32" not in sys.platform),
                env=dict(os.environ, **variables),
            )
            for c in cmds
        ]

        # Print the stubs error file to rdkit-stubs/gen_rdkit_stubs.err
        stubs_error_file = (
            build_path / "rdkit" / "build" / "rdkit-stubs" / "gen_rdkit_stubs.err"
        )
        with open(stubs_error_file, "r") as fin:
            print(fin.read(), file=sys.stderr)

        os.chdir(str(cwd))

        # Copy RDKit and additional files to the wheel path
        # Modify RDPaths.py
        sed = "gsed" if sys.platform == "darwin" else "sed"
        call(
            [
                sed,
                "-i",
                "/_share =/c\_share = os.path.dirname(__file__)",  # noqa: W605
                f"{path_site_packages / 'rdkit'/ 'RDPaths.py'}",
            ]
        )

        # RDKit stubs directory
        dir_rdkit_stubs = path_site_packages / "rdkit-stubs"

        # Data directory
        rdkit_data_path = rdkit_install_path / "share" / "RDKit" / "Data"

        # Contrib directory
        rdkit_contrib_path = rdkit_install_path / "share" / "RDKit" / "Contrib"

        # Setuptools searches at this path for files to include
        wheel_path = Path(self.get_ext_fullpath(ext.name)).absolute().parent
        wheel_path.mkdir(exist_ok=True)

        # Copy RDMKit files to .../rdkit directory
        def _logpath(path, names):
            ignore_patterns
            print(f"In directory {path} copy files: {names}", file=sys.stderr)
            return ignore_patterns("*.pyc")(path, names)

        # Copy the RDKit stubs files to the rdkit-stubs wheels path
        copytree(dir_rdkit_stubs, wheel_path / "rdkit-stubs", ignore=_logpath)
        # Copy the Python files
        copytree(path_site_packages / "rdkit", wheel_path / "rdkit", ignore=_logpath)
        # Copy the data directory
        copytree(rdkit_data_path, wheel_path / "rdkit" / "Data", ignore=_logpath)
        # Copy the contrib directory
        copytree(rdkit_contrib_path, wheel_path / "rdkit" / "Contrib", ignore=_logpath)

        # Delete some large files from the Contrib folder
        # that are not necessary for running RDKit
        # See https://github.com/rdkit/rdkit/issues/5601
        _dir = wheel_path / "rdkit" / "Contrib" / "NIBRSubstructureFilters"
        rmtree(str(_dir / "examples"))
        (_dir / "FilterSet_NIBR2019_wPubChemExamples.html").unlink()
        (_dir / "filterExamples.png").unlink()

        _dir = wheel_path / "rdkit" / "Contrib" / "CalcLigRMSD"
        rmtree(str(_dir / "data"))
        rmtree(str(_dir / "figures"))
        (_dir / "Examples_CalcLigRMSD.ipynb").unlink()

        # Copy the license
        copy_file(str(license_file), str(wheel_path / "rdkit"))


setup(
    name="rdkit",
    version=rdkit_tag.replace("Release_", "").replace("_", "."),
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
    install_requires=[
        "numpy",
        "Pillow",
    ],
    ext_modules=[
        RDKit("rdkit", rdkit_tag=rdkit_tag),
    ],
    cmdclass=dict(build_ext=BuildRDKit),
)
