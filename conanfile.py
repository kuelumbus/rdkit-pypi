import os
import sys
from conan import ConanFile
from conan.tools.cmake import CMakeDeps, CMakeToolchain
from conan.tools.env import VirtualRunEnv
from conan.tools.files import copy
from pathlib import Path



class RDKitConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    
    def configure(self):
        # Configure boost options
        self.options["boost/*"].shared = True
        self.options["boost/*"].without_python = False

        # We always need a posix path with forward slashes
        # Because the workflows run on Windows runners with the Git Bash shell,
        # as_posix returns "/" paths
        self.options["boost/*"].python_executable =  Path(sys.executable).as_posix()

        # Platform-specific configurations
        if self.settings.os == "Macos" and self.settings.arch == "armv8":
            # stacktrace does not work on macOS arm64 for some reason
            self.options["boost/*"].without_stacktrace = True
        else:
            self.options["boost/*"].without_stacktrace = False
            
        # Configure Python library linking for wheel building
        if self.settings.os == "Windows":
            self.options["boost/*"].without_python_lib = False
        else:
            self.options["boost/*"].without_python_lib = True

    def requirements(self):
        # Main boost requirement - use modified version
        self.requires("boost/1.85.0@chris/mod_boost")
        # self.requires("boost/1.85.0")
        self.requires("expat/2.7.5")
        
        # Platform-specific requirements
        if self.settings.os == "Macos" and os.environ.get("CIBW_BUILD", "").startswith("cp"):
            # macOS libraries to meet development target
            self.requires("pixman/0.43.4")
            self.requires("cairo/1.18.0") 
            self.requires("libpng/1.6.43")
            self.requires("fontconfig/2.15.0")
            self.requires("freetype/2.13.2")

    def build_requirements(self):
        pass

    def generate(self):
        # Generate CMake dependencies
        deps = CMakeDeps(self)

        # Force Conan to name the generated expat target 'EXPAT::EXPAT' instead of 'expat::expat'
        # introduced in 2026.03.4 because of ChemDraw parser requires this
        deps.set_property("expat", "cmake_target_name", "EXPAT::EXPAT")

        # Fix a bug in conan or rdkit: the boost python/numpy component targets are generated
        # as boost::python{X}{Y} / boost::numpy{X}{Y} (lowercase 'b'), but RDKit's CMakeLists.txt
        # expects Boost::python{X}{Y} / Boost::numpy{X}{Y}. Force the target names to match
        # instead of patching RDKit's CMakeLists.txt.
        py_major, py_minor = sys.version_info.major, sys.version_info.minor
        deps.set_property(
            f"boost::python{py_major}{py_minor}", "cmake_target_name", f"Boost::python{py_major}{py_minor}"
        )
        deps.set_property(
            f"boost::numpy{py_major}{py_minor}", "cmake_target_name", f"Boost::numpy{py_major}{py_minor}"
        )

        # Force Conan to name the generated cairo target 'Cairo::Cairo' instead of 'cairo::cairo'
        # so it matches what RDKit's MolDraw2D CMakeLists.txt expects, instead of patching it.
        deps.set_property("cairo", "cmake_target_name", "Cairo::Cairo")

        deps.generate()
        
        # Generate CMake toolchain
        tc = CMakeToolchain(self)

        # The vendored expatpp links EXPAT::EXPAT only PRIVATEly, so expat's include
        # dir does not propagate to consumers of expatpp.h (which #includes <expat.h>).
        # On Linux/macOS this is masked by expat.h living in default system include
        # paths; on MSVC it breaks the ChemDraw build (C1083). Inject the include dir
        # globally so headers always match the Conan expat we link against.
        expat_inc = self.dependencies["expat"].cpp_info.includedirs[0].replace("\\", "/")
        tc.extra_cxxflags.append("-I" + expat_inc)
        tc.extra_cflags.append("-I" + expat_inc)

        tc.generate()
        
        # Generate virtual run environment
        env = VirtualRunEnv(self)
        env.generate()
