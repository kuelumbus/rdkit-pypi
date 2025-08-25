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
        deps.generate()
        
        # Generate CMake toolchain
        tc = CMakeToolchain(self)
        tc.generate()
        
        # Generate virtual run environment
        env = VirtualRunEnv(self)
        env.generate()
