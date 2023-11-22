from conan import ConanFile
from conan.tools.cmake import CMake

class Conan(ConanFile):
    name = "StochasticPackageQuery"
    version = "1.0.0"

    # Define your package's settings and options if necessary
    settings = "os", "compiler", "build_type", "arch"

    # Define your package requirements
    def requirements(self):
        self.requires("boost/1.83.0")
        self.requires("highs/1.6.0")
        self.requires("libpq/15.4")
        self.requires("fmt/10.1.1")
        self.requires("pcg-cpp/cci.20220409")
        self.requires("timsort/2.1.0")
        self.requires("gmp/6.3.0")

    # Define your generators
    generators = "CMakeDeps", "CMakeToolchain"

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
