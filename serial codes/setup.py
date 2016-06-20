from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize("c_functions.pyx"),
    include_dirs=[numpy.get_include()]
)

# setup(
#     ext_modules = cythonize("increment_pixel_radiance.pyx"),
#     include_dirs=[numpy.get_include()]
# )

# setup(
#     ext_modules = cythonize("hello_world.pyx"),
#     include_dirs = [numpy.get_include()]
# )