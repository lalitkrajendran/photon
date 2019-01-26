Instructions for running the software.
- Compile parallel_ray_tracing.cu under cuda_codes (see dependencies below).
- Use create_simulation_parameters.py to create a set of parameters in a python dictionary and save to a .mat file.
- Use batch run simulation to read parameter file (.mat) and perform the ray tracing.

Dependencies for Python:
- numpy
- scipy
- matplotlib
- libtiff
- ctypes

Dependencies for CUDA:
- bz2 (https://www.sourceware.org/bzip2/)
- png (http://www.libpng.org/pub/png/libpng.html)
- teem (http://teem.sourceforge.net)
- cubic_interpolation_cuda (https://github.com/DannyRuijters/CubicInterpolationCUDA)

