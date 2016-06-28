from setuptools import setup



setup(
    name="camera_simulation_package",
    version="1.0",
    author="Lalit Rajendran",
    author_email="lrajendr@purdue.edu",
    description=("A set of ray tracing codes that can render synthetic PIV/BOS 
      images for user-defined optical arrangements"),
    url="https://github.com/lrajendr/camera_simulation.git",
    packages=["scipy","numpy","matplotlib","skimage","sklearn",
      "ctypes","nrrd","progressbar"]
    )

