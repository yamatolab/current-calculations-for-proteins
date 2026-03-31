"""
For pip install installation.
"""

import os
import fnmatch

import setuptools
import sys
import numpy
import pybind11
from pathlib import Path
import subprocess

sys.path.append("./curp/")
from _version import __version__

def ext_modules(config, _dir):
    """Fetch f90 files in src and automatically create an extension"""
    pattern = "*.f90"
    pattern_cpp = "*.cpp"
    
    try:
        conda_prefix = Path(os.environ["CONDA_PREFIX"])
    except KeyError:
        conda_prefix = Path("/usr")
    MPI_DIR = os.environ.get("MPI_DIR", conda_prefix)
    NETCDF_DIR = os.environ.get("NETCDF_DIR", conda_prefix)
    GRAPHVIZ_DIR = os.environ.get("GRAPHVIZ_DIR", conda_prefix)
    EIGEN_DIR = os.environ.get("EIGEN3_INCLUDE_DIR", conda_prefix)

    # Typical library/include paths:
    mpi_inc = os.path.join(MPI_DIR, "include")
    mpi_lib = os.path.join(MPI_DIR, "lib")
    netcdf_inc = os.path.join(NETCDF_DIR, "include")
    netcdf_lib = os.path.join(NETCDF_DIR, "lib")
    graphviz_inc = os.path.join(GRAPHVIZ_DIR, "include")
    graphviz_lib = os.path.join(GRAPHVIZ_DIR, "lib/x86_64-linux-gnu/graphviz")
    eigen_inc = os.path.join(EIGEN_DIR, "include/eigen3")
    if not os.path.isdir(eigen_inc):
        raise RuntimeError(
            f"Eigen3 include directory not found: {eigen_inc}. ")
    pybind11_inc = subprocess.check_output(
        ["python3", "-m", "pybind11", "--includes"], universal_newlines=True
    ).strip()
    
    extra_f90_compile_args = [
        "-O1", 
        "-fopenmp",
        f"-I{mpi_inc}",
        f"-I{netcdf_inc}",
        f"-I{graphviz_inc}",
    ]
    extra_link_args = [
        "-lgomp",
        f"-L{mpi_lib}", "-lmpi",
        f"-L{netcdf_lib}", "-lnetcdf",
        f"-L{graphviz_lib}", "-lgraphviz",
    ]
    
    if os.path.isdir(_dir):
        # add fortran scripts
        for root, dirs, files in os.walk(_dir):
            match = fnmatch.filter(files, pattern)
            for name in match:
                f90_file = os.path.join(root, name)
                ext_name = os.path.splitext(f90_file)[0].replace("/", ".")
                config.add_extension(ext_name,
                                     [f90_file],
                                     f2py_options=["--quiet"],
                                     extra_f90_compile_args=["-O1", "-fopenmp"],
                                     extra_link_args=["-lgomp"],
                                    )
        # add cpp scripts
        for root, dirs, files in os.walk(_dir):
            match = fnmatch.filter(files, pattern_cpp)
            for name in match:
                cpp_file = os.path.join(root, name)
                ext_name = os.path.splitext(cpp_file)[0].replace("/", ".")
                include_dirs = [numpy.get_include(), pybind11.get_include(), eigen_inc]
                config.add_extension(ext_name,
                                     [cpp_file],
                                     include_dirs=include_dirs,
                                     libraries=["mpi", "netcdf"],
                                     extra_compile_args=["-O3", "-Wall","-std=c++14", "-fopenmp",
                                                         "-fPIC", pybind11_inc],
                                     extra_link_args=["-lgomp"],
                                    )

with open('README.rst', 'r', encoding='utf-8') as summary:
    LONG_DESCRIPTION = summary.read()

def run_setup():
    """Setup"""
    from numpy.distutils.core import setup
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration(None, '', '')
    ext_modules(config, "curp")
    config.add_data_files(os.path.join("curp", "LICENSE-short.txt"))
    setup(
        name="Curp",
        version=__version__,
        author="Yamato's Lab",
        author_email="yamato@nagoya-u.jp",
        description="Inter-residue Current calculation in Proteins from MD \
            trajectory",
        long_description=LONG_DESCRIPTION,
        long_description_content_type='text/x-rst',
        python_requires=">3.5,<3.7",
        url=("https://github.com/yamatolab/current-calculations-for-proteins"),
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Operating System :: MacOS",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Fortran",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Chemistry"
            ],

        install_requires=["numpy>=1.11.2,<1.17",
                          "nose==1.3.7",
                          "mpi4py>=1.2",
                          "pygraphviz>1.2,<1.6",
                          "netcdf4>=1.4.2,<1.7",
                          "pybind11>2.11,<2.13"],
        
        setup_requires = ["numpy>1.11.2,<1.17",
                          "pybind11>2.11,<2.13"],
        
        extras_require={
            "dev": ["benchmarker>=4.0,<5",]
        },

        packages=setuptools.find_packages(),

        package_data={'curp':['volume/random20.pdb.gz']},
        
        entry_points={
            "console_scripts": [
                "curp = curp.console:main",
                ]
            },
        **config.todict()
        )

if __name__ == '__main__':
    run_setup()
