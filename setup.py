"""
For pip install installation.
"""

import os
import fnmatch

import setuptools

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def ext_modules(config, _dir):
    """Fetch f90 files in src and automatically create an extension"""
    pattern = "*.f90"
    if os.path.isdir(_dir):
        for root, dirs, files in os.walk(_dir):
            match = fnmatch.filter(files, pattern)
            for name in match:
                f90_file = os.path.join(root, name)
                ext_name = os.path.splitext(f90_file)[0].replace("/", ".")
                config.add_extension(ext_name,
                                     [f90_file],
                                     f2py_options=['--quiet']
                                    )

def run_setup():
    """Setup"""
    config = Configuration(None, '', '')
    ext_modules(config, "curp")
    config.add_data_files(os.path.join("curp", "LICENSE-short.txt"))
    setup(
        name="curp",
        version="1.2dev1",
        author="Yamato's Lab",
        author_email="yamato@nagoya-u.jp",
        description="Inter-residue Current calculation in Proteins from MD \
            trajectory",
        url=("https://gitlab.com/yamato97/current-calculations-for-proteins"),
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "Operating System :: MacOS",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Fortran",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Physics"
            ],

        install_requires=["numpy>=1.11.2",
                          "nose",
                          "mpi4py>=2.0",
                          "benchmarker",
                          "setproctitle",
                          "epydoc",
                          "pygraphviz",
                          "netCDF4>=1.2.4"],
        #packages=["curp",
        #          "curp.current",
        #          "curp.dynamics",
        #          "curp.forcefield",
        #          "curp.parser",
        #          "curp.twobody",
        #          "curp.table",
        #          "curp.volume",
        #          ],

        packages=setuptools.find_packages(),

        python_requires=">=2.7, <3.0",
        entry_points={
            "console_scripts": [
                "curp = curp.curp:main",
                "cal_tc = curp.script.cal_tc:main",
                "conv_trj = curp.script.conv_trj:main",
                "graph_een = curp.script.graph_een:main"
                ]
            },
        **config.todict()
        )

if __name__ == '__main__':
    run_setup()
