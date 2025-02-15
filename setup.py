"""
For pip install installation.
"""

import os
import fnmatch

import setuptools
import sys

sys.path.append("./curp/")
from _version import __version__

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
                                     f2py_options=["--quiet"],
                                     extra_f90_compile_args=["-O1", "-fopenmp"],
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
                          "netcdf4>=1.4.2,<1.7"],
        
        setup_requires = ["numpy>1.11.2,<1.17"],
        
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
