import setuptools
from numpy.distutils.core import setup, Extension
import sys
import os
import fnmatch

ext = Extension("lib_current",
                sources=["src/current/lib_current.f90"])

def ext_modules():
    """Fetch f90 files in src and create automatically an extension"""
    curp_dir = sys.path[0]
    ext_modules = []
    os.chdir(curp_dir)
    _dir = "src"
    pattern = "*.f90"
    if os.path.isdir(_dir):
        for root, dirs, files in os.walk(_dir):
            fs = fnmatch.filter(files, pattern)
            for name in fs:
                f90_file = os.path.join(root, name)
                ext_name = os.path.splitext(f90_file)[0].replace("/",".")
                ext_modules.append(Extension(ext_name,
                                             [f90_file],
                                             # f2py_options*['--quiet']
                                             ))
    return ext_modules



def run_setup():
    ext_modules()
    setup(
        name="CURP",
        version="1.2",
        author="Yamato's Lab",
        author_email="yamato@nagoya-u.jp",
        description="Current calculations between protein residues from MD trajectory",
        url=("https://gitlab.com/yamato97/current-calculations-for-proteins"),
        classifiers=[
            "Development Status :: ",
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
        install_requires=['numpy'],
        ext_modules=ext_modules()
        )

if __name__ == '__main__':
    run_setup()
