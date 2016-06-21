import glob
import os

from setuptools import setup, find_packages

MPINT_DIR = os.path.dirname(os.path.abspath(__file__))

setup(
    name="mdtools",
    version="0.0.0",
    install_requires=["pymatgen>=3.6.1"],
    extras_require={"babel": ["openbabel", "pybel"]},
    author="Kiran Mathew, Brandon Wood",
    author_email="kmathew@lbl.gov",
    description=("Molecular Dynamics Helper Tools"),
    license="MIT",
    url="https://github.com/matk86/MDTools",
    packages=find_packages(),
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    scripts=glob.glob(os.path.join(MPINT_DIR, "scripts", "*"))
)
