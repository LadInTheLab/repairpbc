#!/usr/bin/env python3
"""
Setup script for RepairPBC tool.
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="repairpbc",
    version="1.0.0",
    author="Andrew Brodrick",
    description="A tool to repair periodic boundary conditions in GROMACS trajectory files",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/andrewbrodrick/repairpbc",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.6",
    install_requires=read_requirements(),
    entry_points={
        "console_scripts": [
            "repairpbc=repair_pbc:main",
        ],
    },
    keywords="gromacs, molecular dynamics, trajectory, pbc, periodic boundary conditions",
    project_urls={
        "Bug Reports": "https://github.com/LadInTheLab/repairpbc/issues",
        "Source": "https://github.com/LadInTheLab/repairpbc",
    },
) 