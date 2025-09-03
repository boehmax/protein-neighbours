#!/usr/bin/env python3
"""Setup script for protein-neighbours Python package."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="protein-neighbours",
    version="0.2.0",
    author="Maximilian Böhm",
    author_email="maximilian.bohm@kemi.uu.se",
    description="A Python package for analyzing the genomic environment of proteins",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/boehmax/protein-neighbours",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "biopython>=1.79",
        "pyyaml>=5.4.0",
        "scipy>=1.7.0",
        "jinja2>=3.0.0",
        "click>=8.0.0",
        "tqdm>=4.60.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.12",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.900",
        ],
        "eggnog": [
            "eggnog-mapper>=2.1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "protein-neighbours=protein_neighbours.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "protein_neighbours": [
            "config/*.yaml",
            "templates/*.html",
        ],
    },
)