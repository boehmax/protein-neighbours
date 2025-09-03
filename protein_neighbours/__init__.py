"""
Protein Neighbours - A Python package for analyzing protein genomic environments.

This package provides comprehensive tools for analyzing the genomic neighborhood 
of proteins, including neighbor identification, annotation, and visualization.
"""

__version__ = "0.2.0"
__author__ = "Maximilian Böhm"
__email__ = "maximilian.bohm@kemi.uu.se"

from .utils import load_config, setup_logging
from .io import read_protein_assembly_data, read_clades
from .analysis import analyze_neighbor_types, generate_summary_statistics

__all__ = [
    "load_config",
    "setup_logging", 
    "read_protein_assembly_data",
    "read_clades",
    "analyze_neighbor_types",
    "generate_summary_statistics",
]