"""
Utility functions for protein neighbours analysis.

This module provides configuration loading, logging setup, and other utility functions.
"""

import os
import logging
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime


def load_config(config_file: str = "config/config.yaml", 
               override_params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_file: Path to the configuration YAML file
        override_params: Dictionary of parameters that override the config file values
        
    Returns:
        Dictionary containing the configuration parameters
        
    Raises:
        FileNotFoundError: If config file is not found
        yaml.YAMLError: If YAML parsing fails
    """
    # Check if config file exists, otherwise use default
    if not os.path.exists(config_file):
        logging.warning(f"Configuration file {config_file} not found. Using default configuration.")
        config_file = os.path.join(os.path.dirname(__file__), "config", "default_config.yaml")
        
        if not os.path.exists(config_file):
            raise FileNotFoundError("Default configuration file not found.")
    
    # Load the configuration from YAML file
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Override configuration parameters if provided
    if override_params:
        config = override_config(config, override_params)
    
    # Set date if not specified
    if not config.get('analysis', {}).get('date'):
        config.setdefault('analysis', {})['date'] = datetime.now().strftime("%Y-%m-%d")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.join(config['paths']['output_dir'], config['analysis']['date'])
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the configuration for reproducibility
    save_config(config, output_dir)
    
    return config


def override_config(config: Dict[str, Any], overrides: Dict[str, Any]) -> Dict[str, Any]:
    """
    Override configuration parameters with provided values.
    
    Args:
        config: Original configuration dictionary
        overrides: Dictionary of override parameters
        
    Returns:
        Updated configuration dictionary
    """
    def update_nested_dict(d: Dict[str, Any], u: Dict[str, Any]) -> Dict[str, Any]:
        for k, v in u.items():
            if isinstance(v, dict):
                d[k] = update_nested_dict(d.get(k, {}), v)
            else:
                d[k] = v
        return d
    
    return update_nested_dict(config.copy(), overrides)


def save_config(config: Dict[str, Any], output_dir: str) -> None:
    """
    Save configuration to output directory for reproducibility.
    
    Args:
        config: Configuration dictionary to save
        output_dir: Directory to save the configuration file
    """
    config_path = os.path.join(output_dir, "analysis_config.yaml")
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    logging.info(f"Configuration saved to: {config_path}")


def setup_logging(config: Dict[str, Any]) -> str:
    """
    Set up logging system based on configuration.
    
    Args:
        config: Configuration dictionary containing logging settings
        
    Returns:
        Path to the log file
    """
    # Set output directory
    output_dir = os.path.join(config['paths']['output_dir'], config['analysis']['date'])
    os.makedirs(output_dir, exist_ok=True)
    
    # Set log file
    log_file = os.path.join(output_dir, config['logging']['file'])
    
    # Configure logging
    log_level = getattr(logging, config['logging']['level'].upper())
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    # Initialize log file
    logging.info("Starting protein neighborhood analysis")
    logging.info(f"Log level set to: {config['logging']['level']}")
    
    return log_file


def create_output_file_path(basepairs: int, max_neighbors: int, 
                          date: str, config: Dict[str, Any]) -> str:
    """
    Create output file path for neighbor data.
    
    Args:
        basepairs: Maximum acceptable intergenic space
        max_neighbors: Maximum number of neighbors to identify
        date: Date string for output directory
        config: Configuration dictionary
        
    Returns:
        Full path to the output file
    """
    output_dir = os.path.join(config['paths']['output_dir'], date)
    filename = f"all_neighbours_bp{basepairs}_n{max_neighbors}.csv"
    return os.path.join(output_dir, filename)


def log_input_file_info(config: Dict[str, Any]) -> None:
    """
    Log information about input files.
    
    Args:
        config: Configuration dictionary
    """
    base_dir = config['paths']['base_dir']
    
    # Core input files
    input_files = [
        os.path.join(base_dir, config['files']['proteins']),
        os.path.join(base_dir, config['files']['assemblies']),
        os.path.join(base_dir, config['files']['protein_assembly']),
    ]
    
    # Representative files if they exist
    if config['files'].get('representative_files'):
        rep_files = [
            os.path.join(base_dir, config['files']['representative_files']['ipg']),
            os.path.join(base_dir, config['files']['representative_files']['pdb']),
            os.path.join(base_dir, config['files']['representative_files']['cluster']),
        ]
        input_files.extend(rep_files)
    
    # Log file information
    for file_path in input_files:
        if os.path.exists(file_path):
            size = os.path.getsize(file_path)
            mtime = datetime.fromtimestamp(os.path.getmtime(file_path))
            logging.info(f"File: {file_path}, Size: {size} bytes, Modified: {mtime}")
        else:
            logging.warning(f"File not found: {file_path}")