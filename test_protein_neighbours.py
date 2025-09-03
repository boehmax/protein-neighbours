"""Test suite for protein-neighbours Python package."""

import pytest
import pandas as pd
import tempfile
import os
from unittest.mock import patch, MagicMock

# Import modules to test
from protein_neighbours.utils import load_config, setup_logging
from protein_neighbours.io import read_protein_assembly_data
from protein_neighbours.analysis import analyze_neighbor_types


def test_load_config():
    """Test configuration loading."""
    # Create a temporary config file
    config_content = """
    paths:
      base_dir: "test_data"
      output_dir: "test_output"
    analysis:
      basepairs: 300
      max_neighbors: 15
    """
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(config_content)
        config_file = f.name
    
    try:
        config = load_config(config_file)
        assert config['paths']['base_dir'] == 'test_data'
        assert config['analysis']['basepairs'] == 300
        assert 'date' in config['analysis']  # Should be auto-added
    finally:
        os.unlink(config_file)


def test_read_protein_assembly_data():
    """Test protein assembly data reading."""
    # Create temporary CSV files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test protein file
        protein_file = os.path.join(temp_dir, 'proteins.csv')
        with open(protein_file, 'w') as f:
            f.write("WP_123456789.1\nWP_987654321.1\n")
        
        # Create test assembly file
        assembly_file = os.path.join(temp_dir, 'assemblies.csv')
        with open(assembly_file, 'w') as f:
            f.write("GCA_123456789.1\nGCA_987654321.1\n")
        
        # Create test protein-assembly file
        protein_assembly_file = os.path.join(temp_dir, 'protein_assembly.csv')
        with open(protein_assembly_file, 'w') as f:
            f.write("GCA_123456789.1,WP_123456789.1\nGCA_987654321.1,WP_987654321.1\n")
        
        # Test reading
        data = read_protein_assembly_data(
            protein_file='proteins.csv',
            assembly_file='assemblies.csv',
            protein_assembly_file='protein_assembly.csv',
            base_path=temp_dir,
            interactive=False
        )
        
        assert len(data['protein']) == 2
        assert len(data['assembly']) == 2
        assert len(data['protein_assembly']) == 2
        assert 'protein_id' in data['protein'].columns
        assert 'assembly_id' in data['assembly'].columns


def test_analyze_neighbor_types():
    """Test neighbor type analysis."""
    # Create test COG data
    cog_data = pd.DataFrame({
        'Description': ['Protein A', 'Protein B', 'Protein A'],
        'COG_category': ['J', 'K', 'J']
    })
    
    with tempfile.TemporaryDirectory() as temp_dir:
        result = analyze_neighbor_types(cog_data, temp_dir)
        
        assert result is not None
        assert len(result) == 2  # Two unique combinations
        assert 'count' in result.columns
        assert result.loc[result['Description'] == 'Protein A', 'count'].iloc[0] == 2


def test_package_import():
    """Test that the main package can be imported."""
    import protein_neighbours
    assert hasattr(protein_neighbours, '__version__')


if __name__ == "__main__":
    pytest.main([__file__])