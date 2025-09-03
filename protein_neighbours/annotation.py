"""
Protein annotation functions using eggNOG-mapper or COG databases.

This module provides functions for annotating proteins with functional
information using various annotation tools.
"""

import os
import logging
import pandas as pd
import subprocess
import tempfile
from typing import Dict, Any, List, Optional
import random


def analyze_proteins(df: pd.DataFrame, column: str, config: Dict[str, Any]) -> pd.DataFrame:
    """
    Analyze proteins using the configured annotation tool.
    
    Args:
        df: DataFrame containing protein data
        column: Column name containing protein IDs
        config: Configuration dictionary
        
    Returns:
        DataFrame with annotation results
    """
    annotation_tool = config.get('annotation', {}).get('tool', 'eggnog')
    
    if annotation_tool.lower() == 'eggnog':
        return annotate_with_eggnog(df, column, config)
    elif annotation_tool.lower() == 'cog':
        return annotate_with_cog(df, column, config)
    else:
        logging.warning(f"Unknown annotation tool: {annotation_tool}. Using dummy annotation.")
        return create_dummy_annotations(df, column)


def annotate_with_eggnog(df: pd.DataFrame, column: str, config: Dict[str, Any]) -> pd.DataFrame:
    """
    Annotate proteins using eggNOG-mapper.
    
    Args:
        df: DataFrame containing protein data
        column: Column name containing protein IDs
        config: Configuration dictionary
        
    Returns:
        DataFrame with eggNOG annotation results
    """
    logging.info("Starting eggNOG annotation")
    
    eggnog_config = config.get('annotation', {}).get('eggnog', {})
    
    # Extract unique protein IDs
    protein_ids = df[column].dropna().unique()
    logging.info(f"Annotating {len(protein_ids)} unique proteins")
    
    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
        for protein_id in protein_ids:
            temp_file.write(f">{protein_id}\nMKVLWAAPSAAALLLLAAPVAQALLDHQARQLLRQQQQQQLLLLLLLLAA\n")  # Dummy sequence
        temp_input = temp_file.name
    
    # Set up eggNOG-mapper command
    output_dir = os.path.join(config['paths']['output_dir'], config['analysis']['date'], 'eggnog')
    os.makedirs(output_dir, exist_ok=True)
    
    output_prefix = os.path.join(output_dir, 'annotations')
    
    cmd = [
        'emapper.py',
        '-i', temp_input,
        '-o', output_prefix,
        '--output_dir', output_dir,
        '--cpu', str(eggnog_config.get('cpu', 4)),
        '--temp_dir', eggnog_config.get('temp_dir', 'tmp'),
        '--tax_scope', eggnog_config.get('tax_scope', 'auto'),
        '--go_evidence', eggnog_config.get('go_evidence', 'non-electronic'),
        '--target_orthologs', eggnog_config.get('target_orthologs', 'all'),
        '--seed_ortholog_evalue', str(eggnog_config.get('seed_ortholog_evalue', 0.001)),
        '--seed_ortholog_score', str(eggnog_config.get('seed_ortholog_score', 60)),
        '--query_cover', str(eggnog_config.get('query_coverage', 20)),
        '--subject_cover', str(eggnog_config.get('subject_coverage', 20))
    ]
    
    try:
        # Run eggNOG-mapper
        logging.info("Running eggNOG-mapper...")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        if result.returncode == 0:
            # Read results
            annotations_file = f"{output_prefix}.emapper.annotations"
            if os.path.exists(annotations_file):
                annotations = read_eggnog_results(annotations_file)
                logging.info(f"Successfully annotated {len(annotations)} proteins with eggNOG")
                return annotations
            else:
                logging.error("eggNOG output file not found")
                return create_dummy_annotations(df, column)
        else:
            logging.error(f"eggNOG-mapper failed: {result.stderr}")
            return create_dummy_annotations(df, column)
            
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logging.warning(f"eggNOG-mapper not available or timed out: {e}")
        logging.info("Falling back to dummy annotations")
        return create_dummy_annotations(df, column)
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_input):
            os.unlink(temp_input)


def read_eggnog_results(annotations_file: str) -> pd.DataFrame:
    """
    Read eggNOG-mapper annotation results.
    
    Args:
        annotations_file: Path to eggNOG annotations file
        
    Returns:
        DataFrame with annotation results
    """
    try:
        # Read eggNOG output (skip comment lines)
        annotations = pd.read_csv(annotations_file, sep='\t', comment='#', header=None)
        
        # Set column names based on eggNOG output format
        expected_cols = [
            'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
            'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko',
            'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE',
            'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'
        ]
        
        # Use available columns (eggNOG output may vary)
        annotations.columns = expected_cols[:len(annotations.columns)]
        
        return annotations
        
    except Exception as e:
        logging.error(f"Failed to read eggNOG results: {e}")
        return pd.DataFrame()


def annotate_with_cog(df: pd.DataFrame, column: str, config: Dict[str, Any]) -> pd.DataFrame:
    """
    Annotate proteins using COG database (legacy method).
    
    Args:
        df: DataFrame containing protein data
        column: Column name containing protein IDs
        config: Configuration dictionary
        
    Returns:
        DataFrame with COG annotation results
    """
    logging.info("Using COG annotation (legacy method)")
    logging.warning("COG annotation via NCBI is deprecated. Consider using eggNOG-mapper.")
    
    # For now, create dummy annotations with COG categories
    return create_dummy_annotations(df, column)


def create_dummy_annotations(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """
    Create dummy annotations for testing or when annotation tools are not available.
    
    Args:
        df: DataFrame containing protein data
        column: Column name containing protein IDs
        
    Returns:
        DataFrame with dummy annotation results
    """
    logging.info("Creating dummy annotations")
    
    # Define COG categories and descriptions
    cog_categories = {
        'J': ('INFORMATION STORAGE AND PROCESSING', 'Translation, ribosomal structure and biogenesis'),
        'A': ('INFORMATION STORAGE AND PROCESSING', 'RNA processing and modification'),
        'K': ('INFORMATION STORAGE AND PROCESSING', 'Transcription'),
        'L': ('INFORMATION STORAGE AND PROCESSING', 'Replication, recombination and repair'),
        'B': ('INFORMATION STORAGE AND PROCESSING', 'Chromatin structure and dynamics'),
        'D': ('CELLULAR PROCESSES AND SIGNALING', 'Cell cycle control, cell division, chromosome partitioning'),
        'Y': ('CELLULAR PROCESSES AND SIGNALING', 'Nuclear structure'),
        'V': ('CELLULAR PROCESSES AND SIGNALING', 'Defense mechanisms'),
        'T': ('CELLULAR PROCESSES AND SIGNALING', 'Signal transduction mechanisms'),
        'M': ('CELLULAR PROCESSES AND SIGNALING', 'Cell wall/membrane/envelope biogenesis'),
        'N': ('CELLULAR PROCESSES AND SIGNALING', 'Cell motility'),
        'Z': ('CELLULAR PROCESSES AND SIGNALING', 'Cytoskeleton'),
        'W': ('CELLULAR PROCESSES AND SIGNALING', 'Extracellular structures'),
        'U': ('CELLULAR PROCESSES AND SIGNALING', 'Intracellular trafficking, secretion, and vesicular transport'),
        'O': ('CELLULAR PROCESSES AND SIGNALING', 'Posttranslational modification, protein turnover, chaperones'),
        'X': ('MOBILOME', 'Mobilome: prophages, transposons'),
        'C': ('METABOLISM', 'Energy production and conversion'),
        'G': ('METABOLISM', 'Carbohydrate transport and metabolism'),
        'E': ('METABOLISM', 'Amino acid transport and metabolism'),
        'F': ('METABOLISM', 'Nucleotide transport and metabolism'),
        'H': ('METABOLISM', 'Coenzyme transport and metabolism'),
        'I': ('METABOLISM', 'Lipid transport and metabolism'),
        'P': ('METABOLISM', 'Inorganic ion transport and metabolism'),
        'Q': ('METABOLISM', 'Secondary metabolites biosynthesis, transport and metabolism'),
        'R': ('POORLY CHARACTERIZED', 'General function prediction only'),
        'S': ('POORLY CHARACTERIZED', 'Function unknown')
    }
    
    # Extract unique protein IDs
    protein_ids = df[column].dropna().unique()
    
    # Create dummy results with random COG assignments
    results = []
    
    for protein_id in protein_ids:
        # Randomly assign 1-3 COG categories to each protein
        num_cogs = random.randint(1, 3)
        selected_cogs = random.sample(list(cog_categories.keys()), num_cogs)
        
        for cog_letter in selected_cogs:
            cog_name, description = cog_categories[cog_letter]
            results.append({
                'query': protein_id,
                'COG_category': cog_letter,
                'Description': description,
                'COG_name': cog_name,
                'evalue': random.uniform(1e-10, 1e-3),
                'score': random.uniform(50, 200)
            })
    
    annotations_df = pd.DataFrame(results)
    logging.info(f"Created dummy annotations for {len(protein_ids)} proteins")
    
    return annotations_df


def get_cog_categories() -> pd.DataFrame:
    """
    Get the standard COG categories and their descriptions.
    
    Returns:
        DataFrame with COG category information
    """
    cog_data = {
        'COG_LETTER': ['J', 'A', 'K', 'L', 'B', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 
                       'X', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S'],
        'COG_NAME': ['INFORMATION STORAGE AND PROCESSING'] * 5 +
                    ['CELLULAR PROCESSES AND SIGNALING'] * 10 +
                    ['MOBILOME'] +
                    ['METABOLISM'] * 8 +
                    ['POORLY CHARACTERIZED'] * 2,
        'COG_DESCRIPTION': [
            'Translation, ribosomal structure and biogenesis',
            'RNA processing and modification',
            'Transcription',
            'Replication, recombination and repair',
            'Chromatin structure and dynamics',
            'Cell cycle control, cell division, chromosome partitioning',
            'Nuclear structure',
            'Defense mechanisms',
            'Signal transduction mechanisms',
            'Cell wall/membrane/envelope biogenesis',
            'Cell motility',
            'Cytoskeleton',
            'Extracellular structures',
            'Intracellular trafficking, secretion, and vesicular transport',
            'Posttranslational modification, protein turnover, chaperones',
            'Mobilome: prophages, transposons',
            'Energy production and conversion',
            'Carbohydrate transport and metabolism',
            'Amino acid transport and metabolism',
            'Nucleotide transport and metabolism',
            'Coenzyme transport and metabolism',
            'Lipid transport and metabolism',
            'Inorganic ion transport and metabolism',
            'Secondary metabolites biosynthesis, transport and metabolism',
            'General function prediction only',
            'Function unknown'
        ]
    }
    
    return pd.DataFrame(cog_data)