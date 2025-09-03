"""
Input/Output functions for protein neighbor analysis.

This module contains functions for reading and writing data files,
including protein and assembly data, clades, and representative proteins.
"""

import os
import logging
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import glob


def read_protein_assembly_data(protein_file: str = "proteins.csv",
                             assembly_file: str = "assm_accs.csv",
                             protein_assembly_file: str = "assm_accs_protein.csv",
                             base_path: str = "data",
                             interactive: bool = True,
                             protein_of_interest: Optional[str] = None) -> Dict[str, Any]:
    """
    Read protein and assembly data from specified files.
    
    Args:
        protein_file: Path to the protein file (CSV with protein IDs)
        assembly_file: Path to the assembly file (CSV with assembly IDs)
        protein_assembly_file: Path to the protein assembly file (maps proteins to assemblies)
        base_path: Base path to the data directory
        interactive: Whether to prompt the user for input
        protein_of_interest: Optional protein ID to use
        
    Returns:
        Dictionary containing the protein, assembly, and protein assembly DataFrames
        
    Raises:
        FileNotFoundError: If required files are not found
    """
    logging.info(f"Reading protein and assembly data from: {base_path}")
    
    # Construct full paths
    protein_path = os.path.join(base_path, protein_file)
    assembly_path = os.path.join(base_path, assembly_file) 
    protein_assembly_path = os.path.join(base_path, protein_assembly_file)
    
    # Check if files exist
    for path, name in [(protein_path, "Protein"), (assembly_path, "Assembly"), 
                       (protein_assembly_path, "Protein-assembly mapping")]:
        if not os.path.exists(path):
            logging.error(f"{name} file not found: {path}")
            raise FileNotFoundError(f"{name} file not found: {path}")
    
    # Read files
    try:
        protein = pd.read_csv(protein_path, header=None)
        logging.info(f"Read {len(protein)} proteins from {protein_path}")
    except Exception as e:
        logging.error(f"Failed to read protein file: {e}")
        raise
    
    try:
        assembly = pd.read_csv(assembly_path, header=None)
        logging.info(f"Read {len(assembly)} assemblies from {assembly_path}")
    except Exception as e:
        logging.error(f"Failed to read assembly file: {e}")
        raise
        
    try:
        protein_assembly = pd.read_csv(protein_assembly_path, header=None)
        logging.info(f"Read {len(protein_assembly)} protein-assembly mappings from {protein_assembly_path}")
    except Exception as e:
        logging.error(f"Failed to read protein-assembly mapping file: {e}")
        raise
    
    # Set column names
    protein.columns = ['protein_id']
    assembly.columns = ['assembly_id']
    protein_assembly.columns = ['assembly_id', 'protein_id']
    
    # Handle protein of interest
    if not protein_of_interest and interactive:
        print("Available proteins:")
        print(protein['protein_id'].head(10).to_string(index=False))
        if len(protein) > 10:
            print(f"... and {len(protein) - 10} more")
        protein_of_interest = input("Enter protein of interest (or press Enter to skip): ").strip()
        if not protein_of_interest:
            protein_of_interest = None
    
    return {
        'protein': protein,
        'assembly': assembly,
        'protein_assembly': protein_assembly,
        'protein_of_interest': protein_of_interest
    }


def read_clades(base_path: str = "data", 
               clade_dir: str = "clades",
               pattern: str = "^[C]") -> pd.DataFrame:
    """
    Read clade assignment files.
    
    Args:
        base_path: Base path to the data directory
        clade_dir: Directory containing clade files
        pattern: Pattern to match clade files
        
    Returns:
        DataFrame with protein IDs and their clade assignments
    """
    clade_path = os.path.join(base_path, clade_dir)
    
    if not os.path.exists(clade_path):
        logging.warning(f"Clade directory not found: {clade_path}")
        return pd.DataFrame(columns=['protein_id', 'clade', 'PIGI'])
    
    # Find clade files matching pattern
    clade_files = glob.glob(os.path.join(clade_path, f"{pattern}*"))
    
    if not clade_files:
        logging.warning(f"No clade files found matching pattern: {pattern}")
        return pd.DataFrame(columns=['protein_id', 'clade', 'PIGI'])
    
    logging.info(f"Found {len(clade_files)} clade files")
    
    all_clades = []
    for i, clade_file in enumerate(clade_files):
        try:
            # Read clade file
            clade_data = pd.read_csv(clade_file, header=None, names=['protein_id'])
            clade_data['clade'] = chr(ord('A') + i)  # Assign clade letters A, B, C, etc.
            clade_data['PIGI'] = clade_data['protein_id']  # PIGI same as protein_id
            all_clades.append(clade_data)
            logging.debug(f"Processed clade file: {os.path.basename(clade_file)}")
        except Exception as e:
            logging.error(f"Failed to read clade file {clade_file}: {e}")
    
    if all_clades:
        result = pd.concat(all_clades, ignore_index=True)
        logging.info(f"Processed {len(result)} entries from clade files")
        return result
    else:
        return pd.DataFrame(columns=['protein_id', 'clade', 'PIGI'])


def read_representatives(base_path: str = "data",
                        rep_path: str = "representatives",
                        ipg_file: str = "ipg_representative.txt",
                        pdb_file: str = "pdb_representative.txt", 
                        cluster_file: str = "cluster_representative.txt") -> pd.DataFrame:
    """
    Read representative protein files.
    
    Args:
        base_path: Base path to the data directory
        rep_path: Path to representatives directory
        ipg_file: IPG representative file name
        pdb_file: PDB representative file name
        cluster_file: Cluster representative file name
        
    Returns:
        DataFrame with representative protein information
    """
    rep_dir = os.path.join(base_path, rep_path)
    
    if not os.path.exists(rep_dir):
        logging.warning(f"Representatives directory not found: {rep_dir}")
        return pd.DataFrame(columns=['alias', 'representative'])
    
    representatives = []
    
    # Read each representative file if it exists
    for filename, file_type in [(ipg_file, 'ipg'), (pdb_file, 'pdb'), (cluster_file, 'cluster')]:
        file_path = os.path.join(rep_dir, filename)
        if os.path.exists(file_path):
            try:
                rep_data = pd.read_csv(file_path, sep='\t', header=None)
                if len(rep_data.columns) >= 2:
                    rep_data.columns = ['alias', 'representative']
                    rep_data['type'] = file_type
                    representatives.append(rep_data)
                    logging.info(f"Read {len(rep_data)} {file_type} representatives")
                else:
                    logging.warning(f"Invalid format in {file_path}")
            except Exception as e:
                logging.error(f"Failed to read {file_path}: {e}")
        else:
            logging.warning(f"Representative file not found: {file_path}")
    
    if representatives:
        result = pd.concat(representatives, ignore_index=True)
        logging.info(f"Total representatives loaded: {len(result)}")
        return result
    else:
        return pd.DataFrame(columns=['alias', 'representative', 'type'])


def read_annotations(date: str, base_path: str = "output") -> Optional[pd.DataFrame]:
    """
    Read manual annotation file if it exists.
    
    Args:
        date: Date string for the output directory
        base_path: Base path to output directory
        
    Returns:
        DataFrame with manual annotations or None if file doesn't exist
    """
    annotation_file = os.path.join(base_path, date, "types_of_neighbours_annotated.csv")
    
    if os.path.exists(annotation_file):
        try:
            annotations = pd.read_csv(annotation_file)
            logging.info(f"Read {len(annotations)} manual annotations")
            
            # Replace NA values if needed
            annotations = annotations.fillna("Unknown")
            
            # Rename columns to expected names
            if len(annotations.columns) >= 4:
                annotations.columns = ['COG_NAME', 'COG_LETTER', 'N', 'ANNOTATION']
            
            return annotations
        except Exception as e:
            logging.error(f"Failed to read annotation file: {e}")
            return None
    else:
        logging.info("No manual annotation file found")
        return None


def save_data(data: pd.DataFrame, filename: str, output_dir: str) -> None:
    """
    Save DataFrame to CSV file.
    
    Args:
        data: DataFrame to save
        filename: Name of the output file
        output_dir: Output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    data.to_csv(output_path, index=False)
    logging.info(f"Saved data to: {output_path}")