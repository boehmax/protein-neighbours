"""
Neighbor identification functions for protein genomic analysis.

This module contains functions for identifying neighboring proteins
in genomic data, parsing GFF files, and extracting attributes.
"""

import os
import logging
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import re


def get_attribute_field(attributes: str, field: str, attr_sep: str = ";") -> Optional[str]:
    """
    Extract a specific field from the attribute column of a GFF file.
    
    Args:
        attributes: The attribute column as a string
        field: The field to extract
        attr_sep: The separator used in the attribute column
        
    Returns:
        The extracted field value or None if not found
    """
    if not attributes or pd.isna(attributes):
        return None
        
    # Split the attributes
    attrs = attributes.split(attr_sep)
    
    # Look for the desired field
    for attr in attrs:
        if '=' in attr:
            key, value = attr.split('=', 1)
            if key.strip() == field:
                return value.strip()
    
    return None


def parse_gff_file(gff_file: str) -> pd.DataFrame:
    """
    Parse a GFF3 file and return a DataFrame.
    
    Args:
        gff_file: Path to the GFF3 file
        
    Returns:
        DataFrame containing GFF3 data
        
    Raises:
        FileNotFoundError: If GFF file is not found
    """
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")
    
    try:
        # Read GFF file
        gff_data = []
        with open(gff_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        gff_data.append({
                            'seqid': parts[0],
                            'source': parts[1],
                            'type': parts[2],
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'score': parts[5] if parts[5] != '.' else None,
                            'strand': parts[6],
                            'phase': parts[7] if parts[7] != '.' else None,
                            'attributes': parts[8]
                        })
        
        df = pd.DataFrame(gff_data)
        
        # Extract ID and Name from attributes
        if not df.empty:
            df['ID'] = df['attributes'].apply(lambda x: get_attribute_field(x, 'ID'))
            df['Name'] = df['attributes'].apply(lambda x: get_attribute_field(x, 'Name'))
            
            # Clean ID field (remove prefix if present)
            df['ID'] = df['ID'].apply(lambda x: re.sub(r'.*-', '', x) if x else x)
        
        logging.debug(f"Successfully parsed GFF file: {gff_file}")
        return df
        
    except Exception as e:
        logging.error(f"Failed to parse GFF file {gff_file}: {e}")
        raise


def get_protein_neighbors_from_gff(gff_file: str, protein_id: str, 
                                 basepairs: int = 300, max_neighbors: int = 15,
                                 overlap: int = 50) -> pd.DataFrame:
    """
    Get neighboring proteins from a GFF3 file for a specific protein.
    
    Args:
        gff_file: Path to the GFF3 file
        protein_id: Target protein ID
        basepairs: Maximum acceptable intergenic space
        max_neighbors: Maximum number of neighbors to identify
        overlap: Acceptable overlap of two adjacent genes
        
    Returns:
        DataFrame containing neighbor information
    """
    try:
        # Parse GFF file
        df = parse_gff_file(gff_file)
        
        if df.empty:
            logging.warning(f"Empty GFF file: {gff_file}")
            return pd.DataFrame()
        
        # Find the target protein
        target_protein = df[df['ID'] == protein_id]
        if target_protein.empty:
            logging.warning(f"Protein ID {protein_id} not found in GFF file: {gff_file}")
            return pd.DataFrame()
        
        target = target_protein.iloc[0]
        
        # Filter for CDS features on the same sequence
        cds_features = df[(df['type'] == 'CDS') & (df['seqid'] == target['seqid'])].copy()
        
        if cds_features.empty:
            logging.warning(f"No CDS features found in {gff_file}")
            return pd.DataFrame()
        
        # Calculate distances and identify neighbors
        neighbors = []
        target_start, target_end = target['start'], target['end']
        
        for _, feature in cds_features.iterrows():
            if feature['ID'] == protein_id:
                continue  # Skip the target protein itself
                
            # Calculate distance
            if feature['end'] < target_start:
                # Upstream neighbor
                distance = target_start - feature['end']
                position = 'upstream'
            elif feature['start'] > target_end:
                # Downstream neighbor
                distance = feature['start'] - target_end
                position = 'downstream'
            else:
                # Overlapping
                distance = 0
                position = 'overlapping'
            
            # Check if within distance threshold
            if distance <= basepairs or position == 'overlapping':
                neighbors.append({
                    'molecule': feature['seqid'],
                    'gene': feature['ID'],
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': feature['strand'],
                    'distance': distance,
                    'position': position,
                    'Name': feature.get('Name', ''),
                    'is_neighbour': True
                })
        
        # Sort by distance and limit to max_neighbors
        neighbors.sort(key=lambda x: x['distance'])
        neighbors = neighbors[:max_neighbors]
        
        # Create DataFrame
        neighbor_df = pd.DataFrame(neighbors)
        
        # Add target protein info
        if not neighbor_df.empty:
            target_info = {
                'molecule': target['seqid'],
                'gene': protein_id,
                'start': target['start'],
                'end': target['end'],
                'strand': target['strand'],
                'distance': 0,
                'position': 'target',
                'Name': target.get('Name', ''),
                'is_neighbour': False
            }
            target_df = pd.DataFrame([target_info])
            neighbor_df = pd.concat([neighbor_df, target_df], ignore_index=True)
        
        logging.debug(f"Found {len(neighbors)} neighbors for protein {protein_id}")
        return neighbor_df
        
    except Exception as e:
        logging.error(f"Error getting neighbors for protein {protein_id}: {e}")
        return pd.DataFrame()


def collect_all_neighbors(protein_assembly: pd.DataFrame, basepairs: int = 300,
                         max_neighbors: int = 15, base_path: str = "data",
                         overlap: int = 50) -> pd.DataFrame:
    """
    Collect neighbor information for all proteins.
    
    Args:
        protein_assembly: DataFrame with protein-assembly mappings
        basepairs: Maximum acceptable intergenic space  
        max_neighbors: Maximum number of neighbors to identify
        base_path: Base path to data directory
        overlap: Acceptable overlap of two adjacent genes
        
    Returns:
        DataFrame containing all neighbor information
    """
    all_neighbors = []
    total_proteins = len(protein_assembly)
    
    logging.info(f"Processing {total_proteins} proteins for neighbor identification")
    
    for i, (_, row) in enumerate(protein_assembly.iterrows()):
        # Show progress
        if i % max(1, total_proteins // 20) == 0 or i == total_proteins - 1:
            progress_pct = round((i + 1) / total_proteins * 100)
            logging.info(f"Processing protein {i + 1} of {total_proteins} ({progress_pct}%)")
        
        # Construct GFF file path
        gff_file = os.path.join(base_path, 'ncbi_dataset/data', 
                               row['assembly_id'], 'genomic.gff')
        
        # Get neighbors for this protein
        neighbors = get_protein_neighbors_from_gff(
            gff_file, row['protein_id'], basepairs, max_neighbors, overlap
        )
        
        if not neighbors.empty:
            # Add assembly information
            neighbors['assembly'] = row['assembly_id']
            neighbors['PIGI'] = row['protein_id']
            all_neighbors.append(neighbors)
    
    if all_neighbors:
        result = pd.concat(all_neighbors, ignore_index=True)
        logging.info(f"Collected neighbor information for {len(result)} proteins")
        return result
    else:
        logging.warning("No neighbor information collected")
        return pd.DataFrame()


def collect_all_protein_info(protein_assembly: pd.DataFrame, 
                           base_path: str = "data") -> pd.DataFrame:
    """
    Collect general protein information from GFF files.
    
    Args:
        protein_assembly: DataFrame with protein-assembly mappings
        base_path: Base path to data directory
        
    Returns:
        DataFrame containing protein information
    """
    all_proteins = []
    total_proteins = len(protein_assembly)
    
    logging.info(f"Collecting protein information for {total_proteins} proteins")
    
    for i, (_, row) in enumerate(protein_assembly.iterrows()):
        # Show progress
        if i % max(1, total_proteins // 20) == 0 or i == total_proteins - 1:
            progress_pct = round((i + 1) / total_proteins * 100)
            logging.info(f"Processing protein {i + 1} of {total_proteins} ({progress_pct}%)")
        
        # Construct GFF file path
        gff_file = os.path.join(base_path, 'ncbi_dataset/data',
                               row['assembly_id'], 'genomic.gff')
        
        try:
            # Parse GFF and find the protein
            df = parse_gff_file(gff_file)
            if not df.empty:
                protein_info = df[df['ID'] == row['protein_id']]
                if not protein_info.empty:
                    protein = protein_info.iloc[0]
                    protein_data = {
                        'molecule': protein['seqid'],
                        'gene': protein['ID'],
                        'start': protein['start'],
                        'end': protein['end'],
                        'strand': protein['strand'],
                        'Name': protein.get('Name', ''),
                        'assembly': row['assembly_id'],
                        'PIGI': row['protein_id'],
                        'is_neighbour': False
                    }
                    all_proteins.append(protein_data)
        except Exception as e:
            logging.warning(f"Could not process protein {row['protein_id']}: {e}")
    
    if all_proteins:
        result = pd.DataFrame(all_proteins)
        logging.info(f"Collected information for {len(result)} proteins")
        return result
    else:
        logging.warning("No protein information collected")
        return pd.DataFrame()