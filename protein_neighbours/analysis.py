"""
Analysis functions for protein neighbor data.

This module contains functions for statistical analysis and summary generation
of protein neighborhood data.
"""

import os
import logging
import pandas as pd
import numpy as np
from typing import Dict, Any, List, Optional, Tuple
from collections import Counter


def analyze_neighbor_types(cog_data: pd.DataFrame, output_dir: str) -> Optional[pd.DataFrame]:
    """
    Analyze neighbor types and their counts.
    
    Args:
        cog_data: DataFrame containing COG annotation data
        output_dir: Directory to save output files
        
    Returns:
        DataFrame with neighbor type analysis or None if analysis fails
    """
    os.makedirs(output_dir, exist_ok=True)
    logging.info("Analyzing neighbor types and counts")
    
    # Check if cog_data is valid
    if cog_data is None or cog_data.empty:
        logging.warning("No COG data provided. Cannot analyze neighbor types.")
        return None
    
    # Check required columns
    required_cols = ['Description', 'COG_category']
    if not all(col in cog_data.columns for col in required_cols):
        logging.warning(f"Required columns {required_cols} not found in COG data.")
        logging.warning(f"Available columns: {list(cog_data.columns)}")
        return None
    
    try:
        # Calculate types of neighbors and their counts
        types_of_neighbors = (cog_data
                            .groupby(['Description', 'COG_category'])
                            .size()
                            .reset_index(name='count')
                            .sort_values('count', ascending=False))
        
        logging.info(f"Found {len(types_of_neighbors)} unique neighbor types")
        
        # Save neighbor types data
        output_file = os.path.join(output_dir, 'types_of_neighbors.csv')
        types_of_neighbors.to_csv(output_file, index=False)
        logging.info(f"Saved neighbor types data to: {output_file}")
        
        # Print instruction for manual annotation
        logging.info("Please check the output folder for neighbor types and annotate them manually if needed")
        
        return types_of_neighbors
        
    except Exception as e:
        logging.error(f"Failed to analyze neighbor types: {e}")
        return None


def generate_summary_statistics(combined_data: pd.DataFrame, output_dir: str) -> Dict[str, Any]:
    """
    Generate comprehensive summary statistics of the analysis.
    
    Args:
        combined_data: Combined DataFrame with all analysis data
        output_dir: Directory to save output files
        
    Returns:
        Dictionary containing summary statistics
    """
    os.makedirs(output_dir, exist_ok=True)
    logging.info("Generating summary statistics")
    
    results = {}
    
    try:
        # Count proteins and neighbors
        results['total_proteins'] = len(combined_data[~combined_data['is_neighbour']]['PIGI'].unique())
        results['total_neighbors'] = len(combined_data[combined_data['is_neighbour']])
        results['unique_neighbor_types'] = len(combined_data['COG_category'].unique()) if 'COG_category' in combined_data.columns else 0
        
        # Clade information
        if 'clade' in combined_data.columns:
            results['clades'] = combined_data['clade'].unique().tolist()
            clade_counts = (combined_data[~combined_data['is_neighbour']]
                          .groupby('clade')
                          .size()
                          .reset_index(name='count')
                          .sort_values('count', ascending=False))
            results['clade_counts'] = clade_counts
        else:
            results['clades'] = []
            results['clade_counts'] = pd.DataFrame()
        
        # Assembly information
        if 'assembly' in combined_data.columns:
            results['assemblies'] = combined_data['assembly'].unique().tolist()
            results['assembly_count'] = len(results['assemblies'])
            
            # Proteins per assembly
            proteins_per_assembly = (combined_data[~combined_data['is_neighbour']]
                                   .groupby('assembly')
                                   .size()
                                   .reset_index(name='protein_count'))
            results['proteins_per_assembly'] = proteins_per_assembly
        else:
            results['assemblies'] = []
            results['assembly_count'] = 0
            results['proteins_per_assembly'] = pd.DataFrame()
        
        # Neighbors per protein
        if 'PIGI' in combined_data.columns:
            neighbors_per_protein = (combined_data[combined_data['is_neighbour']]
                                    .groupby('PIGI')
                                    .size()
                                    .reset_index(name='neighbor_count'))
            results['neighbors_per_protein'] = neighbors_per_protein
            results['avg_neighbors_per_protein'] = neighbors_per_protein['neighbor_count'].mean()
        else:
            results['neighbors_per_protein'] = pd.DataFrame()
            results['avg_neighbors_per_protein'] = 0
        
        # Neighbor type distribution
        if 'COG_category' in combined_data.columns:
            neighbor_type_dist = (combined_data[combined_data['is_neighbour']]
                                .groupby('COG_category')
                                .size()
                                .reset_index(name='count')
                                .sort_values('count', ascending=False))
            results['neighbor_type_dist'] = neighbor_type_dist
        else:
            results['neighbor_type_dist'] = pd.DataFrame()
        
        # Create summary DataFrame
        summary_df = pd.DataFrame([{
            'Total_Proteins': results['total_proteins'],
            'Total_Neighbors': results['total_neighbors'],
            'Unique_Neighbor_Types': results['unique_neighbor_types'],
            'Number_of_Clades': len(results['clades']),
            'Number_of_Assemblies': results['assembly_count'],
            'Average_Neighbors_per_Protein': results.get('avg_neighbors_per_protein', 0)
        }])
        
        # Save summary statistics
        output_file = os.path.join(output_dir, "summary_statistics.csv")
        summary_df.to_csv(output_file, index=False)
        logging.info(f"Saved summary statistics to: {output_file}")
        
        # Save detailed statistics
        if not results['clade_counts'].empty:
            output_file = os.path.join(output_dir, "clade_distribution.csv")
            results['clade_counts'].to_csv(output_file, index=False)
        
        if not results['neighbors_per_protein'].empty:
            output_file = os.path.join(output_dir, "neighbors_per_protein.csv")
            results['neighbors_per_protein'].to_csv(output_file, index=False)
        
        if not results['neighbor_type_dist'].empty:
            output_file = os.path.join(output_dir, "neighbor_type_distribution.csv")
            results['neighbor_type_dist'].to_csv(output_file, index=False)
        
        logging.info("Successfully generated summary statistics")
        return results
        
    except Exception as e:
        logging.error(f"Failed to generate summary statistics: {e}")
        return {}


def combine_and_analyze(neighbors_data: pd.DataFrame, 
                       cog_data: Optional[pd.DataFrame],
                       clade_assign: pd.DataFrame,
                       neighbor_annotations: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Combine neighbor data with COG annotations and clade assignments.
    
    Args:
        neighbors_data: DataFrame with neighbor information
        cog_data: DataFrame with COG annotation data
        clade_assign: DataFrame with clade assignments
        neighbor_annotations: Optional manual annotations
        
    Returns:
        Combined DataFrame with all analysis data
    """
    logging.info("Combining data and analyzing")
    
    # Start with neighbors data
    combined_df = neighbors_data.copy()
    
    # Add clade information
    if not clade_assign.empty and 'PIGI' in combined_df.columns:
        combined_df = combined_df.merge(
            clade_assign[['PIGI', 'clade']], 
            on='PIGI', 
            how='left'
        )
        logging.info("Added clade information")
    
    # Add COG annotation information
    if cog_data is not None and not cog_data.empty:
        # Merge on gene/protein ID
        if 'gene' in combined_df.columns and 'query' in cog_data.columns:
            combined_df = combined_df.merge(
                cog_data[['query', 'COG_category', 'Description']],
                left_on='gene',
                right_on='query',
                how='left'
            )
            logging.info("Added COG annotation information")
    
    # Add manual annotations if provided
    if neighbor_annotations is not None and not neighbor_annotations.empty:
        # This would typically merge on COG categories or descriptions
        logging.info("Manual annotations available for integration")
    
    # Fill missing values
    combined_df = combined_df.fillna('Unknown')
    
    logging.info(f"Combined dataset has {len(combined_df)} rows")
    return combined_df


def calculate_correlation_matrix(data: pd.DataFrame, clades: List[str]) -> np.ndarray:
    """
    Calculate correlation matrix of clade co-occurrence.
    
    Args:
        data: DataFrame with assembly and clade information
        clades: List of clade names
        
    Returns:
        Correlation matrix as numpy array
    """
    if data.empty or 'assembly' not in data.columns or 'clade' not in data.columns:
        logging.warning("Insufficient data for correlation matrix calculation")
        return np.array([])
    
    # Create pivot table of assemblies vs clades
    pivot_table = data.pivot_table(
        index='assembly',
        columns='clade', 
        aggfunc='size',
        fill_value=0
    )
    
    # Calculate correlation matrix
    correlation_matrix = pivot_table.corr().values
    
    logging.info("Calculated correlation matrix")
    return correlation_matrix


def analyze_protein_distribution(combined_data: pd.DataFrame) -> Dict[str, Any]:
    """
    Analyze protein distribution across clades and assemblies.
    
    Args:
        combined_data: Combined DataFrame with all analysis data
        
    Returns:
        Dictionary with distribution analysis results
    """
    results = {}
    
    if 'clade' in combined_data.columns and 'assembly' in combined_data.columns:
        # Proteins per clade
        proteins_per_clade = (combined_data[~combined_data['is_neighbour']]
                            .groupby('clade')
                            .size()
                            .to_dict())
        results['proteins_per_clade'] = proteins_per_clade
        
        # Assemblies per clade
        assemblies_per_clade = (combined_data[~combined_data['is_neighbour']]
                              .groupby('clade')['assembly']
                              .nunique()
                              .to_dict())
        results['assemblies_per_clade'] = assemblies_per_clade
        
        # Proteins per assembly
        proteins_per_assembly = (combined_data[~combined_data['is_neighbour']]
                               .groupby('assembly')
                               .size()
                               .describe()
                               .to_dict())
        results['proteins_per_assembly_stats'] = proteins_per_assembly
        
        logging.info("Analyzed protein distribution")
    
    return results