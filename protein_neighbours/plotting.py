"""
Visualization functions for protein neighbor analysis.

This module provides functions for creating plots and visualizations
of protein neighborhood data using matplotlib and seaborn.
"""

import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Any, List, Optional, Tuple
import matplotlib.patches as patches


def plot_neighbors_per_clade(combined_data: pd.DataFrame, 
                           exclude_unknown_clade: bool = False,
                           exclude_unknown_cog: bool = False,
                           output_path: Optional[str] = None,
                           plot_count_codh: bool = False,
                           width: int = 30, height: int = 15) -> None:
    """
    Plot the distribution of neighbors by clade.
    
    Args:
        combined_data: Combined DataFrame with all analysis data
        exclude_unknown_clade: Whether to exclude unknown clades
        exclude_unknown_cog: Whether to exclude unknown COGs
        output_path: Optional output subdirectory
        plot_count_codh: Whether to plot CODH counts
        width: Plot width in inches
        height: Plot height in inches
    """
    logging.info("Plotting neighbors per clade")
    
    # Set output directory
    current_date = pd.Timestamp.now().strftime("%Y-%m-%d")
    if output_path is None:
        output_dir = os.path.join("output", current_date)
    else:
        output_dir = os.path.join("output", current_date, output_path)
    
    # Add subdirectory for excluded unknowns if needed
    if exclude_unknown_clade or exclude_unknown_cog:
        output_dir = os.path.join(output_dir, "exclude_unknowns")
    
    # Add subdirectory for CODH count if needed
    if plot_count_codh:
        output_dir = os.path.join(output_dir, "count_codh")
    
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Prepare data for plotting
        plot_data = combined_data.copy()
        
        # Filter data if requested
        if exclude_unknown_clade and 'clade' in plot_data.columns:
            plot_data = plot_data[plot_data['clade'] != 'Unknown']
        
        if exclude_unknown_cog and 'COG_category' in plot_data.columns:
            plot_data = plot_data[plot_data['COG_category'] != 'Unknown']
        
        if plot_data.empty:
            logging.warning("No data remaining after filtering")
            return
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(width, height))
        
        if 'clade' in plot_data.columns and 'COG_category' in plot_data.columns:
            # Create a heatmap-style plot
            pivot_data = plot_data.groupby(['clade', 'COG_category']).size().unstack(fill_value=0)
            
            sns.heatmap(pivot_data, annot=True, fmt='d', ax=ax, cmap='viridis')
            ax.set_title('Neighbor Distribution by Clade and COG Category')
            ax.set_xlabel('COG Category')
            ax.set_ylabel('Clade')
        else:
            # Simple bar plot
            if 'clade' in plot_data.columns:
                clade_counts = plot_data['clade'].value_counts()
                clade_counts.plot(kind='bar', ax=ax)
                ax.set_title('Distribution by Clade')
                ax.set_xlabel('Clade')
                ax.set_ylabel('Count')
            else:
                logging.warning("No clade information available for plotting")
                return
        
        plt.tight_layout()
        
        # Save the plot
        plot_file = os.path.join(output_dir, "neighbor_distribution_by_clade.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Saved plot to: {plot_file}")
        
        # Save prepared data
        if 'clade' in plot_data.columns:
            neighbor_count = plot_data.groupby(['clade', 'COG_category']).size().reset_index(name='count')
            neighbor_count.to_csv(os.path.join(output_dir, "neighbor_count_per_clade.csv"), index=False)
        
    except Exception as e:
        logging.error(f"Failed to create neighbor distribution plot: {e}")


def create_correlation_matrix(data: pd.DataFrame, clades: List[str], 
                            output_dir: str = "output",
                            width: int = 10, height: int = 8,
                            suppress_output: bool = False) -> Optional[plt.Figure]:
    """
    Create a correlation matrix plot showing clade co-occurrence.
    
    Args:
        data: DataFrame with assembly and clade information
        clades: List of clade names
        output_dir: Output directory for saving plots
        width: Plot width in inches
        height: Plot height in inches
        suppress_output: Whether to suppress plot display
        
    Returns:
        Figure object if not suppressed, None otherwise
    """
    logging.info("Creating correlation matrix")
    
    if data.empty or 'assembly' not in data.columns or 'clade' not in data.columns:
        logging.warning("Insufficient data for correlation matrix")
        return None
    
    try:
        # Create pivot table
        pivot_table = data.pivot_table(
            index='assembly',
            columns='clade',
            aggfunc='size',
            fill_value=0
        )
        
        # Calculate correlation matrix
        correlation_matrix = pivot_table.corr()
        
        # Create plot
        fig, ax = plt.subplots(figsize=(width, height))
        
        # Create heatmap
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0,
                   square=True, ax=ax, cbar_kws={'label': 'Correlation'})
        
        ax.set_title('Clade Co-occurrence Correlation Matrix')
        plt.tight_layout()
        
        # Save plot
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, 'correlation_matrix.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logging.info(f"Saved correlation matrix plot to: {output_file}")
        
        # Save matrix as CSV
        output_file_csv = os.path.join(output_dir, 'correlation_matrix.csv')
        correlation_matrix.to_csv(output_file_csv)
        logging.info(f"Saved correlation matrix data to: {output_file_csv}")
        
        if not suppress_output:
            return fig
        else:
            plt.close(fig)
            return None
            
    except Exception as e:
        logging.error(f"Failed to create correlation matrix: {e}")
        return None


def create_clade_histograms(data: pd.DataFrame, output_dir: str = "output",
                          clade_colors: Optional[List[str]] = None) -> None:
    """
    Create histograms for each clade showing protein distribution.
    
    Args:
        data: DataFrame with protein and clade information
        output_dir: Output directory for saving plots
        clade_colors: Optional list of colors for clades
    """
    logging.info("Creating clade histograms")
    
    if data.empty or 'clade' not in data.columns:
        logging.warning("No clade data available for histograms")
        return
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Get unique clades
        clades = data['clade'].unique()
        clades = [c for c in clades if c != 'Unknown']  # Remove unknown
        
        if not clades:
            logging.warning("No valid clades found for histograms")
            return
        
        # Set up colors
        if clade_colors is None:
            clade_colors = plt.cm.Set3(np.linspace(0, 1, len(clades)))
        
        # Create subplots
        n_clades = len(clades)
        n_cols = min(3, n_clades)
        n_rows = (n_clades + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = [axes]
        elif n_rows == 1 or n_cols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        for i, clade in enumerate(clades):
            clade_data = data[data['clade'] == clade]
            
            if 'assembly' in clade_data.columns:
                # Count proteins per assembly for this clade
                proteins_per_assembly = clade_data.groupby('assembly').size()
                
                axes[i].hist(proteins_per_assembly, bins=20, alpha=0.7, 
                           color=clade_colors[i % len(clade_colors)])
                axes[i].set_title(f'Clade {clade} - Proteins per Assembly')
                axes[i].set_xlabel('Number of Proteins')
                axes[i].set_ylabel('Number of Assemblies')
        
        # Hide unused subplots
        for i in range(len(clades), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(output_dir, 'clade_histograms.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Saved clade histograms to: {output_file}")
        
    except Exception as e:
        logging.error(f"Failed to create clade histograms: {e}")


def plot_neighbor_types(cog_data: pd.DataFrame, output_dir: str = "output") -> None:
    """
    Plot the distribution of neighbor types.
    
    Args:
        cog_data: DataFrame with COG annotation data
        output_dir: Output directory for saving plots
    """
    logging.info("Plotting neighbor types")
    
    if cog_data is None or cog_data.empty:
        logging.warning("No COG data available for neighbor types plot")
        return
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Count neighbor types
        if 'COG_category' in cog_data.columns:
            type_counts = cog_data['COG_category'].value_counts()
            
            # Create plot
            fig, ax = plt.subplots(figsize=(12, 8))
            
            type_counts.plot(kind='bar', ax=ax)
            ax.set_title('Distribution of Neighbor Types (COG Categories)')
            ax.set_xlabel('COG Category')
            ax.set_ylabel('Count')
            ax.tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            
            # Save plot
            output_file = os.path.join(output_dir, 'types_of_neighbors.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logging.info(f"Saved neighbor types plot to: {output_file}")
        
    except Exception as e:
        logging.error(f"Failed to create neighbor types plot: {e}")


def plot_genomic_context(neighbor_data: pd.DataFrame, protein_id: str,
                        output_dir: str = "output") -> None:
    """
    Plot the genomic context around a protein of interest.
    
    Args:
        neighbor_data: DataFrame with neighbor information
        protein_id: Target protein ID
        output_dir: Output directory for saving plots
    """
    logging.info(f"Plotting genomic context for protein: {protein_id}")
    
    # Filter data for the protein of interest
    protein_neighbors = neighbor_data[neighbor_data['PIGI'] == protein_id]
    
    if protein_neighbors.empty:
        logging.warning(f"No neighbor data found for protein: {protein_id}")
        return
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Sort by genomic position
        if 'start' in protein_neighbors.columns:
            protein_neighbors = protein_neighbors.sort_values('start')
        
        # Create plot
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Plot each gene as a rectangle
        for _, gene in protein_neighbors.iterrows():
            start = gene.get('start', 0)
            end = gene.get('end', start + 1000)
            strand = gene.get('strand', '+')
            
            # Color based on whether it's the target protein or neighbor
            if gene.get('is_neighbour', True):
                color = 'lightblue'
            else:
                color = 'red'  # Target protein
            
            # Draw gene rectangle
            rect = patches.Rectangle((start, 0), end - start, 1, 
                                   linewidth=1, edgecolor='black', 
                                   facecolor=color, alpha=0.7)
            ax.add_patch(rect)
            
            # Add gene label
            ax.text((start + end) / 2, 0.5, gene.get('gene', ''), 
                   ha='center', va='center', fontsize=8, rotation=90)
            
            # Add arrow for strand direction
            if strand == '+':
                ax.arrow(end - 50, 0.8, 40, 0, head_width=0.1, 
                        head_length=30, fc='black', ec='black')
            else:
                ax.arrow(start + 50, 0.8, -40, 0, head_width=0.1, 
                        head_length=30, fc='black', ec='black')
        
        ax.set_xlim(protein_neighbors['start'].min() - 500, 
                   protein_neighbors['end'].max() + 500)
        ax.set_ylim(-0.2, 1.2)
        ax.set_xlabel('Genomic Position (bp)')
        ax.set_title(f'Genomic Context of Protein {protein_id}')
        ax.set_yticks([])
        
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(output_dir, f'genomic_context_{protein_id}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Saved genomic context plot to: {output_file}")
        
    except Exception as e:
        logging.error(f"Failed to create genomic context plot: {e}")