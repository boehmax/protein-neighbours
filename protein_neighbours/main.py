#!/usr/bin/env python3
"""
Main script for Protein Neighborhood Analysis.

This script coordinates the analysis of the genomic environment of proteins.
It loads configuration, sets up logging, and runs the analysis pipeline.
"""

import os
import sys
import logging
import click
from datetime import datetime
from typing import Dict, Any, Optional

# Add the package to Python path for development
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from protein_neighbours.utils import load_config, setup_logging, create_output_file_path
from protein_neighbours.io import (
    read_protein_assembly_data, read_clades, read_representatives, read_annotations, save_data
)
from protein_neighbours.neighbors import collect_all_neighbors, collect_all_protein_info
from protein_neighbours.annotation import analyze_proteins
from protein_neighbours.analysis import (
    analyze_neighbor_types, generate_summary_statistics, combine_and_analyze
)
from protein_neighbours.plotting import (
    plot_neighbors_per_clade, create_correlation_matrix, create_clade_histograms,
    plot_neighbor_types, plot_genomic_context
)


@click.command()
@click.option('--config', '-c', default="config/config.yaml", 
              help='Path to configuration YAML file')
@click.option('--interactive/--no-interactive', default=True,
              help='Whether to run in interactive mode')
@click.option('--protein-of-interest', '-p', default=None,
              help='Protein of interest to analyze')
@click.option('--output-dir', '-o', default=None,
              help='Override output directory')
def main(config: str, interactive: bool, protein_of_interest: Optional[str], 
         output_dir: Optional[str]) -> Dict[str, Any]:
    """
    Main function to run the protein neighborhood analysis.
    
    Args:
        config: Path to the configuration YAML file
        interactive: Whether to run in interactive mode
        protein_of_interest: Optional protein ID to analyze
        output_dir: Optional override for output directory
        
    Returns:
        Dictionary containing the analysis results
    """
    # Prepare override parameters
    override_params = {}
    if output_dir:
        override_params['paths'] = {'output_dir': output_dir}
    if protein_of_interest:
        override_params['analysis'] = {'protein_of_interest': protein_of_interest}
    
    # Load configuration
    try:
        config_dict = load_config(config, override_params)
    except Exception as e:
        print(f"Error loading configuration: {e}")
        sys.exit(1)
    
    # Set up logging
    try:
        setup_logging(config_dict)
        logging.info("Starting protein neighborhood analysis")
        logging.info(f"Using configuration file: {config}")
    except Exception as e:
        print(f"Error setting up logging: {e}")
        sys.exit(1)
    
    try:
        # Create output directory for the current date
        current_date = config_dict['analysis']['date']
        output_directory = os.path.join(config_dict['paths']['output_dir'], current_date)
        os.makedirs(output_directory, exist_ok=True)
        
        # Create output file paths
        output_file_neighbors = create_output_file_path(
            basepairs=config_dict['analysis']['basepairs'],
            max_neighbors=config_dict['analysis']['max_neighbors'],
            date=current_date,
            config=config_dict
        )
        output_file_proteins = os.path.join(output_directory, 'all_protein_info.csv')
        
        # Read protein and assembly data
        logging.info("Reading protein and assembly data")
        protein_data = read_protein_assembly_data(
            protein_file=config_dict['files']['proteins'],
            assembly_file=config_dict['files']['assemblies'],
            protein_assembly_file=config_dict['files']['protein_assembly'],
            base_path=config_dict['paths']['base_dir'],
            interactive=interactive,
            protein_of_interest=protein_of_interest
        )
        
        protein_assembly = protein_data['protein_assembly']
        target_protein = protein_data['protein_of_interest']
        
        # Generate protein alias data from representative files
        logging.info("Reading protein representatives")
        if config_dict['files'].get('representative_files'):
            rep_files = config_dict['files']['representative_files']
            protein_aliases = read_representatives(
                base_path=config_dict['paths']['base_dir'],
                rep_path=os.path.dirname(rep_files['ipg']),
                ipg_file=os.path.basename(rep_files['ipg']),
                pdb_file=os.path.basename(rep_files['pdb']),
                cluster_file=os.path.basename(rep_files['cluster'])
            )
        else:
            protein_aliases = None
        
        # Check if output files exist and read them if available
        if os.path.exists(output_file_neighbors):
            logging.info("Reading existing neighbor data from file")
            import pandas as pd
            all_neighbors = pd.read_csv(output_file_neighbors)
            
            if os.path.exists(output_file_proteins):
                logging.info("Reading existing protein data from file")
                all_proteins = pd.read_csv(output_file_proteins)
            else:
                logging.info("Collecting protein information")
                all_proteins = collect_all_protein_info(protein_assembly, 
                                                       config_dict['paths']['base_dir'])
        else:
            # Collect neighbor information
            logging.info("Collecting neighbor information")
            all_neighbors = collect_all_neighbors(
                protein_assembly,
                basepairs=config_dict['analysis']['basepairs'],
                max_neighbors=config_dict['analysis']['max_neighbors'],
                base_path=config_dict['paths']['base_dir'],
                overlap=config_dict['analysis']['overlap']
            )
            
            logging.info("Collecting protein information")
            all_proteins = collect_all_protein_info(protein_assembly,
                                                   config_dict['paths']['base_dir'])
            
            # Save the results
            logging.info("Saving neighbor and protein data")
            if not all_neighbors.empty:
                save_data(all_neighbors, os.path.basename(output_file_neighbors), 
                         os.path.dirname(output_file_neighbors))
            if not all_proteins.empty:
                save_data(all_proteins, os.path.basename(output_file_proteins),
                         os.path.dirname(output_file_proteins))
        
        # Combine all data
        import pandas as pd
        all_data = pd.concat([all_neighbors, all_proteins], ignore_index=True)
        
        # Plot genomic context if a valid protein of interest is provided
        if target_protein and not all_neighbors.empty:
            if protein_aliases is None or target_protein in protein_aliases['alias'].values:
                logging.info(f"Plotting genomic context for protein: {target_protein}")
                plot_genomic_context(all_neighbors, target_protein, output_directory)
            else:
                logging.warning('No valid protein of interest given, skipping genomic context plot.')
        
        # Get clade information
        logging.info("Reading clade information")
        clades = read_clades(
            base_path=config_dict['paths']['base_dir'],
            clade_dir=config_dict['files']['clade_dir'],
            pattern=config_dict['files']['clade_pattern']
        )
        
        # Annotate proteins
        logging.info(f"Annotating proteins with {config_dict['annotation']['tool']}")
        annotation_results = analyze_proteins(
            df=all_neighbors,
            column='gene',
            config=config_dict
        )
        
        # Analyze neighbor types
        logging.info("Analyzing neighbor types")
        analyze_neighbor_types(annotation_results, output_directory)
        
        # Read manual annotations if available
        logging.info("Reading any manual annotations")
        annotated_neighbors = read_annotations(current_date, config_dict['paths']['output_dir'])
        
        # Combine and analyze data
        logging.info("Combining data and generating analysis")
        combined_df = combine_and_analyze(
            neighbors_data=all_data,
            cog_data=annotation_results,
            clade_assign=clades,
            neighbor_annotations=annotated_neighbors
        )
        
        # Generate summary statistics
        logging.info("Generating summary statistics")
        summary_stats = generate_summary_statistics(combined_df, output_directory)
        
        # Generate visualizations
        logging.info("Generating visualization plots")
        
        # Standard plots
        plot_neighbors_per_clade(
            combined_df,
            exclude_unknown_clade=config_dict['visualization']['exclude_unknown_clade'],
            exclude_unknown_cog=config_dict['visualization']['exclude_unknown_cog']
        )
        
        plot_neighbors_per_clade(combined_df)
        
        # Plots with CODH count
        plot_neighbors_per_clade(
            combined_df,
            exclude_unknown_clade=config_dict['visualization']['exclude_unknown_clade'],
            exclude_unknown_cog=config_dict['visualization']['exclude_unknown_cog'],
            plot_count_codh=True
        )
        
        plot_neighbors_per_clade(combined_df, plot_count_codh=True)
        
        # Correlation matrix
        if 'clade' in combined_df.columns and 'assembly' in combined_df.columns:
            logging.info("Generating correlation matrix")
            correlation_data = combined_df[['assembly', 'clade']].drop_duplicates()
            create_correlation_matrix(
                correlation_data,
                combined_df['clade'].unique().tolist(),
                output_directory
            )
        
        # Clade histograms
        if 'clade' in combined_df.columns:
            logging.info("Generating clade histograms")
            protein_data_for_hist = combined_df[['PIGI', 'assembly', 'clade']].drop_duplicates()
            create_clade_histograms(
                protein_data_for_hist,
                output_directory,
                config_dict['visualization'].get('clade_colors')
            )
        
        # Neighbor type plot
        if not annotation_results.empty:
            plot_neighbor_types(annotation_results, output_directory)
        
        # Generate plots with manual annotations if available
        if annotated_neighbors is not None and not annotated_neighbors.empty:
            logging.info("Generating plots with manual annotations")
            # Add manual annotation column to combined data
            combined_annotated = combined_df.copy()
            if 'ANNOTATION' in annotated_neighbors.columns:
                combined_annotated['COG_LETTER'] = combined_annotated.get('ANNOTATION', 'Unknown')
            
            plot_neighbors_per_clade(
                combined_annotated,
                exclude_unknown_clade=True,
                exclude_unknown_cog=True,
                output_path="annotated_neighbors"
            )
        
        # Generate HTML report
        logging.info("Analysis pipeline completed successfully")
        generate_html_report(combined_df, config_dict, output_directory)
        
        logging.info(f"Analysis complete. Results available in {output_directory}")
        
        # Return results
        return {
            'combined_df': combined_df,
            'all_neighbors': all_neighbors,
            'all_proteins': all_proteins,
            'clades': clades,
            'annotation_results': annotation_results,
            'summary_stats': summary_stats,
            'config': config_dict
        }
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        raise


def generate_html_report(combined_df, config: Dict[str, Any], output_dir: str) -> Optional[str]:
    """
    Generate an HTML report of the analysis.
    
    Args:
        combined_df: Combined analysis DataFrame
        config: Configuration dictionary
        output_dir: Output directory
        
    Returns:
        Path to generated report or None if failed
    """
    try:
        from jinja2 import Template
        
        # Create a simple HTML report template
        template_str = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Protein Neighborhood Analysis Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                h1, h2 { color: #333; }
                table { border-collapse: collapse; width: 100%; }
                th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                th { background-color: #f2f2f2; }
                .summary { background-color: #f9f9f9; padding: 20px; margin: 20px 0; }
            </style>
        </head>
        <body>
            <h1>Protein Neighborhood Analysis Report</h1>
            <p>Generated on: {{ date }}</p>
            
            <div class="summary">
                <h2>Summary Statistics</h2>
                <p><strong>Total Proteins:</strong> {{ total_proteins }}</p>
                <p><strong>Total Neighbors:</strong> {{ total_neighbors }}</p>
                <p><strong>Number of Clades:</strong> {{ num_clades }}</p>
                <p><strong>Number of Assemblies:</strong> {{ num_assemblies }}</p>
            </div>
            
            <h2>Configuration Used</h2>
            <pre>{{ config_yaml }}</pre>
            
            <h2>Output Files</h2>
            <ul>
                <li>neighbor_distribution_by_clade.png - Neighbor distribution visualization</li>
                <li>correlation_matrix.png - Clade co-occurrence correlation</li>
                <li>clade_histograms.png - Protein distribution per clade</li>
                <li>summary_statistics.csv - Summary statistics</li>
                <li>types_of_neighbors.csv - Neighbor type analysis</li>
            </ul>
        </body>
        </html>
        """
        
        template = Template(template_str)
        
        # Prepare data for template
        total_proteins = len(combined_df[~combined_df['is_neighbour']]) if 'is_neighbour' in combined_df.columns else len(combined_df)
        total_neighbors = len(combined_df[combined_df['is_neighbour']]) if 'is_neighbour' in combined_df.columns else 0
        num_clades = len(combined_df['clade'].unique()) if 'clade' in combined_df.columns else 0
        num_assemblies = len(combined_df['assembly'].unique()) if 'assembly' in combined_df.columns else 0
        
        import yaml
        config_yaml = yaml.dump(config, default_flow_style=False)
        
        # Render template
        html_content = template.render(
            date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            total_proteins=total_proteins,
            total_neighbors=total_neighbors,
            num_clades=num_clades,
            num_assemblies=num_assemblies,
            config_yaml=config_yaml
        )
        
        # Save report
        report_file = os.path.join(output_dir, "analysis_report.html")
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        logging.info(f"Analysis report generated: {report_file}")
        return report_file
        
    except Exception as e:
        logging.error(f"Failed to generate HTML report: {e}")
        return None


if __name__ == "__main__":
    main()