# Protein Genomic Environment Analysis

A comprehensive R package for analyzing the genomic neighborhood of proteins. This package examines neighboring proteins, annotates them using eggNOG-mapper, and visualizes the results to provide insights into genomic context and potential functional relationships.

## Overview

This package explores and analyzes the genomic environment of proteins by:

1. Identifying neighboring proteins within a specified distance
2. Annotating proteins using eggNOG-mapper with COG classification
3. Grouping proteins by clades or other classifications
4. Generating visualizations to understand distribution patterns
5. Creating reports for sharing and documentation

## Installation

### Prerequisites

- R (>= 4.0.0)
- Required R packages (will be installed automatically):
  - tidyverse
  - ape
  - ggplot2
  - gggenes
  - RColorBrewer
  - yaml
  - logger (optional, for enhanced logging)
  - rmarkdown (optional, for report generation)

### External Dependencies

- **eggNOG-mapper**: For protein annotation (replaces the previously used COGclassifier)
  - Installation: [eggNOG-mapper Documentation](http://eggnog-mapper.embl.de/)
  - `pip install eggnog-mapper`

### Installing the Package

```r
# Install from GitHub
devtools::install_github("yourusername/protein-neighbours")

# Or install from local directory
install.packages("path/to/protein-neighbours", repos = NULL, type = "source")
```

## Quick Start

```r
# Load the package
library(proteinNeighbours)

# Run the analysis with default configuration
results <- main()

# Run with custom configuration file
results <- main("path/to/your/config.yaml")

# Override specific parameters
results <- main(override_params = list(
  "analysis.basepairs" = 500,
  "analysis.max_neighbors" = 20,
  "annotation.tool" = "eggnog"
))
```

## Configuration

The package uses a YAML configuration file for all parameters. You can customize the analysis by:

1. Editing the default config file at `config/config.yaml`
2. Creating your own config file and specifying it with `main("path/to/config.yaml")`
3. Overriding specific parameters with the `override_params` argument

### Key Configuration Parameters

```yaml
# Analysis parameters
analysis:
  basepairs: 300                     # Distance in base pairs to consider for neighbors
  max_neighbors: 15                  # Maximum number of neighbors to identify
  protein_of_interest: null          # Default to null, can be set at runtime

# Annotation settings
annotation:
  tool: "eggnog"                     # Tool to use for annotation: "eggnog" or "cog"
  eggnog:
    db_dir: "data/eggnog_db"         # Path to eggNOG database
    cpu: 4                           # Number of CPUs to use
```

See the full `config.yaml` file for all available options.

## Data Requirements

### Required Files

Place the following files in the directory specified by `paths.base_dir` (default: `data/`):

- `proteins.csv`: List of protein accession numbers
- `assm_accs.csv`: List of assembly accession numbers
- `assm_accs_protein.csv`: Mapping between proteins and assemblies

### GFF Files

The package expects GFF3 files for each assembly in:

```
{paths.base_dir}/ncbi_dataset/data/{assembly_id}/genomic.gff
```

You can download these files using NCBI Datasets:

```bash
datasets download genome accession --inputfile assm_accs.csv --include gff3
```

### Optional Files

#### Clade Files

Create a `clades` folder with files for each clade:

```
{paths.base_dir}/clades/Clade*.txt
```

#### Representative Files

For protein alias mapping:

```
{paths.base_dir}/representatives/ipg_representative.txt
{paths.base_dir}/representatives/pdb_representative.txt
{paths.base_dir}/representatives/cluster_representative.txt
```

## Outputs

The package generates outputs in a directory named with the current date (or specified date):

```
output/{date}/
├── all_neighbours_bp{basepairs}_n{max_neighbors}.csv   # Neighbor information
├── all_protein_info.csv                               # Protein information
├── analysis_config.yaml                               # Configuration used
├── input_file_info.csv                                # Input file metadata
├── analysis_report.html                               # HTML report
├── eggnog/                                            # eggNOG annotation results
├── types_of_neighbours.csv                            # Neighbor type counts
└── various plot files (PNG, SVG)                      # Visualizations
```

### Visualization Examples

- Distribution of neighbors by clade
- Correlation matrix of clade co-occurrence
- Histograms of protein distribution
- Annotation type distribution

## Advanced Usage

### Manual Annotation

1. Run the initial analysis
2. Examine `output/{date}/types_of_neighbours.csv`
3. Create an annotated version named `types_of_neighbours_annotated.csv` with your classification in the 4th column
4. Re-run the analysis to incorporate your annotations

### Creating Custom Plots

The package provides several plotting functions you can use directly:

```r
# Plot neighbors for a specific protein
plot_neighbours(all_neighbours, "YP_123456.1")

# Plot neighbors per clade with custom options
plot_neighbours_per_clade(combined_df, exclude_unknown_clade = TRUE, plot_count_codh = TRUE)

# Create a correlation matrix
make_correlation_matrix(df, clade_vector)

# Create histograms for clades
create_clade_histograms2(df)
```

## Troubleshooting

### Common Issues

- **Missing GFF files**: Ensure assembly GFF files are in the correct location
- **Protein not found**: Check protein accession numbers and use representative mapping
- **eggNOG-mapper errors**: Verify eggNOG installation and database files
- **Memory issues**: Reduce the number of proteins or assemblies being analyzed

### Logging

The package uses the logger package for detailed logging. To enable debugging:

```r
# In your config.yaml
logging:
  level: "DEBUG"  # Options: DEBUG, INFO, WARNING, ERROR
  file: "protein_neighbors.log"
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Original development by Maximilian Böhm
- eggNOG-mapper for protein annotation
- NCBI for genome data and protein information

---

For more detailed documentation, see the vignettes and function documentation within the package.
