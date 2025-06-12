# Protein Genomic Environment Analysis

This project provides a comprehensive **R package and pipeline** for analyzing the genomic neighborhood of proteins. You can use it as a standard R package (recommended for reproducibility and modularity), or simply run the main script directly if you prefer a quick, script-based workflow.

---

## Installation

### Option 1: Install as an R Package (Recommended)

1. **Clone the repository:**
    ```sh
    git clone https://github.com/boehmax/protein-neighbours.git
    cd protein-neighbours
    ```

2. **Open R (or RStudio) in the project directory and run:**
    ```r
    # If you don't have devtools, install it first:
    # install.packages("devtools")
    devtools::install()
    ```

3. **Load the package:**
    ```r
    library(proteinNeighbours)
    ```

4. **Run the pipeline:**
    ```r
    results <- main()
    ```

---

### Option 2: Run as a Standalone Script

If you don’t want to install the package, you can simply source and run the main script:

1. **Clone the repository and set your working directory to the project folder.**

2. **In R:**
    ```r
    source("main.R")
    results <- main()
    ```

---

### Prerequisites

- **R** (>= 4.0.0)
- Required R packages (will be installed automatically if you use the package):
    - tidyverse
    - ape
    - ggplot2
    - gggenes
    - RColorBrewer
    - yaml
    - logger (optional, for enhanced logging)
    - rmarkdown (optional, for report generation)
- **eggNOG-mapper** (external, for protein annotation)
    - [eggNOG-mapper Documentation](http://eggnog-mapper.embl.de/)
    - Install with: `pip install eggnog-mapper`

---

You can now follow the rest of the instructions below for configuration, data requirements, and usage.  
**Whether you use the package or just the script, the workflow and outputs are the same!**

## Overview

This package explores and analyzes the genomic environment of proteins by:

1. Identifying neighboring proteins within a specified distance
2. Annotating proteins using eggNOG-mapper with COG classification
3. Grouping proteins by clades or other classifications
4. Generating visualizations to understand distribution patterns
5. Creating reports for sharing and documentation

## Running the Pipeline

1. **Source the main script:**

    ```r
    source("main.R")
    ```

2. **Run the analysis with default configuration:**

    ```r
    results <- main()
    ```

3. **Run with a custom configuration file:**

    ```r
    results <- main("path/to/your/config.yaml")
    ```

4. **Override specific parameters directly:**

    ```r
    results <- main(override_params = list(
      "analysis.basepairs" = 500,
      "analysis.max_neighbors" = 20,
      "annotation.tool" = "eggnog"
    ))
    ```

---

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

### Required CSV Files

Place the following files in the directory specified by `paths.base_dir` (default: `data/`):

- `proteins.csv`: List of protein accession numbers
- `assm_accs.csv`: List of assembly accession numbers
- `assm_accs_protein.csv`: Mapping between proteins and assemblies

### Required Genome Files

The package expects GFF3 files for each assembly in:

```
{paths.base_dir}/ncbi_dataset/data/{assembly_id}/genomic.gff
```

You can download these files using NCBI Datasets:

```bash
datasets download genome accession --inputfile assm_accs.csv --include gff3
```
or check out my other repository: https://github.com/boehmax/protein-to-genome which I built this pipline on.

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

These files are tables that mapp accession numbers with other identifiers. 
- 'ipg_representative.txt' you will get from boehmax/protein-to-genome, and makes sure that accessions from identical proteins are connected via one overarchning accesion number
- 'pdb_representative.txt' you would need to wirte manually, to make sure that pdb ids are mapped to NCBI accession numbers
- 'cluster_representative.txt' could be generated after running CD-hit on a protein file, to make sure that accession numbers from a cluster all are connected via one overarching accesion number (a script to generate this file is part of boehmax/protein-per-organism, https://github.com/boehmax/protein-per-organism/blob/main/cluster_alias.sh)

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
5. If done on different date, change folder name to appropriate date.

### Creating Custom Plots

After a succesfull run the package provides several plotting functions you can use directly:

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
