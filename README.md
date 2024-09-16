# Genomic Environment Analysis

**Work in Progress**

This project is currently under active development and is not yet complete. Features and functionality may change, and the documentation may not be fully up to date.

This script is designed to explore and analyze the genomic environment of a given protein. It performs the following tasks:

1. Loads necessary libraries.
2. Sources other scripts.
3. Reads in data.
4. Performs analysis.
5. Generates plots.

## Libraries Used
- `tidyverse`
- `ape`
- `dplyr`
- `gggenes`
- `ggplot2`
- `RColorBrewer`

## Scripts Sourced
- `01_open.R`
- `02_clean.R`
- `03_functions.R`

## Main Function
The main function of the script performs the following steps:

1. Reads protein and assembly data.
2. Generates protein alias data.
3. Collects neighboring proteins.
4. Checks if results are saved.
5. Plots the neighbors.
6. Reads clade information.
7. Reads cluster domain information.
8. Gets the types of neighbors and their counts.
9. Re-imports the types of neighbors after manual annotation.
10. Combines and plots the data.


## Files to Supply

Ensure the following files are placed in the directory specified by the `PATH` variable:

- `proteins.csv`: Contains a list of all proteins of interest for neighbor analysis.
- `assm_accs.csv`: Contains genome assembly accession numbers, corresponding to folder names in the `PATH/ncbi_dataset/data/` directory.
- `assm_accs_protein.csv`: A CSV file with genome assembly accession numbers in the first column and protein accession numbers in the second column.
- `ncbi_dataset/data/GENOME_ACCESSION/genomic.gff`: Contains downloaded and extracted genomic information from the NCBI database in `.gff` format. This file can be obtained using the following command:
    ```bash
    PATH=PATH/TO/YOUR/ncbi/datasets
    $PATH download genome accession --inputfile assm_accs.csv --include gff3
    ```
## Optional Files

### Clade Files
- **`PATH/clades/Clade*.txt`**: Create a `clades` folder where you save lists of proteins. Each file should correspond to a specific clade (or any other grouping you want to analyze).

### Representative Files
- **`PATH/representatives/ipg_representative.csv`**
- **`PATH/representatives/pdb_representative.csv`**
- **`PATH/representatives/cluster_representative.csv`**

During the analysis of the primary list of proteins, it might happen that the protein you provided as input (PIGI...protein I gave as input) no longer corresponds to the protein accession numbers in the genomic data you downloaded. It is therefore important to collect files to cross-reference how the numbers changed to avoid losing too much data. All these `.csv` files should have the following format:

- **`accession`**: The current name of the protein through the analysis.
- **`IPGI`**: The protein accession you provided as input, which resulted in the `accession`.
- **`identity`**: Specifies how the amino acid sequence of the PIGI relates to that of the `accession`. For the `ipg_` and `pdb_` files, this can be omitted since it is expected to be 1. For the `cluster` files, this can be anything between 0.9 to 1 (depending on what you specified in your CD-Hit run as the cutoff).

### COG files
- PATH/classification/COG/cog_*hitdata.txt: File(s) created by XXX webserver to generate annotation of genes.

## Parameters

- **BASEPAIRS**: The maximum distance in base pairs between two proteins to consider them as neighbors. Default is `300`.
- **MAX_NEIGHBORS**: The maximum number of neighboring proteins to consider. Default is `15`.
- **protein_of_interest**: The specific protein to focus on for closer analysis. This parameter is optional.
- **PATH**: The directory path to your input files. Default is `'data'`.


## Outputs
The outputs are saved in a folder named with the current date. If you want to run multiple analyses, please rename the previously created folder before running the script again to avoid overwriting the output files. 

The script generates multiple CSV files containing information about the neighboring proteins and other relevant data. It also produces plots showing the number of proteins per type per clade. 

For more concise plots, you can specify short names for proteins in the `output/YYYY-MM-DD/types_of_neighbours_annotated.csv` file. The program will initially export a `output/YYYY-MM-DD/types_of_neighbours.csv` file. You can then annotate this file in the third column with names or abbreviations as needed and save it as `types_of_neighbours_annotated.csv`.
 

## Usage

1. Ensure all required data files are available and placed in the specified directories and all required R libraries are installed.
2. Open `00_main.R` and modify the `main()` function at the bottom of the script to include your desired parameters.
3. Run the script. After execution, a `types_of_neighbours.csv` file will be generated.
4. (Optional) Annotate the `types_of_neighbours.csv` file as needed and rerun the script. The rerun will be faster since the neighboring proteins have already been calculated and will be reloaded.

Happy analyzing!


