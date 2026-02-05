# Targets Pipeline Guide

## Introduction

This project uses the [targets](https://docs.ropensci.org/targets/) R package for pipeline management. The targets framework provides:

- **Automatic caching**: Only re-runs steps when their inputs change
- **Dependency tracking**: Automatically determines what needs updating
- **Parallel execution**: Runs independent steps simultaneously
- **Reproducibility**: Tracks all dependencies and versions
- **Visualization**: Interactive pipeline graphs

## Quick Start

### Running the Pipeline

```r
# Load targets
library(targets)

# Run the entire pipeline
tar_make()

# View results
tar_read(combined_df)
```

### Visualizing the Pipeline

```r
# Interactive network diagram
tar_visnetwork()

# Static graph
tar_network()
```

## Pipeline Structure

The pipeline is defined in `_targets.R` and consists of 21 targets:

1. **config**: Load configuration from YAML
2. **logging_setup**: Initialize logging system
3. **output_dir**: Create output directory
4. **output_files**: Define file paths
5. **protein_assembly_data**: Read input data
6. **protein_alias**: Read protein representatives
7. **all_neighbours**: Collect neighbor information (cached)
8. **all_protein**: Collect protein information (cached)
9. **all_data**: Combine neighbor and protein data
10. **neighbor_plot**: Plot neighbors for protein of interest
11. **clades**: Read clade information
12. **annotation_results**: Annotate proteins with eggNOG/COG
13. **neighbor_types**: Analyze neighbor types
14. **annotated_neighbours**: Read manual annotations
15. **combined_df**: Combine and plot all data
16. **standard_plots**: Generate standard visualizations
17. **correlation_matrix**: Create correlation matrix
18. **clade_histograms**: Generate clade histograms
19. **annotated_plots**: Generate plots with manual annotations
20. **analysis_report**: Create HTML report
21. **analysis_summary**: Final summary and outputs

## Common Tasks

### Check Pipeline Status

```r
# See what's up to date
tar_progress()

# Check which targets are outdated
tar_outdated()

# View pipeline metadata
tar_meta()
```

### Load Results

```r
# Load specific targets
config <- tar_read(config)
neighbours <- tar_read(all_neighbours)
combined_df <- tar_read(combined_df)

# Or use the convenience function
library(proteinNeighbours)
df <- load_target("combined_df")
```

### Force Re-execution

```r
# Invalidate specific targets
tar_invalidate(c("annotation_results"))

# Then re-run
tar_make()

# Or clean everything and start fresh
tar_destroy()
tar_make()
```

## Parallel Execution

### Local Parallel Execution

```r
# Using future package
library(future)
plan(multisession, workers = 4)
tar_make_future()

# Or use the convenience function
library(proteinNeighbours)
run_targets_parallel(workers = 4)
```

## Debugging

### Interactive Debugging

```r
# Load the workspace for a specific target
tar_workspace(annotation_results)

# Load all dependencies of a target
tar_load(all_neighbours)
tar_load(config)

# Now you can run the target code interactively
```

## Additional Resources

- [Targets Manual](https://books.ropensci.org/targets/)
- [Targets Reference](https://docs.ropensci.org/targets/)
- [GitHub Issues](https://github.com/ropensci/targets/issues)
