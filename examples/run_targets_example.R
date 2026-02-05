#!/usr/bin/env Rscript
#' Example Script: Running the Targets Pipeline
#'
#' This script demonstrates how to use the targets pipeline for protein
#' neighborhood analysis. It can be used as a template for your own analyses.
#'
#' Usage:
#'   Rscript examples/run_targets_example.R
#'
#' Or source it in R:
#'   source("examples/run_targets_example.R")

# Load required libraries
library(targets)

cat("=== Protein Neighborhood Analysis with Targets ===\n\n")

# Check if we're in the right directory
if (!file.exists("_targets.R")) {
  stop("Error: _targets.R not found. Please run this script from the project root directory.")
}

# 1. Visualize the pipeline structure
cat("Step 1: Visualizing pipeline structure...\n")
cat("   Run tar_visnetwork() in an interactive R session to see the dependency graph.\n\n")

# 2. Check which targets are outdated (initially, all will be outdated)
cat("Step 2: Checking which targets need to be run...\n")
outdated <- tar_outdated()
cat("   Outdated targets:", length(outdated), "\n")
if (length(outdated) > 0) {
  cat("   ", paste(head(outdated, 10), collapse = ", "), "\n")
  if (length(outdated) > 10) {
    cat("   ... and", length(outdated) - 10, "more\n")
  }
}
cat("\n")

# 3. Run the pipeline
cat("Step 3: Running the pipeline...\n")
cat("   This may take some time depending on your data size.\n")
cat("   Note: This example script doesn't actually run tar_make() to avoid\n")
cat("   requiring data files. To run the pipeline, uncomment the line below:\n\n")
cat("   tar_make()\n\n")

# Uncomment the following line to actually run the pipeline:
# tar_make()

# 4. Show how to check progress
cat("Step 4: Checking pipeline progress...\n")
cat("   After running tar_make(), use:\n")
cat("   - tar_progress() to see execution status\n")
cat("   - tar_meta() to see detailed metadata\n\n")

# 5. Show how to load results
cat("Step 5: Loading results...\n")
cat("   After successful execution, load results with:\n")
cat("   - config <- tar_read(config)\n")
cat("   - combined_df <- tar_read(combined_df)\n")
cat("   - all_neighbours <- tar_read(all_neighbours)\n")
cat("   - results <- tar_read(analysis_summary)\n\n")

# 6. Show convenience functions
cat("Step 6: Using convenience functions...\n")
cat("   The package provides helper functions:\n")
cat("   - run_targets_pipeline() - Run the pipeline\n")
cat("   - visualize_pipeline() - Visualize dependencies\n")
cat("   - load_target('target_name') - Load specific results\n")
cat("   - get_pipeline_progress() - Check progress\n")
cat("   - clean_targets_cache() - Reset the cache\n\n")

# 7. Example of selective execution
cat("Step 7: Selective execution...\n")
cat("   To run only specific targets:\n")
cat("   tar_make(names = c('config', 'protein_assembly_data'))\n\n")

# 8. Example of parallel execution
cat("Step 8: Parallel execution...\n")
cat("   To run with parallel workers:\n")
cat("   library(future)\n")
cat("   plan(multisession, workers = 4)\n")
cat("   tar_make_future()\n\n")
cat("   Or use the convenience function:\n")
cat("   library(proteinNeighbours)\n")
cat("   run_targets_parallel(workers = 4)\n\n")

cat("=== Example Complete ===\n")
cat("\nFor more information, see:\n")
cat("  - README.md for quick start guide\n")
cat("  - docs/TARGETS_GUIDE.md for detailed documentation\n")
cat("  - Run ?targets::tar_make for targets package help\n")
