# Targets Quick Reference

## Essential Commands

### Running the Pipeline
```r
library(targets)
tar_make()                              # Run the entire pipeline
tar_make(names = c("config", "data"))   # Run specific targets
tar_make_future(workers = 4)            # Run with parallel execution
```

### Visualizing
```r
tar_visnetwork()      # Interactive dependency graph
tar_glimpse()         # Text-based dependency view
tar_manifest()        # List all targets
```

### Checking Status
```r
tar_outdated()        # Which targets need updating
tar_progress()        # Execution status
tar_meta()            # Detailed metadata
tar_sitrep()          # Situation report
```

### Loading Results
```r
tar_read(combined_df)         # Load a specific target
tar_load(all_neighbours)      # Load into workspace
tar_load_everything()         # Load all targets
```

### Cache Management
```r
tar_invalidate(c("target1", "target2"))  # Mark targets as outdated
tar_delete(c("target1", "target2"))      # Delete specific targets
tar_prune()                              # Remove obsolete targets
tar_destroy()                            # Delete entire cache (use with caution!)
```

### Debugging
```r
tar_workspace(target_name)    # Load workspace for debugging
tar_traceback(target_name)    # View error traceback
```

## Convenience Functions (proteinNeighbours package)

```r
library(proteinNeighbours)

run_targets_pipeline()              # Run pipeline with nice defaults
visualize_pipeline()                # Visualize dependencies
load_target("combined_df")          # Load target with error handling
get_pipeline_progress()             # Check progress
clean_targets_cache()               # Clean cache with confirmation
run_targets_parallel(workers = 4)   # Run with parallel workers
```

## Common Workflows

### First Time Run
```r
tar_make()              # Run everything
tar_visnetwork()        # See what was built
tar_read(combined_df)   # Load results
```

### Incremental Update
```r
# Edit config/config.yaml
tar_outdated()          # Check what will re-run
tar_make()              # Run updated targets only
```

### After Code Changes
```r
tar_invalidate(everything())  # Mark all outdated
tar_make()                    # Re-run everything
```

### Selective Re-run
```r
tar_invalidate(c("annotation_results"))  # Mark specific targets
tar_make()                               # Re-run affected targets
```

### Debugging Issues
```r
tar_meta(fields = error, complete_only = TRUE)  # Check errors
tar_workspace(problematic_target)               # Load workspace
tar_load_everything()                           # Load all dependencies
# Debug interactively...
```

## Tips

- Run `tar_visnetwork()` early to understand dependencies
- Use `tar_outdated()` before `tar_make()` to preview changes
- Keep `_targets/` in `.gitignore` - it can get large
- Use `tar_prune()` regularly to clean up old targets
- For long-running pipelines, use `reporter = "verbose"` to see progress
- Store large results with `format = "qs"` for speed

## Common Options in _targets.R

```r
tar_option_set(
  packages = c("dplyr", "ggplot2"),     # Required packages
  format = "rds",                        # Storage format (rds, qs, feather, etc.)
  error = "continue",                    # Continue on errors
  memory = "transient",                  # Unload from memory ASAP
  garbage_collection = TRUE              # Run gc() between targets
)
```

## Environment Variables

```bash
export TAR_PROJECT = "my_project"      # Use named project
export TAR_ASK = "false"               # Skip confirmations
```

## Getting Help

```r
?tar_make          # Help for specific function
?targets::         # Browse all targets functions
tar_config_get()   # View current configuration
```
