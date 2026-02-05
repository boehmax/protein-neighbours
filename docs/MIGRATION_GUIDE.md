# Migration Guide: Moving to Targets Pipeline

This guide helps existing users of the protein-neighbours pipeline migrate to the new targets-based workflow.

## Do I Need to Migrate?

**No!** The traditional `main()` function still works exactly as before. However, migrating to targets provides significant benefits:

- ⚡ **Faster**: Only re-runs changed steps
- 🔍 **Transparent**: See exactly what's happening
- 🔄 **Reproducible**: Complete dependency tracking
- ⚙️ **Parallel**: Run independent steps simultaneously
- 🐛 **Debuggable**: Inspect intermediate results anytime

## Quick Migration (5 minutes)

If you're currently using:
```r
source("main.R")
results <- main()
```

Switch to:
```r
library(targets)
tar_make()
results <- tar_read(analysis_summary)
```

That's it! Everything else works the same.

## Detailed Migration Steps

### Step 1: Understanding Your Current Workflow

Your current workflow probably looks like this:

```r
# Traditional workflow
source("main.R")
results <- main()

# Or with custom config
results <- main("config/my_config.yaml")

# Or with parameter overrides
results <- main(override_params = list(
  "analysis.basepairs" = 500
))
```

### Step 2: Understanding the New Workflow

The targets workflow separates configuration from execution:

```r
# New workflow
library(targets)

# 1. Configure (edit config/config.yaml)
# 2. Run
tar_make()

# 3. Load results
results <- tar_read(analysis_summary)
combined_df <- tar_read(combined_df)
```

### Step 3: Migrating Your Configuration

Your existing `config/config.yaml` works without changes! The targets pipeline uses the same configuration system.

If you had multiple config files:
```
config/
├── config.yaml              # Default
├── experiment1.yaml         # Custom config 1
└── experiment2.yaml         # Custom config 2
```

You have two options:

**Option A: Switch configs manually**
```r
# Copy your config
file.copy("config/experiment1.yaml", "config/config.yaml", overwrite = TRUE)
tar_make()
```

**Option B: Modify _targets.R (advanced)**
```r
# In _targets.R, change:
tar_target(
  config,
  load_config(Sys.getenv("CONFIG_FILE", "config/config.yaml"))
)

# Then use:
Sys.setenv(CONFIG_FILE = "config/experiment1.yaml")
tar_make()
```

### Step 4: Migrating Your Analysis Scripts

If you have custom scripts that call `main()`:

**Before:**
```r
# my_analysis.R
source("main.R")
results <- main()

# Do custom analysis
my_plot <- custom_analysis(results$combined_df)
```

**After:**
```r
# my_analysis.R
library(targets)
tar_make()

# Load what you need
combined_df <- tar_read(combined_df)

# Do custom analysis
my_plot <- custom_analysis(combined_df)
```

### Step 5: Migrating Batch Processing

If you process multiple datasets:

**Before:**
```r
configs <- c("config1.yaml", "config2.yaml", "config3.yaml")
results <- list()

for (cfg in configs) {
  results[[cfg]] <- main(cfg)
}
```

**After:**
```r
configs <- c("config1.yaml", "config2.yaml", "config3.yaml")

for (cfg in configs) {
  file.copy(cfg, "config/config.yaml", overwrite = TRUE)
  tar_make()
  
  # Save results
  saveRDS(tar_read(analysis_summary), 
          paste0("results_", basename(cfg), ".rds"))
  
  # Clean for next run
  tar_destroy()
}
```

### Step 6: Understanding Caching

The targets pipeline caches intermediate results:

```r
# First run: Everything executes
tar_make()  # Takes 20 minutes

# No changes: Nothing executes
tar_make()  # Takes 2 seconds

# Config change: Only affected targets execute
# Edit config/config.yaml, change basepairs: 300 -> 500
tar_make()  # Takes 10 minutes (only downstream work)
```

To force a complete re-run:
```r
tar_destroy()  # Delete cache
tar_make()     # Run everything
```

### Step 7: Migrating Custom Functions

If you source additional custom functions:

**Before:**
```r
source("my_custom_functions.R")
source("main.R")
results <- main()
```

**After:**
Add to `_targets.R` near the top:
```r
# After the existing source() calls:
source("my_custom_functions.R")
```

Then run normally:
```r
tar_make()
```

### Step 8: Accessing Intermediate Results

One major benefit of targets is easy access to intermediate results.

**Before:** You couldn't easily access intermediate results
```r
results <- main()
# All you get is the final results
```

**After:** Load any intermediate step
```r
tar_make()

# Load any target
config <- tar_read(config)
neighbours <- tar_read(all_neighbours)
proteins <- tar_read(all_protein)
clades <- tar_read(clades)
annotations <- tar_read(annotation_results)
combined <- tar_read(combined_df)
```

## Common Migration Scenarios

### Scenario 1: Interactive Analysis

**Before:**
```r
source("main.R")
results <- main()
head(results$combined_df)
```

**After:**
```r
tar_make()
combined_df <- tar_read(combined_df)
head(combined_df)
```

### Scenario 2: Reproducible Research

**Before:**
```r
# run_analysis.R
source("main.R")
results <- main("config/paper_config.yaml")
saveRDS(results, "results_for_paper.rds")
```

**After:**
```r
# run_analysis.R
library(targets)
file.copy("config/paper_config.yaml", "config/config.yaml")
tar_make()
saveRDS(tar_read(analysis_summary), "results_for_paper.rds")
```

### Scenario 3: Parameter Exploration

**Before:**
```r
basepairs <- c(100, 200, 300, 400, 500)
results <- list()

for (bp in basepairs) {
  results[[as.character(bp)]] <- main(
    override_params = list("analysis.basepairs" = bp)
  )
}
```

**After:**
```r
library(yaml)
basepairs <- c(100, 200, 300, 400, 500)

for (bp in basepairs) {
  # Update config
  config <- read_yaml("config/config.yaml")
  config$analysis$basepairs <- bp
  write_yaml(config, "config/config.yaml")
  
  # Run pipeline
  tar_make()
  
  # Save results
  saveRDS(tar_read(analysis_summary), 
          sprintf("results_bp%d.rds", bp))
  
  # Clean for next iteration
  tar_destroy()
}
```

### Scenario 4: Debugging

**Before:** Hard to debug intermediate steps
```r
# Had to add print statements and re-run
```

**After:** Load and inspect any step
```r
tar_make()  # Run until error

# Load everything up to the error
tar_load_everything()

# Now you can debug interactively
ls()  # See all loaded objects
debug(problematic_function)
```

## Tips for Successful Migration

### 1. Start Fresh
```r
# Clean start
tar_destroy()
tar_make()
```

### 2. Visualize Dependencies
```r
# See what depends on what
tar_visnetwork()
```

### 3. Check What Will Run
```r
# Before running, check what's outdated
tar_outdated()
```

### 4. Monitor Progress
```r
# In another R session while pipeline runs:
tar_watch()  # Live monitoring
```

### 5. Keep Old Workflow Available
Keep your old scripts! The `main()` function still works:
```r
# If you need the old workflow
source("main.R")
results <- main()
```

## Troubleshooting Migration

### "Target X not found"

Make sure you're in the project directory:
```r
setwd("/path/to/protein-neighbours")
tar_make()
```

### "Config file not found"

The pipeline looks for `config/config.yaml`. Make sure it exists:
```r
file.exists("config/config.yaml")  # Should be TRUE
```

### "Targets not updating after code changes"

Force invalidation:
```r
tar_invalidate(everything())
tar_make()
```

### "Out of disk space"

Clean old targets:
```r
tar_prune()  # Remove obsolete targets
```

Or clean everything:
```r
tar_destroy()  # Delete entire cache
```

## Getting Help

- Check [TARGETS_GUIDE.md](TARGETS_GUIDE.md) for detailed documentation
- See [TARGETS_QUICK_REFERENCE.md](TARGETS_QUICK_REFERENCE.md) for commands
- Run `?targets::tar_make` for help
- Visit https://books.ropensci.org/targets/ for the official manual

## Still Prefer the Old Way?

That's fine! The `main()` function still works:

```r
source("main.R")
results <- main()
```

You can migrate later when you're ready.

## Summary

| Task | Traditional | Targets |
|------|------------|---------|
| Run pipeline | `main()` | `tar_make()` |
| Load results | `results <- main()` | `tar_read(analysis_summary)` |
| Custom config | `main("config.yaml")` | Edit config, then `tar_make()` |
| Check progress | Wait and watch | `tar_progress()` |
| Visualize | N/A | `tar_visnetwork()` |
| Load intermediate | Can't | `tar_read(target_name)` |
| Re-run | Re-run everything | Only changed parts |
| Parallel | No | `tar_make_future()` |

## Next Steps

1. ✅ Read this guide
2. ✅ Try a simple `tar_make()`
3. ✅ Explore with `tar_visnetwork()`
4. ✅ Load some results with `tar_read()`
5. ✅ Read [TARGETS_GUIDE.md](TARGETS_GUIDE.md) for more

Happy analyzing! 🧬
