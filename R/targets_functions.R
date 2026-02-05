#' Targets-Specific Wrapper Functions
#'
#' This file contains wrapper functions and utilities specifically
#' designed for use with the targets pipeline.
#'
#' @author Maximilian Böhm
#' @version 0.2.0

#' Run the complete targets pipeline
#'
#' This is a convenience function that runs the entire targets pipeline.
#' It's equivalent to calling targets::tar_make().
#'
#' @param callr_function The callr function to use for running targets.
#'   Default is NULL which uses targets default (callr::r).
#'   Set to NULL or callr::r_bg for background execution.
#' @param reporter Reporter to use for progress updates.
#'   Options: "verbose", "timestamp", "summary", "forecast"
#' @return Invisibly returns the result of tar_make()
#' @export
#' @examples
#' \dontrun{
#' # Run the pipeline
#' run_targets_pipeline()
#'
#' # Run with verbose output
#' run_targets_pipeline(reporter = "verbose")
#' }
run_targets_pipeline <- function(callr_function = NULL, reporter = "timestamp") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  targets::tar_make(
    callr_function = callr_function,
    reporter = reporter
  )
}

#' Visualize the targets pipeline
#'
#' Creates an interactive network visualization of the pipeline dependency graph.
#'
#' @param targets_only Logical. If TRUE, only show targets (exclude functions).
#' @param outdated Logical. If TRUE, highlight outdated targets.
#' @return A visNetwork htmlwidget object.
#' @export
#' @examples
#' \dontrun{
#' # Visualize the pipeline
#' visualize_pipeline()
#'
#' # Show only outdated targets
#' visualize_pipeline(outdated = TRUE)
#' }
visualize_pipeline <- function(targets_only = TRUE, outdated = FALSE) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  targets::tar_visnetwork(
    targets_only = targets_only,
    label = c("time", "size", "branches")
  )
}

#' Get pipeline progress summary
#'
#' Returns a summary of pipeline execution status.
#'
#' @return A data frame with target names, statuses, and other metadata.
#' @export
#' @examples
#' \dontrun{
#' # Check pipeline progress
#' get_pipeline_progress()
#' }
get_pipeline_progress <- function() {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  targets::tar_progress()
}

#' Load a specific target's result
#'
#' Convenience wrapper around tar_read() with better error messages.
#'
#' @param target_name Name of the target to load (as a string or symbol).
#' @return The value stored in the target.
#' @export
#' @examples
#' \dontrun{
#' # Load the combined data frame
#' df <- load_target("combined_df")
#'
#' # Load the configuration
#' config <- load_target("config")
#' }
load_target <- function(target_name) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  tryCatch({
    targets::tar_read(!!rlang::sym(target_name))
  }, error = function(e) {
    stop(paste0(
      "Failed to load target '", target_name, "'.\n",
      "Make sure the pipeline has been run with run_targets_pipeline() first.\n",
      "Original error: ", e$message
    ))
  })
}

#' Clean the targets cache
#'
#' Removes the _targets/ directory to start fresh.
#'
#' @param ask Logical. If TRUE, asks for confirmation before cleaning.
#' @return Invisibly returns TRUE if successful.
#' @export
#' @examples
#' \dontrun{
#' # Clean the cache
#' clean_targets_cache()
#'
#' # Clean without asking
#' clean_targets_cache(ask = FALSE)
#' }
clean_targets_cache <- function(ask = TRUE) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  if (ask) {
    response <- readline(prompt = "Are you sure you want to delete the targets cache? (yes/no): ")
    if (tolower(response) != "yes") {
      message("Cache cleaning cancelled.")
      return(invisible(FALSE))
    }
  }
  
  targets::tar_destroy(ask = FALSE)
  message("Targets cache cleaned successfully.")
  invisible(TRUE)
}

#' Invalidate specific targets
#'
#' Marks specific targets as outdated so they will be re-run.
#'
#' @param target_names Character vector of target names to invalidate.
#' @return Invisibly returns TRUE if successful.
#' @export
#' @examples
#' \dontrun{
#' # Invalidate the annotation step
#' invalidate_targets("annotation_results")
#'
#' # Invalidate multiple targets
#' invalidate_targets(c("annotation_results", "combined_df"))
#' }
invalidate_targets <- function(target_names) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  targets::tar_invalidate(matches(target_names))
  message(paste("Invalidated targets:", paste(target_names, collapse = ", ")))
  invisible(TRUE)
}

#' Check which targets are outdated
#'
#' Returns a character vector of target names that need to be updated.
#'
#' @return Character vector of outdated target names.
#' @export
#' @examples
#' \dontrun{
#' # Check outdated targets
#' outdated <- get_outdated_targets()
#' print(outdated)
#' }
get_outdated_targets <- function() {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  targets::tar_outdated()
}

#' Run targets pipeline with parallel execution
#'
#' Runs the pipeline using parallel workers for independent targets.
#'
#' @param workers Number of parallel workers to use.
#' @param reporter Reporter to use for progress updates.
#' @return Invisibly returns the result of tar_make_future().
#' @export
#' @examples
#' \dontrun{
#' # Run with 4 parallel workers
#' run_targets_parallel(workers = 4)
#' }
run_targets_parallel <- function(workers = 2, reporter = "timestamp") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for parallel execution. Install with: install.packages('future')")
  }
  
  # Set up parallel backend
  future::plan(future::multisession, workers = workers)
  
  # Run pipeline
  result <- targets::tar_make_future(
    reporter = reporter
  )
  
  # Reset to sequential
  future::plan(future::sequential)
  
  invisible(result)
}

#' Generate targets pipeline manifest
#'
#' Creates a data frame describing all targets in the pipeline.
#'
#' @return A data frame with target information.
#' @export
#' @examples
#' \dontrun{
#' # Get pipeline manifest
#' manifest <- get_pipeline_manifest()
#' print(manifest)
#' }
get_pipeline_manifest <- function() {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  targets::tar_manifest()
}

#' Create a targets pipeline report
#'
#' Generates a detailed HTML report of the pipeline execution.
#'
#' @param output_file Path to save the HTML report.
#' @return Path to the generated report file.
#' @export
#' @examples
#' \dontrun{
#' # Generate pipeline report
#' create_pipeline_report("pipeline_report.html")
#' }
create_pipeline_report <- function(output_file = "targets_report.html") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required. Please install it with: install.packages('targets')")
  }
  
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for report generation.")
  }
  
  # Create a temporary Rmd file for the report
  temp_rmd <- tempfile(fileext = ".Rmd")
  
  rmd_content <- '---
title: "Targets Pipeline Report"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(targets)
```

## Pipeline Overview

```{r}
tar_manifest()
```

## Pipeline Progress

```{r}
tar_progress()
```

## Pipeline Network

```{r}
tar_visnetwork()
```

## Outdated Targets

```{r}
tar_outdated()
```
'
  
  writeLines(rmd_content, temp_rmd)
  
  rmarkdown::render(
    input = temp_rmd,
    output_file = output_file,
    quiet = TRUE
  )
  
  message(paste("Pipeline report generated:", output_file))
  invisible(output_file)
}
