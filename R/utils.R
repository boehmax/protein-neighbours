#' Load configuration from YAML file
#'
#' This function loads configuration parameters from a YAML file.
#' If a configuration file is not provided, it uses the default configuration.
#'
#' @param config_file Path to the configuration YAML file
#' @param override_params Named list of parameters that override the config file values
#' @return A list containing the configuration parameters
#' @importFrom yaml read_yaml
#' @importFrom utils str
#' @export
load_config <- function(config_file = "config/config.yaml", override_params = NULL) {
  # Check if yaml package is installed
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is needed for configuration. Please install it.")
  }
  
  # Check if config file exists, otherwise use default
  if (!file.exists(config_file)) {
    warning(paste("Configuration file", config_file, "not found. Using default configuration."))
    config_file <- system.file("config", "default_config.yaml", package = "proteinNeighbours")
    
    # If default config is not found, throw an error
    if (!file.exists(config_file)) {
      stop("Default configuration file not found. Please reinstall the package or provide a valid configuration file.")
    }
  }
  
  # Load the configuration from YAML file
  config <- yaml::read_yaml(config_file)
  
  # Override configuration parameters if provided
  if (!is.null(override_params)) {
    config <- override_config(config, override_params)
  }
  
  # Set date if not specified
  if (is.null(config$analysis$date)) {
    config$analysis$date <- format(Sys.Date(), "%Y-%m-%d")
  }
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(config$paths$output_dir, config$analysis$date)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the configuration for reproducibility
  save_config(config, output_dir)
  
  return(config)
}

#' Override configuration parameters
#'
#' This function overrides the configuration with user-provided parameters.
#'
#' @param config The configuration list
#' @param override_params Named list of parameters that override the config
#' @return The modified configuration list
#' @keywords internal
override_config <- function(config, override_params) {
  # Loop through override parameters and update config
  for (param_name in names(override_params)) {
    # Handle nested parameters with dot notation (e.g., "analysis.basepairs")
    if (grepl("\\.", param_name)) {
      parts <- strsplit(param_name, "\\.")[[1]]
      
      # Navigate to the correct level in the config
      current_level <- config
      for (i in 1:(length(parts) - 1)) {
        if (is.null(current_level[[parts[i]]])) {
          current_level[[parts[i]]] <- list()
        }
        current_level <- current_level[[parts[i]]]
      }
      
      # Set the parameter value
      current_level[[parts[length(parts)]]] <- override_params[[param_name]]
    } else {
      # For top-level parameters
      config[[param_name]] <- override_params[[param_name]]
    }
  }
  
  return(config)
}

#' Save configuration for reproducibility
#'
#' This function saves the current configuration to a YAML file in the output directory.
#'
#' @param config The configuration list
#' @param output_dir The output directory
#' @importFrom yaml write_yaml
#' @keywords internal
save_config <- function(config, output_dir) {
  # Save the configuration to the output directory
  yaml::write_yaml(config, file.path(output_dir, "analysis_config.yaml"))
}

#' Log input file information for reproducibility
#'
#' This function logs information about input files for reproducibility.
#'
#' @param config The configuration list
#' @return A data frame with input file information
#' @export
pn_input_files <- function(config) {
  # Create a summary data frame for input files
  input_files <- c(
    file.path(config$paths$base_dir, config$files$proteins),
    file.path(config$paths$base_dir, config$files$assemblies),
    file.path(config$paths$base_dir, config$files$protein_assembly)
  )
  
  # Add representative files if they exist
  if (!is.null(config$files$representative_files)) {
    rep_files <- c(
      file.path(config$paths$base_dir, config$files$representative_files$ipg),
      file.path(config$paths$base_dir, config$files$representative_files$pdb),
      file.path(config$paths$base_dir, config$files$representative_files$cluster)
    )
    input_files <- c(input_files, rep_files)
  }
  
  # Get file information
  file_info <- data.frame(
    File = input_files,
    Exists = file.exists(input_files),
    Size_bytes = sapply(input_files, function(f) ifelse(file.exists(f), file.size(f), NA)),
    Modified = sapply(input_files, function(f) ifelse(file.exists(f), 
                                                      format(file.mtime(f), "%Y-%m-%d %H:%M:%S"), NA)),
    stringsAsFactors = FALSE
  )
  
  # Add number of lines for text files
  file_info$Lines <- sapply(input_files, function(f) {
    if (file.exists(f) && file.size(f) < 10e6) {  # Only read if file is smaller than 10MB
      tryCatch(length(readLines(f)), error = function(e) NA)
    } else {
      NA
    }
  })
  
  # Save file info to the output directory
  output_dir <- file.path(config$paths$output_dir, config$analysis$date)
  write.csv(file_info, file.path(output_dir, "input_file_info.csv"), row.names = FALSE)
  
  # Log file information
  pn_info("Input file information saved to input_file_info.csv")
  for (i in 1:nrow(file_info)) {
    if (file_info$Exists[i]) {
      pn_info(sprintf("File: %s, Size: %d bytes, Modified: %s", 
                      file_info$File[i], file_info$Size_bytes[i], file_info$Modified[i]))
    } else {
      pn_warn(sprintf("File not found: %s", file_info$File[i]))
    }
  }
  
  return(file_info)
}

#' Set up standalone logging
#'
#' This function sets up a standalone logging system with no external dependencies
#'
#' @param config The configuration list
#' @return The log file path
#' @export
pn_setup_logging <- function(config) {
  # Set output directory
  output_dir <- file.path(config$paths$output_dir, config$analysis$date)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set log file
  log_file <- file.path(output_dir, config$logging$file)
  
  # Store log level and file as options (not in global environment)
  options(pn_log_level = toupper(config$logging$level))
  options(pn_log_file = log_file)
  
  # Initialize log file
  write(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "[INFO]", 
              "Starting protein neighborhood analysis"), log_file)
  write(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "[INFO]", 
              "Configuration loaded"), log_file, append = TRUE)
  
  # Return log file path
  message("Logging to: ", log_file)
  return(log_file)
}

#' Get the numeric priority of a log level
#'
#' @param level The log level string
#' @return Numeric priority (higher = more severe)
#' @keywords internal
pn_log_level_priority <- function(level) {
  level_priorities <- c("DEBUG" = 1, "INFO" = 2, "WARN" = 3, "WARNING" = 3, 
                        "ERROR" = 4, "FATAL" = 5)
  return(level_priorities[toupper(level)])
}

#' Log a message with a specific level
#'
#' @param level The log level (DEBUG, INFO, WARN, ERROR, FATAL)
#' @param ... The message components to log
#' @keywords internal
pn_log_message <- function(level, ...) {
  # Get current log level from options
  current_level <- getOption("pn_log_level", "INFO")
  log_file <- getOption("pn_log_file", NULL)
  
  # Check if message should be logged based on level
  if (pn_log_level_priority(level) >= pn_log_level_priority(current_level)) {
    # Format the message
    msg <- paste(..., collapse = " ")
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full_msg <- paste(timestamp, "[", level, "]", msg)
    
    # Print to console
    cat(full_msg, "\n")
    
    # Write to log file if available
    if (!is.null(log_file)) {
      write(full_msg, log_file, append = TRUE)
    }
  }
}

#' Log an info message
#' 
#' @param ... The message components to log
#' @export
pn_info <- function(...) {
  pn_log_message("INFO", ...)
}

#' Log a warning message
#' 
#' @param ... The message components to log
#' @export
pn_warn <- function(...) {
  pn_log_message("WARN", ...)
}

#' Log an error message
#' 
#' @param ... The message components to log
#' @export
pn_error <- function(...) {
  pn_log_message("ERROR", ...)
}

#' Log a debug message
#' 
#' @param ... The message components to log
#' @export
pn_debug <- function(...) {
  pn_log_message("DEBUG", ...)
}

# Find the script path properly and make it executable

#' Get the path to a script file
#'
#' This function finds the path to a script file, trying multiple locations.
#'
#' @param script_name The name of the script file
#' @param default_locations Additional locations to check
#' @return The full path to the script file, or NULL if not found
#' @keywords internal
find_script <- function(script_name, default_locations = c()) {
  # First, check if it's in the package's scripts directory
  script_path <- NULL
  
  # Try system.file first (works when installed as a package)
  potential_path <- system.file("scripts", script_name, package = "proteinNeighbours")
  if (file.exists(potential_path)) {
    script_path <- potential_path
  }
  
  # Try relative to the current directory (development mode)
  if (is.null(script_path)) {
    potential_paths <- c(
      file.path("inst", "scripts", script_name),
      file.path("scripts", script_name),
      script_name
    )
    
    # Add any additional locations
    potential_paths <- c(potential_paths, file.path(default_locations, script_name))
    
    # Check each potential path
    for (path in potential_paths) {
      if (file.exists(path)) {
        script_path <- normalizePath(path)
        break
      }
    }
  }
  
  if (is.null(script_path)) {
    warning("Could not find script: ", script_name)
    return(NULL)
  }
  
  # Make the script executable if it's not already
  if (.Platform$OS.type != "windows") {
    system(paste("chmod +x", shQuote(script_path)))
  }
  
  return(script_path)
}

#' Rounding function
#'
#' @param x The number to be rounded.
#' @param accuracy The accuracy to which the number should be rounded.
#' @param f The rounding function to use (default is round).
#' @return The rounded number.
#' @export
round_any <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}