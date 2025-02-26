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

#' Setup logging for the analysis
#'
#' This function sets up logging for the analysis.
#'
#' @param config The configuration list
#' @return The configured logger
#' @importFrom logger log_levels log_threshold log_formatter log_appender
#' @keywords internal
setup_logging <- function(config) {
  # Check if logger package is installed
  if (!requireNamespace("logger", quietly = TRUE)) {
    warning("Package 'logger' is needed for logging. Basic logging will be used instead.")
    return(NULL)
  }
  
  # Set log level
  level <- config$logging$level
  logger::log_threshold(logger::log_levels[[level]])
  
  # Set log format
  logger::log_formatter(formatter_with_timestamp)
  
  # Set log file
  log_file <- file.path(config$paths$output_dir, config$analysis$date, config$logging$file)
  logger::log_appender(logger::appender_file(log_file))
  
  # Log the start of the analysis
  logger::log_info("Starting protein neighborhood analysis")
  logger::log_info(paste("Configuration loaded from:", config))
  
  return(logger)
}

#' Custom log formatter with timestamp
#'
#' @param level Log level
#' @param msg Log message
#' @return Formatted log message
#' @keywords internal
formatter_with_timestamp <- function(level, msg) {
  paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " [", level, "] ", msg
  )
}

#' Log input file information for reproducibility
#'
#' This function logs information about input files for reproducibility.
#'
#' @param config The configuration list
#' @importFrom logger log_info
#' @export
log_input_files <- function(config) {
  # Create a summary data frame for input files
  input_files <- c(
    file.path(config$paths$base_dir, config$files$proteins),
    file.path(config$paths$base_dir, config$files$assemblies),
    file.path(config$paths$base_dir, config$files$protein_assembly)
  )
  
  # Add representative files if they exist
  rep_files <- c(
    file.path(config$paths$base_dir, config$files$representative_files$ipg),
    file.path(config$paths$base_dir, config$files$representative_files$pdb),
    file.path(config$paths$base_dir, config$files$representative_files$cluster)
  )
  input_files <- c(input_files, rep_files)
  
  # Get file information
  file_info <- data.frame(
    File = input_files,
    Exists = file.exists(input_files),
    Size_bytes = sapply(input_files, function(f) ifelse(file.exists(f), file.size(f), NA)),
    Modified = sapply(input_files, function(f) ifelse(file.exists(f), format(file.mtime(f), "%Y-%m-%d %H:%M:%S"), NA)),
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
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info("Input file information saved to input_file_info.csv")
    for (i in 1:nrow(file_info)) {
      if (file_info$Exists[i]) {
        logger::log_info(sprintf("File: %s, Size: %d bytes, Modified: %s", 
                                  file_info$File[i], file_info$Size_bytes[i], file_info$Modified[i]))
      } else {
        logger::log_warn(sprintf("File not found: %s", file_info$File[i]))
      }
    }
  }
  
  return(file_info)
}
