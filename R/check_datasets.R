# Ensure necessary packages are listed in DESCRIPTION Imports:
# yaml, readr, utils, stats, methods, tools (if file_path_sans_ext used elsewhere, maybe not here)

#' Check Data Structure, Prepare SUV Data, and Log Issues
#'
#' Loads dictionaries and configuration, checks dataset structures against them,
#' processes SUV variables based on the best dictionary match (renames or subsets
#' columns), and logs operations and warnings to a timestamped file. Warnings are
#' also shown in the console via a calling handler.
#'
#' @param loaded_data List. The output from `load_datasets()`. Expected
#'   to be a named list where each element is a dataset (itself a list
#'   containing data frames like 'data' and 'data_suv_bi').
#' @param files_path Character string. Path to the directory provided by the user
#'   containing supporting files: the SUV dictionary (`dict_suv.csv`) and the
#'   main configuration file (specified by `config_filename`).
#' @param config_filename Character string. Name of the YAML configuration
#'   file expected within `files_path`. Defaults to "config.yaml".
#' @param log_directory Character string or NULL. Path to the **directory** where
#'   timestamped log files should be written. If NULL (default), logging to
#'   file is disabled. The directory will be created if it doesn't exist.
#' @param log_file_prefix Character string. Prefix for the generated log filename.
#'   Defaults to "validation_log_".
#' @param dict_suv_target_col Character string or NULL. Target column name in
#'  `dict_suv.csv` containing the *final desired* variable names for `data_suv_bi`.
#'  Defaults to the first column name of the loaded `dict_suv`.
#' @param dict_suv_potential_cols Character vector or NULL. Names of columns in
#'   `dict_suv.csv` that contain *potential sets* of variable names to match
#'   against the actual `data_suv_bi`. If NULL (default), all columns *except*
#'   `dict_suv_target_col` are used.
#' @param verbose Logical. If TRUE, print some progress messages to console and
#'   enable detailed messages within the log file. Defaults to TRUE. Note: Most
#'   detailed messages go only to the log file when active due to sink settings.
#'
#' @return A list containing:
#'   \describe{
#'     \item{data}{The potentially modified `loaded_data` list. Specifically,
#'       `data_suv_bi` data frames may have columns renamed or subsetted.}
#'     \item{config}{The loaded configuration list.}
#'     \item{validation_report}{A list detailing structural issues found (e.g.,
#'       missing data frames, missing columns relative to best match, columns
#'       discarded/renamed).}
#'     \item{log_file}{The full path to the generated log file, or NULL if
#'       logging was disabled.}
#'   }
#'
#' @details
#' Performs checks, preparation, and logging:
#' 1.  **Logging Setup:** Creates a timestamped log file in `log_directory`.
#'     Uses `sink(..., split=FALSE)` for messages, directing them primarily
#'     to the log file. Warnings are captured via `withCallingHandlers`,
#'     logged to the file, and also allowed to display on the console.
#' 2.  **Config/Dict Loading:** Loads YAML config and `dict_suv.csv` from `files_path`.
#'     Validates essential config/dict structure.
#' 3.  **Dataset Iteration:** Processes each dataset found in `loaded_data`.
#' 4.  **Log Content:** Includes R version, package version, input dataset names,
#'     best SUV match column, missing columns, and discarded/renamed SUV columns.
#' 5.  **`data` Check:** Verifies columns in `data` against `required_data_columns`
#'     specified in the config file.
#' 6.  **`data_suv_bi` Processing:**
#'     *   Finds the column in `dict_suv` (among `dict_suv_potential_cols`)
#'       that has the most column names in common with the actual `data_suv_bi`.
#'     *   Checks for columns expected in that best-matching dictionary column
#'       but missing in the actual `data_suv_bi`.
#'     *   If best match *is* the target column (`dict_suv_target_col`): Subsets
#'       `data_suv_bi` to keep only columns defined in the target column (and
#'       present in the data), logging discarded columns.
#'     *   If best match *is not* the target column: Renames columns in
#'       `data_suv_bi` based on the mapping from the best-match column to the
#'       target column in `dict_suv.csv`, checking for potential name clashes.
#'
#' @export
#' @importFrom utils packageName packageVersion sessionInfo head str
#' @importFrom stats setNames
#' @importFrom yaml read_yaml
#' @importFrom readr read_csv locale
check_datasets <- function(loaded_data,
                                   files_path,
                                   config_filename = "config.yaml",
                                   log_directory = NULL,
                                   log_file_prefix = "validation_log_",
                                   dict_suv_target_col = NULL,
                                   dict_suv_potential_cols = NULL,
                                   verbose = TRUE) {
  
  # --- Logging Setup Variables ---
  log_con <- NULL
  actual_log_path <- NULL
  sink_active <- FALSE # Flag to track if sink was successfully activated
  
  # --- Define the Warning Handler ---
  warning_handler <- function(w) {
    if (sink_active && !is.null(log_con) && isOpen(log_con)) {
      call_info <- ""
      call_obj <- try(conditionCall(w), silent = TRUE)
      if (!inherits(call_obj, "try-error") && !is.null(call_obj)) {
        call_info <- paste(" in", deparse(call_obj, width.cutoff = 60L)[1])
      }
      log_msg <- sprintf("[%s] WARNING%s: %s\n",
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                         call_info,
                         conditionMessage(w))
      cat(log_msg, file = log_con, append = TRUE)
    }
    # Let R's default handler print it to the console.
  }
  
  # --- Setup Cleanup on Exit ---
  on.exit({
    if (sink_active) {
      sink_closed_ok <- FALSE
      try({
        if (sink.number(type = "message") > 0) sink(type = "message")
        sink_closed_ok <- TRUE
      }, silent = TRUE)
      if(!sink_closed_ok && interactive()) warning("on.exit: Attempt to close message sink failed.", call. = FALSE)
    }
    if (!is.null(log_con) && inherits(log_con, "connection") && isOpen(log_con)) {
      try({
        cat(sprintf("[%s] --- Log End ---\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_con, append = TRUE)
        close(log_con)
        if(verbose && interactive()) message("Log file closed: ", actual_log_path)
      }, silent = TRUE)
    }
  }, add = TRUE)
  
  # --- Determine Log File Path & Activate Logging ---
  # (Using the robust version from previous steps)
  if (!is.null(log_directory)) {
    if (!is.character(log_directory) || length(log_directory) != 1) {
      warning("`log_directory` must be a single character string. Disabling file logging.", call. = FALSE)
      log_directory <- NULL
    } else {
      path_exists <- file.exists(log_directory)
      if (path_exists) {
        is_dir <- suppressWarnings(file.info(log_directory)$isdir)
        if (!isTRUE(is_dir)) {
          warning("Provided log path '", log_directory, "' exists but is not a directory. Disabling file logging.", call. = FALSE)
          log_directory <- NULL
        }
      } else {
        if(verbose) message("Log directory not found, attempting to create: ", log_directory)
        tryCatch({
          dir.create(log_directory, recursive = TRUE, showWarnings = TRUE)
          if (!dir.exists(log_directory)) stop("Creation attempt failed.")
        }, error = function(e) {
          warning("Failed to create log directory '", log_directory, "'. Disabling file logging. Error: ", e$message, call. = FALSE)
          log_directory <- NULL
        })
      }
    }
    if (!is.null(log_directory)) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      log_filename <- paste0(log_file_prefix, timestamp, ".log")
      actual_log_path <- file.path(log_directory, log_filename)
      tryCatch({
        log_con <- file(actual_log_path, open = "at", encoding = "UTF-8")
        sink(log_con, type = "message", split = FALSE)
        sink_active <- TRUE
        message("--- Log Start ---")
        message("Logging detailed messages to: ", actual_log_path)
        message("(Warnings generated during processing will also appear in console)")
      }, error = function(e) {
        warning("Failed to open or sink to log file '", actual_log_path, "'. Logging to file disabled. Error: ", e$message, call. = FALSE)
        if(!is.null(log_con) && isOpen(log_con)) try(close(log_con), silent = TRUE)
        log_con <- NULL; actual_log_path <- NULL
        if(sink_active && sink.number("message") > 0) try(sink(type="message"), silent=TRUE)
        sink_active <- FALSE
      })
    }
  } # End logging setup
  
  # --- Log Session Info ---
  if (sink_active && verbose) {
    message("\n--- Session Info ---")
    message("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    si <- utils::sessionInfo()
    message("R Version: ", si$R.version$string)
    message("Platform: ", si$platform)
    pkg_name_env <- try(utils::packageName(), silent = TRUE)
    if (inherits(pkg_name_env, "try-error") || is.null(pkg_name_env)) {
      pkg_name <- "biodiscvr" # Hardcoded fallback - REPLACE if needed
      message("Package Version (", pkg_name, "): Could not determine dynamically.")
    } else {
      pkg_name <- pkg_name_env
      pkg_version_str <- try(as.character(utils::packageVersion(pkg_name)), silent = TRUE)
      if (inherits(pkg_version_str, "try-error")) message("Package Version (", pkg_name, "): Could not determine.")
      else message("Package Version (", pkg_name, "): ", pkg_version_str)
    }
    if(!is.null(loaded_data) && is.list(loaded_data)) message("Input Datasets: ", paste(names(loaded_data), collapse=", "))
    message("--------------------\n")
  }
  
  # --- Input Validation ---
  if (!is.list(loaded_data)) stop("Input 'loaded_data' must be a list.")
  if (!is.character(files_path) || length(files_path) != 1 || !dir.exists(files_path)) {
    stop("'files_path' must be a single string pointing to an existing directory.")
  }
  # Add more checks if needed...
  
  # --- Load Configuration (YAML) ---
  config_path <- file.path(files_path, config_filename)
  config_params <- list()
  if (!file.exists(config_path)) stop("Configuration file not found: ", config_path)
  tryCatch({
    if (sink_active && verbose) message("Loading configuration from: ", config_path)
    config_lines <- readLines(config_path, warn = FALSE)
    if (length(config_lines) > 0 && (nchar(config_lines[[length(config_lines)]]) > 0 && !endsWith(config_lines[[length(config_lines)]], "\n"))) {
      warning("Configuration file '", config_path, "' might be missing a final newline.", call.=FALSE)
    }
    config_params <- yaml::read_yaml(config_path)
    if (sink_active && verbose) message("  -> Loaded ", length(config_params), " config parameters.")
  }, error = function(e) {
    stop("Error parsing YAML configuration file '", config_path, "': ", e$message, call. = FALSE)
  })
  # --- Check for required_data_columns in config ---
  if (!"required_data_columns" %in% names(config_params) ||
      !is.vector(config_params$required_data_columns) ||
      length(config_params$required_data_columns) == 0) {
    stop("Config file ('", config_filename, "') must contain a non-empty list/vector named 'required_data_columns'.")
  }
  expected_data_vars <- unique(config_params$required_data_columns)
  
  # --- Load SUV Dictionary ---
  dict_suv_path <- file.path(files_path, "dict_suv.csv")
  dict_suv <- NULL
  if (!file.exists(dict_suv_path)) stop("SUV dictionary not found: ", dict_suv_path)
  tryCatch({
    if (sink_active && verbose) message("Loading SUV dictionary: ", dict_suv_path)
    dict_suv <- readr::read_csv(dict_suv_path, show_col_types = FALSE, progress = FALSE, locale = readr::locale(encoding = "UTF-8"))
    if (sink_active && verbose) message("  -> Success. Dictionary columns: ", paste(names(dict_suv), collapse=", "))
  }, error = function(e) {
    stop("Failed to load SUV dictionary '", dict_suv_path, "': ", e$message, call. = FALSE)
  })
  
  # --- Set Dictionary Column Defaults & Validate ---
  if (is.null(dict_suv_target_col)) {
    if (ncol(dict_suv) < 1) stop("SUV dictionary has no columns.")
    dict_suv_target_col <- names(dict_suv)[1]
    if (sink_active && verbose) message("Using first column ('", dict_suv_target_col, "') of SUV dictionary as target for renaming/subsetting.")
  } else if (!dict_suv_target_col %in% names(dict_suv)) {
    stop("Specified `dict_suv_target_col` ('", dict_suv_target_col, "') not found in SUV dictionary.")
  }
  
  if (is.null(dict_suv_potential_cols)) {
    if (ncol(dict_suv) < 2) {
      if(sink_active && verbose) message("SUV dictionary has only one column ('", dict_suv_target_col,"'). This will be used for matching.")
      dict_suv_potential_cols <- dict_suv_target_col
    } else {
      dict_suv_potential_cols <- setdiff(names(dict_suv), dict_suv_target_col)
      if (sink_active && verbose) message("Using columns as potential matches in SUV dictionary: ", paste(dict_suv_potential_cols, collapse = ", "))
    }
  } else {
    missing_potential <- setdiff(dict_suv_potential_cols, names(dict_suv))
    if (length(missing_potential) > 0) {
      stop("Specified `dict_suv_potential_cols` not found in SUV dictionary: ", paste(missing_potential, collapse=", "))
    }
    # Ensure target column is not accidentally in potential columns
    dict_suv_potential_cols <- setdiff(dict_suv_potential_cols, dict_suv_target_col)
    if (length(dict_suv_potential_cols) == 0) {
      stop("After validation, no potential match columns remain in `dict_suv_potential_cols`.")
    }
  }
  
  # --- Initialize Validation Report ---
  validation_report <- list()
  
  # --- Prepare Data Copy ---
  data_modified <- loaded_data # Work on a copy
  
  # --- Wrap the main processing loop with the warning handler ---
  if (sink_active && verbose) message("\n--- Starting Dataset Checks & Preparation (Messages below go to log file) ---")
  else if (verbose) cat("\n--- Starting Dataset Checks & Preparation ---\n")
  
  withCallingHandlers({
    for (dataset_name in names(data_modified)) {
      if (sink_active && verbose) message("\nProcessing dataset: ", dataset_name)
      else if(verbose) cat("\nProcessing dataset: ", dataset_name, "\n")
      
      current_dataset_list <- data_modified[[dataset_name]] # Get list for current dataset
      dataset_report <- list()
      
      # --- Check 'data' columns ---
      if ("data" %in% names(current_dataset_list) && !is.null(current_dataset_list$data)) {
        data_df <- current_dataset_list$data # Assign to temporary var
        actual_data_vars <- names(data_df)
        missing_data_vars <- setdiff(expected_data_vars, actual_data_vars)
        
        if (length(missing_data_vars) > 0) {
          msg <- sprintf("Dataset '%s' (data): Missing required columns defined in config ('%s'): %s",
                         dataset_name, config_filename, paste(missing_data_vars, collapse = ", "))
          warning(msg, call. = FALSE) # Handled by warning_handler
          dataset_report$missing_data_columns <- missing_data_vars
        } else {
          if (sink_active && verbose) message("  -> 'data': Columns match required list in config.")
        }
      } else {
        msg <- sprintf("Dataset '%s': 'data' data frame not found or is NULL.", dataset_name)
        warning(msg, call. = FALSE) # Handled by warning_handler
        dataset_report$missing_data_df <- TRUE
      }
      
      # --- Check & Prepare 'data_suv_bi' columns ---
      best_match_col_name <- NULL
      expected_suv_vars_from_match <- NULL # Store expected vars based on best match
      
      if ("data_suv_bi" %in% names(current_dataset_list) && !is.null(current_dataset_list$data_suv_bi)) {
        # Work on the actual data_suv_bi within the list we are modifying
        data_suv_bi_df <- data_modified[[dataset_name]]$data_suv_bi
        actual_suv_vars <- names(data_suv_bi_df)
        
        # --- Find best matching column in dict_suv ---
        best_match_score <- -1
        if (length(dict_suv_potential_cols) == 0){
          # This case might arise if dict only had one column initially
          if(ncol(dict_suv) == 1 && names(dict_suv)[1] == dict_suv_target_col) {
            if(sink_active && verbose) message("  -> 'data_suv_bi': Only target column available in dictionary. Using it for matching.")
            best_match_col_name <- dict_suv_target_col
            best_match_score <- length(intersect(actual_suv_vars, unique(dict_suv[[best_match_col_name]])))
          } else {
            # Should have been caught earlier, but safeguard
            warning("Dataset '", dataset_name, "': No potential columns available in SUV dictionary to perform best match.", call. = FALSE)
          }
        } else {
          # Loop through potential columns
          for (potential_col in dict_suv_potential_cols) {
            potential_vars <- unique(stats::na.omit(dict_suv[[potential_col]]))
            potential_vars <- potential_vars[nzchar(potential_vars)]
            if (length(potential_vars) == 0) next # Skip empty dict columns
            match_count <- length(intersect(actual_suv_vars, potential_vars))
            if (sink_active && verbose) message("    - Comparing vs dict column '", potential_col, "' (", length(potential_vars), " entries): ", match_count, " matches found.")
            if (match_count > best_match_score) {
              best_match_score <- match_count
              best_match_col_name <- potential_col
            }
          }
        } # End finding best match
        
        # --- Process based on best match ---
        if (is.null(best_match_col_name)) {
          msg <- sprintf("Dataset '%s' (data_suv_bi): Could not determine a best matching column in '%s'. Skipping checks and processing.",
                         dataset_name, basename(dict_suv_path))
          warning(msg, call. = FALSE) # Handled by warning_handler
          dataset_report$suv_no_match <- TRUE
        } else {
          if (sink_active && verbose) message("  -> 'data_suv_bi': Best match identified as dictionary column '", best_match_col_name, "' (", best_match_score, " matches).")
          dataset_report$suv_best_match_column <- best_match_col_name
          
          # --- Check missing based on best match ---
          expected_suv_vars_from_match <- unique(stats::na.omit(dict_suv[[best_match_col_name]]))
          expected_suv_vars_from_match <- expected_suv_vars_from_match[nzchar(expected_suv_vars_from_match)]
          missing_suv_vars <- setdiff(expected_suv_vars_from_match, actual_suv_vars)
          
          if (length(missing_suv_vars) > 0) {
            msg <- sprintf("Dataset '%s' (data_suv_bi): Missing expected columns (based on best match '%s'): %s",
                           dataset_name, best_match_col_name, paste(missing_suv_vars, collapse = ", "))
            warning(msg, call. = FALSE) # Handled by warning_handler
            dataset_report$missing_suv_columns <- missing_suv_vars
          } else {
            if (sink_active && verbose) message("  -> 'data_suv_bi': All columns from best match ('", best_match_col_name, "') are present.")
          }
          
          # --- Apply Renaming OR Subsetting ---
          if (best_match_col_name == dict_suv_target_col) {
            # --- Subsetting Logic ---
            if (sink_active && verbose) message("  -> 'data_suv_bi': Best match is the target column ('", dict_suv_target_col, "'). Subsetting columns.")
            # expected_suv_vars_from_match already holds the target vars here
            vars_to_keep <- intersect(expected_suv_vars_from_match, actual_suv_vars)
            vars_to_discard <- setdiff(actual_suv_vars, expected_suv_vars_from_match)
            
            if (length(vars_to_keep) == 0) {
              msg <- sprintf("Dataset '%s' (data_suv_bi): Subsetting based on target column '%s' resulted in NO columns to keep.", dataset_name, dict_suv_target_col)
              warning(msg, call. = FALSE) # Handled by warning_handler
              dataset_report$suv_subset_empty <- TRUE
              # Keep empty dataframe structure
              data_modified[[dataset_name]]$data_suv_bi <- data_suv_bi_df[, vars_to_keep, drop = FALSE]
            } else {
              if (length(vars_to_discard) > 0) {
                if (sink_active && verbose) message("  -> 'data_suv_bi': Discarding columns not in target list ('", dict_suv_target_col, "'): ", paste(vars_to_discard, collapse = ", "))
                dataset_report$discarded_suv_columns <- vars_to_discard
              } else {
                if (sink_active && verbose) message("  -> 'data_suv_bi': All existing columns are defined in the target list. No columns discarded.")
              }
              # *** Perform the subsetting - modify the df IN THE LIST ***
              data_modified[[dataset_name]]$data_suv_bi <- data_suv_bi_df[, vars_to_keep, drop = FALSE]
            }
            
          } else {
            # --- Renaming Logic ---
            if (sink_active && verbose) message("  -> 'data_suv_bi': Best match ('", best_match_col_name, "') differs from target ('", dict_suv_target_col,"'). Applying renaming.")
            mapping_df <- unique(dict_suv[, c(best_match_col_name, dict_suv_target_col)])
            # Filter map based on non-empty/NA source keys
            mapping_df <- mapping_df[!is.na(mapping_df[[best_match_col_name]]) & nzchar(mapping_df[[best_match_col_name]]), ]
            if (anyDuplicated(mapping_df[[best_match_col_name]])) {
              warning("Dataset '", dataset_name, "': Duplicates found in source column '", best_match_col_name, "' of SUV dictionary map. Using first occurrence.", call.=FALSE) # Handled
              mapping_df <- mapping_df[!duplicated(mapping_df[[best_match_col_name]]), ]
            }
            # Ensure target names in map are also valid
            mapping_df <- mapping_df[!is.na(mapping_df[[dict_suv_target_col]]) & nzchar(mapping_df[[dict_suv_target_col]]), ]
            
            renaming_map <- stats::setNames(mapping_df[[dict_suv_target_col]], mapping_df[[best_match_col_name]])
            current_names <- names(data_suv_bi_df) # Get current names of df being modified
            names_to_rename_from <- intersect(current_names, names(renaming_map))
            
            if (length(names_to_rename_from) > 0) {
              new_names_for_these <- renaming_map[names_to_rename_from]
              final_names <- current_names
              match_indices <- match(names_to_rename_from, final_names)
              potential_new_set <- final_names; potential_new_set[match_indices] <- new_names_for_these
              
              if(anyDuplicated(potential_new_set)) {
                warning("Dataset '", dataset_name, "': Renaming SUV columns would create duplicate names. Skipping renaming.", call. = FALSE) # Handled
                dataset_report$suv_rename_skipped_duplicates <- TRUE
              } else {
                final_names[match_indices] <- new_names_for_these
                # *** Apply renaming - modify the df IN THE LIST ***
                names(data_modified[[dataset_name]]$data_suv_bi) <- final_names
                rename_log_pairs <- stats::setNames(new_names_for_these, names_to_rename_from)
                if (sink_active && verbose) {
                  rename_log_text <- paste(paste0("'", names(rename_log_pairs), "' -> '", rename_log_pairs, "'"), collapse = "; ")
                  message("  -> 'data_suv_bi': Renamed columns: ", rename_log_text)
                }
                dataset_report$renamed_suv_columns <- rename_log_pairs
              }
            } else {
              if (sink_active && verbose) message("  -> 'data_suv_bi': No columns needed renaming based on best match->target mapping.")
            }
          } # End Renaming/Subsetting decision
        } # End if best_match_col_name found
      } else {
        msg <- sprintf("Dataset '%s': 'data_suv_bi' data frame not found or is NULL. Skipping SUV checks and processing.", dataset_name)
        warning(msg, call. = FALSE) # Handled by warning_handler
        dataset_report$missing_suv_df <- TRUE
      } # End check for data_suv_bi existence
      
      # --- Store report ---
      if (length(dataset_report) > 0) {
        validation_report[[dataset_name]] <- dataset_report
      }
      # NOTE: data_modified[[dataset_name]] was updated directly in renaming/subsetting
      
    } # End loop through datasets
  }, warning = warning_handler) # End of withCallingHandlers Block
  
  # --- Final message ---
  if (sink_active && verbose) message("\n--- Finished Dataset Checks & Preparation ---")
  else if(verbose) cat("\n--- Finished Dataset Checks & Preparation ---\n")
  
  # --- Return Results ---
  return(list(
    data = data_modified, # Return the potentially modified data list
    config = config_params,
    validation_report = validation_report,
    log_file = actual_log_path
  ))
}