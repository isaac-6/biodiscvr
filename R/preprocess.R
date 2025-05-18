#' Preprocess Loaded Datasets
#'
#' Loads configuration and dictionaries, checks data structures, prepares SUV data
#' (renaming or subsetting columns based on dictionary match), filters datasets
#' based on minimum entries per ID and group-specific inclusion criteria, and
#' logs operations, warnings, and filtering counts.
#'
#' @param loaded_data List. The **raw** output from `load_datasets()`. Expected
#'   to be a named list where each element is a dataset (itself a list
#'   containing data frames like 'data' and 'data_suv_bi').
#' @param files_path Character string. Path to the directory provided by the user
#'   containing supporting files: the SUV dictionary (`dict_suv.csv`) and the
#'   main configuration file (specified by `config_filename`).
#' @param config_filename Character string. Name of the YAML configuration
#'   file expected within `files_path`. Defaults to "config.yaml".
#' @param dict_suv_filename Dictionary file to manage different feature naming conventions
#' @param log_directory Character string or NULL. Path to the **directory** where
#'   timestamped log files should be written. If NULL (default), logging to
#'   file is disabled. The directory will be created if it doesn't exist.
#' @param scandate_column Name of the column in data.csv containing the scan date in Y%m%d format
#' @param log_file_prefix Character string. Prefix for the generated log filename.
#'   Defaults to "preprocess_log_".
#' @param dict_suv_target_col Character string or NULL. Target column name in
#'  `dict_suv.csv` containing the *final desired* variable names for `data_suv_bi`.
#'  Defaults to the first column name of the loaded `dict_suv`.
#' @param dict_suv_potential_cols Character vector or NULL. Names of columns in
#'   `dict_suv.csv` that contain *potential sets* of variable names to match
#'   against the actual `data_suv_bi`. If NULL (default), **all** columns in
#'   `dict_suv` are used for matching.
#' @param verbose Logical. If TRUE, print some progress messages to console and
#'   enable detailed messages within the log file. Defaults to TRUE. Note: Most
#'   detailed messages go only to the log file when active due to sink settings.
#'
#' @return A list containing:
#'   \describe{
#'     \item{data}{The **processed** list of datasets. Data frames within each
#'       dataset will have `data_suv_bi` columns potentially renamed/subsetted
#'       and *all* data frames filtered row-wise based on preprocessing criteria.}
#'     \item{config}{The loaded configuration list.}
#'     \item{validation_report}{A list detailing structural issues found *before*
#'       row filtering (e.g., missing required columns, issues during SUV
#'       renaming/subsetting, discarded columns).}
#'     \item{preprocessing_log}{A list detailing the number of unique IDs and rows
#'       removed during the preprocessing filter steps (min entries, inclusion
#'       criteria) for each dataset.}
#'     \item{log_file}{The full path to the generated detailed text log file,
#'       or NULL if logging was disabled.}
#'   }
#'
#' @details
#' This function performs the entire check, preparation, and preprocessing workflow:
#' 1.  **Logging Setup:** Creates a timestamped log file in `log_directory` if specified.
#'     Uses `sink(..., split=FALSE)` for messages. Warnings are captured via
#'     `withCallingHandlers`, logged, and also displayed on the console.
#' 2.  **Config/Dict Loading:** Loads YAML config and `dict_suv.csv` from `files_path`.
#'     Validates essential configuration.
#' 3.  **Dataset Iteration:** Processes each dataset found in `loaded_data`.
#' 4.  **Initial Checks & SUV Prep (Per Dataset):**
#'     *   Verifies `data` columns against `required_data_columns` in config.
#'     *   Finds the column in `dict_suv` with the most matching columns vs `data_suv_bi`.
#'     *   Checks for missing columns relative to the best match.
#'     *   **Renames OR Subsets** columns in `data_suv_bi` based on best match vs target.
#'     *   Logs findings and actions in `validation_report` and the text log.
#' 5.  **Preprocessing Filters (Per Dataset):**
#'     *   Reads preprocessing settings (`min_entries_per_id`, etc.) and
#'       `inclusion_criteria` from config.
#'     *   Filters individuals based on minimum entries.
#'     *   Filters individuals based on group-specific inclusion criteria (AGE, MMSE, etc.).
#'       Handles missing criteria columns with warnings.
#'     *   Determines the final set of rows to keep based on individuals passing filters.
#'     *   **Applies row filtering** to *all* data frames within the dataset based
#'       on the indices determined from the criteria source file (e.g., 'data'),
#'       ensuring row correspondence. Issues warnings if row counts mismatch.
#'     *   Logs counts of removed IDs/rows at each step in `preprocessing_log` and the text log.
#'
#' @export
#' @importFrom utils packageName packageVersion sessionInfo head str
#' @importFrom stats setNames median quantile IQR sd na.omit aggregate ave
#' @importFrom methods is
#' @importFrom yaml read_yaml
#' @importFrom readr read_csv locale write_csv
#' @importFrom dplyr group_by summarise filter n select all_of pull between distinct %>% across sym
#' @importFrom rlang .data `%||%`
preprocess_data <- function(loaded_data,
                            files_path,
                            config_filename = "config.yaml",
                            dict_suv_filename = "dict_suv.csv",
                            log_directory = NULL,
                            scandate_column = "ScanDate",
                            log_file_prefix = "preprocess_log_",
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
      try(cat(log_msg, file = log_con, append = TRUE), silent=TRUE) # Safely try writing
    }
    # Let R's default handler print it to the console.
  }
  
  # --- Setup Cleanup on Exit ---
  on.exit({
    if (sink_active) {
      sink_closed_ok <- FALSE
      try({ if (sink.number(type = "message") > 0) sink(type = "message") ; sink_closed_ok <- TRUE }, silent = TRUE)
      # Optional: Warn if closing failed (might indicate external interference)
      # if(!sink_closed_ok && interactive()) warning("on.exit: Attempt to close message sink failed.", call. = FALSE)
    }
    if (!is.null(log_con) && inherits(log_con, "connection") && isOpen(log_con)) {
      try({
        cat(sprintf("[%s] --- Log End ---\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_con, append = TRUE)
        close(log_con)
        if(verbose && interactive()) message("Log file closed: ", actual_log_path)
      }, silent = TRUE)
    }
  }, add = TRUE) # Ensures it runs alongside other on.exit calls if any
  
  # --- Determine Log File Path & Activate Logging ---
  # (Using the robust version from previous answers)
  if (!is.null(log_directory)) {
    if (!is.character(log_directory) || length(log_directory) != 1) {
      warning("`log_directory` invalid. Disabling file logging.", call. = FALSE); log_directory <- NULL
    } else {
      path_exists <- base::file.exists(log_directory) # Use base explicitly if needed
      if (path_exists && !isTRUE(suppressWarnings(base::file.info(log_directory)$isdir))) {
        warning("Log path '", log_directory, "' exists but is not a directory. Disabling file logging.", call. = FALSE); log_directory <- NULL
      } else if (!path_exists) {
        if(verbose) message("Creating log directory: ", log_directory)
        tryCatch({ base::dir.create(log_directory, recursive = TRUE, showWarnings = TRUE); if (!base::dir.exists(log_directory)) stop("Creation failed.") },
                 error = function(e) { warning("Failed to create log directory '", log_directory, "'. Error: ", e$message, call. = FALSE); log_directory <- NULL })
      }
    }
    if (!is.null(log_directory)) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S"); log_filename <- paste0(log_file_prefix, timestamp, ".log"); actual_log_path <- base::file.path(log_directory, log_filename)
      tryCatch({
        log_con <- base::file(actual_log_path, open = "at", encoding = "UTF-8"); base::sink(log_con, type = "message", split = FALSE); sink_active <- TRUE
        message("--- Log Start ---"); message("Logging detailed messages to: ", actual_log_path); message("(Warnings generated during processing will also appear in console)")
      }, error = function(e) {
        warning("Failed to open/sink log file '", actual_log_path, "'. Error: ", e$message, call. = FALSE)
        if(!is.null(log_con) && isOpen(log_con)) try(close(log_con), silent = TRUE); log_con <- NULL; actual_log_path <- NULL
        if(sink_active && sink.number("message") > 0) try(sink(type="message"), silent=TRUE); sink_active <- FALSE
      })
    }
  } # End logging setup
  
  # --- Log Session Info ---
  if (sink_active && verbose) {
    try({ # Wrap logging just in case session info functions fail
      message("\n--- Session Info ---")
      message("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
      si <- utils::sessionInfo()
      message("R Version: ", si$R.version$version.string)
      message("Platform: ", si$platform)
      pkg_name_env <- try(utils::packageName(), silent = TRUE)
      if (inherits(pkg_name_env, "try-error") || is.null(pkg_name_env)) {
        pkg_name <- "biodiscvr" # Hardcoded fallback - !!! REPLACE if your package name is different !!!
        message("Package Version (", pkg_name, "): Could not determine dynamically.")
      } else {
        pkg_name <- pkg_name_env
        pkg_version_str <- try(as.character(utils::packageVersion(pkg_name)), silent = TRUE)
        if (inherits(pkg_version_str, "try-error")) message("Package Version (", pkg_name, "): Could not determine.")
        else message("Package Version (", pkg_name, "): ", pkg_version_str)
      }
      if(!is.null(loaded_data) && is.list(loaded_data)) message("Input Datasets: ", paste(names(loaded_data), collapse=", "))
      message("--------------------\n")
    }, silent=TRUE)
  }
  
  # --- Initial Input Validation ---
  if (!is.list(loaded_data)) stop("Input 'loaded_data' must be a list.")
  if (!is.character(files_path) || length(files_path) != 1 || !dir.exists(files_path)) {
    stop("'files_path' must be a single string pointing to an existing directory.")
  }
  
  # --- Load Configuration (YAML) ---
  config_path <- base::file.path(files_path, config_filename)
  config <- list() # Initialize config
  if (!base::file.exists(config_path)) stop("Configuration file not found: ", config_path)
  tryCatch({
    if (sink_active && verbose) message("Loading configuration from: ", config_path)
    config <- yaml::read_yaml(config_path) # Assign loaded config
    if (sink_active && verbose) message("  -> Loaded ", length(config), " top-level config parameters.")
  }, error = function(e) { stop("Error parsing YAML config file '", config_path, "': ", e$message, call. = FALSE) })
  
  # --- Validate Essential Config Sections ---
  required_config_sections <- c("required_data_columns", "preprocessing", "inclusion_criteria")
  missing_sections <- setdiff(required_config_sections, names(config))
  if(length(missing_sections)>0) stop("Config file missing required sections: ", paste(missing_sections, collapse=", "))
  if (!is.vector(config$required_data_columns) || length(config$required_data_columns) == 0) stop("Config 'required_data_columns' invalid.")
  expected_data_vars <- unique(config$required_data_columns)
  if (!is.list(config$preprocessing)) stop("Config 'preprocessing' section must be a list.")
  if (!is.list(config$inclusion_criteria)) stop("Config 'inclusion_criteria' section must be a list.")
  
  # --- Load SUV Dictionary ---
  # dict_suv_path <- base::file.path(files_path, "dict_suv.csv")
  dict_suv_path <- file.path(files_path, dict_suv_filename)
  dict_suv <- NULL
  if (!base::file.exists(dict_suv_path)) stop("SUV dictionary not found: ", dict_suv_path)
  tryCatch({
    if (sink_active && verbose) message("Loading SUV dictionary: ", dict_suv_path)
    dict_suv <- readr::read_csv(dict_suv_path, show_col_types = FALSE, progress = FALSE, locale = readr::locale(encoding = "UTF-8"))
    if (sink_active && verbose) message("  -> Success. Dictionary columns: ", paste(names(dict_suv), collapse=", "))
  }, error = function(e) { stop("Failed to load SUV dictionary '", dict_suv_path, "': ", e$message, call. = FALSE) })
  
  # --- Set Dictionary Column Defaults & Validate ---
  if (is.null(dict_suv_target_col)) {
    if (ncol(dict_suv) < 1) stop("SUV dictionary has no columns."); dict_suv_target_col <- names(dict_suv)[1]
    if (sink_active && verbose) message("Using first column ('", dict_suv_target_col, "') as target for SUV prep.")
  } else if (!dict_suv_target_col %in% names(dict_suv)) stop("Specified `dict_suv_target_col` not found.")
  if (is.null(dict_suv_potential_cols)) { dict_suv_potential_cols <- names(dict_suv) }
  else { missing_potential <- setdiff(dict_suv_potential_cols, names(dict_suv)); if(length(missing_potential)>0) stop("Specified `dict_suv_potential_cols` not found: ", paste(missing_potential, collapse=", "))}
  if (length(dict_suv_potential_cols) == 0) stop("No potential columns for SUV matching.")
  if (sink_active && verbose) message("Using columns for SUV matching: ", paste(dict_suv_potential_cols, collapse=", "))
  
  
  # --- Initialize Output Structures ---
  data_original <- loaded_data # Keep original for row index reference
  data_processed <- loaded_data # Start with a copy to modify
  validation_report <- list()  # For checks before filtering
  preprocessing_log <- list() # For filter counts
  
  # --- Extract Preprocessing Config ---
  preproc_params <- config$preprocessing # Already checked it's a list
  id_col <- preproc_params$id_column %||% "RID"
  min_entries <- preproc_params$min_entries_per_id %||% 1
  crit_source_file <- preproc_params$criteria_source_file %||% "data"
  inclusion_criteria_config <- config$inclusion_criteria # Already checked it's a list
  group_defs <- list( CU = list(DX_val = 0L, criteria = inclusion_criteria_config$CU %||% list()),
                      CI = list(DX_val = 1L, criteria = inclusion_criteria_config$CI %||% list()) )
  
  # --- Wrap the main processing loop ---
  if (sink_active && verbose) message("\n--- Starting Dataset Checks, Preparation & Preprocessing ---")
  else if (verbose) cat("\n--- Starting Dataset Checks, Preparation & Preprocessing ---\n")
  
  withCallingHandlers({
    for (dset_name in names(data_processed)) {
      
      if (sink_active && verbose) message("\nProcessing dataset: ", dset_name)
      else if(verbose) cat("\nProcessing dataset: ", dset_name, "\n")
      
      # Get pointers to original and modifiable lists/dataframes for this dataset
      current_dataset_list_orig <- data_original[[dset_name]]
      dataset_report <- list()
      dataset_preproc_log <- list(ids_start = NA_integer_, rows_start = NA_integer_, n_ids_removed_min_entries=0L) # Initialize specific log entry
      
      # =============================================================
      # Phase 1: Initial Checks and SUV Data Preparation
      # =============================================================
      if (sink_active && verbose) message("  Phase 1: Initial Checks & SUV Preparation")
      
      # --- 1a. Check source dataframes ---
      data_df_orig <- current_dataset_list_orig[["data"]]
      suv_df_orig <- current_dataset_list_orig[["data_suv_bi"]]
      can_do_checks <- TRUE
      can_process_suv <- TRUE
      
      if (!is.data.frame(data_df_orig)) {
        msg <- sprintf("Dataset '%s': 'data' not found or not a dataframe. Cannot perform checks or preprocessing.", dset_name); warning(msg, call.=FALSE)
        dataset_report$missing_data_df <- TRUE; can_do_checks <- FALSE; can_process_suv <- FALSE # Cannot process SUV if base data missing
      }
      if (can_do_checks && !is.data.frame(suv_df_orig)) {
        if("data_suv_bi" %in% names(current_dataset_list_orig)) msg <- sprintf("Dataset '%s': 'data_suv_bi' not a dataframe. Skipping SUV processing.", dset_name)
        else msg <- sprintf("Dataset '%s': 'data_suv_bi' not found. Skipping SUV processing.", dset_name)
        warning(msg, call.=FALSE); dataset_report$missing_suv_df <- TRUE; can_process_suv <- FALSE
      }
      if (can_process_suv && nrow(suv_df_orig) == 0) {
        warning(sprintf("Dataset '%s': 'data_suv_bi' is empty. Skipping SUV processing.", dset_name), call.=FALSE)
        dataset_report$empty_suv_df <- TRUE; can_process_suv <- FALSE
      }
      
      # # --- 1b. Check required 'data' columns ---
      # if(can_do_checks) {
      #   
      #   date_col_name <- scandate_column %||% "ScanDate"
      #   data_df_orig[[date_col_name]] <- "ScanDate"
      #   
      #   actual_data_vars <- names(data_df_orig)
      #   missing_data_vars <- setdiff(expected_data_vars, actual_data_vars)
      #   if (length(missing_data_vars) > 0) {
      #     msg <- sprintf("Dataset '%s' (data): Missing required columns (from config): %s", dset_name, paste(missing_data_vars, collapse = ", "))
      #     warning(msg, call. = FALSE); dataset_report$missing_data_columns <- missing_data_vars
      #   } else { if (sink_active && verbose) message("  -> 'data': Required columns present.") }
      # }
      

      
      # --- 1d. Process 'data_suv_bi' (Best Match, Rename/Subset) ---
      if (can_process_suv) {
        actual_suv_vars <- names(suv_df_orig)
        best_match_col_name <- NULL; best_match_score <- -1
        
        # Find best match
        for (potential_col in dict_suv_potential_cols) {
          potential_vars <- unique(stats::na.omit(dict_suv[[potential_col]]))
          potential_vars <- potential_vars[nzchar(potential_vars)]; if (length(potential_vars) == 0) next
          match_count <- length(intersect(actual_suv_vars, potential_vars))
          if (sink_active && verbose) message("    - Comparing vs dict col '", potential_col, "': ", match_count, " matches.")
          if (match_count > best_match_score) { best_match_score <- match_count; best_match_col_name <- potential_col }
        }
        
        if (is.null(best_match_col_name)) {
          msg <- sprintf("Dataset '%s' (data_suv_bi): Could not find suitable best matching column in dictionary '%s'.", dset_name, basename(dict_suv_path))
          warning(msg, call. = FALSE); dataset_report$suv_no_match <- TRUE
        } else {
          if (sink_active && verbose) message("  -> 'data_suv_bi': Best match column: '", best_match_col_name, "'.")
          dataset_report$suv_best_match_column <- best_match_col_name
          expected_suv_vars_from_match <- unique(stats::na.omit(dict_suv[[best_match_col_name]]))
          expected_suv_vars_from_match <- expected_suv_vars_from_match[nzchar(expected_suv_vars_from_match)]
          missing_suv_vars <- setdiff(expected_suv_vars_from_match, actual_suv_vars)
          if (length(missing_suv_vars) > 0) {
            msg <- sprintf("Dataset '%s' (data_suv_bi): Missing columns relative to best match '%s': %s", dset_name, best_match_col_name, paste(missing_suv_vars, collapse = ", "))
            warning(msg, call. = FALSE); dataset_report$missing_suv_columns <- missing_suv_vars
          } else { if (sink_active && verbose) message("  -> 'data_suv_bi': All columns from best match present.") }
          
          # Rename OR Subset (modifies data_processed directly)
          if (best_match_col_name == dict_suv_target_col) {
            if (sink_active && verbose) message("  -> 'data_suv_bi': Subsetting columns based on target ('", dict_suv_target_col, "').")
            vars_to_keep <- intersect(expected_suv_vars_from_match, actual_suv_vars)
            vars_to_discard <- setdiff(actual_suv_vars, expected_suv_vars_from_match)
            if (length(vars_to_keep) == 0) {
              msg <- sprintf("Dataset '%s' (data_suv_bi): Subsetting resulted in NO columns.", dset_name); warning(msg, call.=FALSE)
              dataset_report$suv_subset_empty <- TRUE
            }
            data_processed[[dset_name]]$data_suv_bi <- data_processed[[dset_name]]$data_suv_bi[, vars_to_keep, drop = FALSE]
            if(length(vars_to_discard) > 0){ dataset_report$discarded_suv_columns <- vars_to_discard; if(sink_active && verbose) message("    - Discarded: ", paste(vars_to_discard, collapse=", "))}
            else{ if(sink_active && verbose) message("    - No columns discarded.")}
          } else {
            if (sink_active && verbose) message("  -> 'data_suv_bi': Renaming based on map '", best_match_col_name, "' -> '", dict_suv_target_col, "'.")
            mapping_df <- unique(dict_suv[, c(best_match_col_name, dict_suv_target_col)])
            mapping_df <- mapping_df[!is.na(mapping_df[[best_match_col_name]]) & nzchar(mapping_df[[best_match_col_name]]), ]
            if (anyDuplicated(mapping_df[[best_match_col_name]])) { warning("Dataset '", dset_name, "': Duplicates in source column '", best_match_col_name, "' for renaming.", call.=FALSE); mapping_df <- mapping_df[!duplicated(mapping_df[[best_match_col_name]]), ] }
            mapping_df <- mapping_df[!is.na(mapping_df[[dict_suv_target_col]]) & nzchar(mapping_df[[dict_suv_target_col]]), ]
            renaming_map <- stats::setNames(mapping_df[[dict_suv_target_col]], mapping_df[[best_match_col_name]])
            current_names <- names(data_processed[[dset_name]]$data_suv_bi)
            names_to_rename_from <- intersect(current_names, names(renaming_map))
            if(length(names_to_rename_from)>0){
              new_names_for_these <- renaming_map[names_to_rename_from]; final_names <- current_names; match_indices <- match(names_to_rename_from, final_names)
              potential_new_set <- final_names; potential_new_set[match_indices] <- new_names_for_these
              if(anyDuplicated(potential_new_set)) { warning("Dataset '", dset_name, "': Renaming would create duplicate names. Skipping.", call. = FALSE); dataset_report$suv_rename_skipped_duplicates <- TRUE}
              else {
                final_names[match_indices] <- new_names_for_these; names(data_processed[[dset_name]]$data_suv_bi) <- final_names
                rename_log_pairs <- stats::setNames(new_names_for_these, names_to_rename_from); dataset_report$renamed_suv_columns <- rename_log_pairs
                if(sink_active && verbose){ rename_log_text <- paste(paste0("'", names(rename_log_pairs), "'->'", rename_log_pairs, "'"), collapse = "; "); message("    - Renamed: ", rename_log_text)}
              }
            } else { if(sink_active && verbose) message("    - No columns required renaming.") }
          } # End rename/subset
        } # End if best match found
      } # End if can_process_suv
      
      # =============================================================
      # Phase 2: Preprocessing Filters (Row/ID based)
      # =============================================================
      if (sink_active && verbose) message("  Phase 2: Preprocessing Filters")
      
      # --- Use the potentially modified data_processed as the source for filtering criteria ---
      crit_df_for_filtering <- data_processed[[dset_name]][[crit_source_file]]
      can_do_filtering <- TRUE # Assume true initially
      
      # Check if criteria source file is valid for filtering
      if (!is.data.frame(crit_df_for_filtering)) {
        warning(sprintf("Dataset '%s': Criteria source file '%s' invalid or missing after initial checks. Skipping row filters.", dset_name, crit_source_file), call.=FALSE)
        can_do_filtering <- FALSE
      } else if (!id_col %in% names(crit_df_for_filtering)) {
        warning(sprintf("Dataset '%s': ID column '%s' missing in '%s'. Skipping filters.", dset_name, id_col, crit_source_file), call.=FALSE)
        can_do_filtering <- FALSE
      }
      
      if (can_do_filtering) {
        ids_before_any_filter <- unique(crit_df_for_filtering[[id_col]])
        n_ids_start <- length(ids_before_any_filter)
        n_rows_start <- nrow(crit_df_for_filtering)
        dataset_preproc_log[[dset_name]]$ids_start <- n_ids_start
        dataset_preproc_log[[dset_name]]$rows_start <- n_rows_start
        if (verbose) message(sprintf("  - Initial for filtering: %d rows, %d unique IDs.", n_rows_start, n_ids_start))
        
        # --- Filter 1: Minimum Entries ---
        ids_after_min_entries <- ids_before_any_filter
        if (min_entries > 1) {
          entry_counts <- crit_df_for_filtering %>% dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) %>% dplyr::summarise(n_entries = dplyr::n(), .groups = "drop")
          ids_meeting_min_entries <- entry_counts %>% dplyr::filter(.data$n_entries >= min_entries) %>% dplyr::pull(!!rlang::sym(id_col))
          ids_removed_by_min_entries <- setdiff(ids_before_any_filter, ids_meeting_min_entries)
          n_ids_removed_min_entries <- length(ids_removed_by_min_entries)
          ids_after_min_entries <- ids_meeting_min_entries
          dataset_preproc_log[[dset_name]]$n_ids_removed_min_entries <- n_ids_removed_min_entries
          if (verbose) message(sprintf("  - Filter (min_entries >= %d): Removed %d IDs. Kept %d.", min_entries, n_ids_removed_min_entries, length(ids_after_min_entries)))
          if(length(ids_after_min_entries) == 0) can_do_filtering <- FALSE # Stop if no IDs left
        } else { if (verbose) message("  - Filter (min_entries <= 1): Step skipped.") }
        
        # --- Filter 2: Inclusion Criteria ---
        final_ids_to_keep <- ids_after_min_entries # Initialize with IDs passing min_entries
        if(can_do_filtering) { # Check if we still have IDs and can proceed
          can_do_group_filter <- "DX" %in% names(crit_df_for_filtering)
          if (!can_do_group_filter) { warning(sprintf("Dataset '%s': DX column missing in '%s'. Skipping group inclusion filters.", dset_name, crit_source_file), call.=FALSE) }
          else {
            ids_passed_inclusion_list <- list(); ids_processed_inclusion <- character(0)
            crit_df_after_min_entries <- crit_df_for_filtering %>% dplyr::filter(.data[[id_col]] %in% ids_after_min_entries) # Work with relevant subset
            
            for(group_name in names(group_defs)) {
              group_info <- group_defs[[group_name]]; dx_value <- group_info$DX_val; criteria <- group_info$criteria
              group_crit_df_initial <- crit_df_after_min_entries %>% dplyr::filter(.data$DX == dx_value)
              group_ids_initial <- unique(group_crit_df_initial[[id_col]]); ids_processed_inclusion <- union(ids_processed_inclusion, group_ids_initial)
              n_group_ids_initial <- length(group_ids_initial); dataset_preproc_log[[dset_name]][[paste0("n_ids_initial_", group_name)]] <- n_group_ids_initial
              if(n_group_ids_initial == 0) { if(verbose) message(sprintf("  - Filter (Inclusion %s): No individuals initially.", group_name)); ids_passed_inclusion_list[[group_name]] <- character(0); next }
              if (verbose) message(sprintf("  - Filter (Inclusion %s): Applying criteria to %d individuals...", group_name, n_group_ids_initial))
              
              available_criteria_vars <- intersect(names(criteria), names(group_crit_df_initial))
              missing_crit_cols <- setdiff(names(criteria), available_criteria_vars); if(length(missing_crit_cols)>0) warning(sprintf("Dataset '%s': Group '%s' missing columns for criteria: %s.", dset_name, group_name, paste(missing_crit_cols, collapse=", ")), call.=FALSE)
              
              if (length(available_criteria_vars) > 0) {
                filtered_group_df <- group_crit_df_initial
                for (crit_var in available_criteria_vars) {
                  crit_range <- criteria[[crit_var]]; if (!is.numeric(crit_range) || length(crit_range)!=2){ warning(sprintf("Dataset '%s': Invalid range for criterion '%s', group '%s'. Skipping.", dset_name, crit_var, group_name), call.=FALSE); next }
                  min_val <- crit_range[1]; max_val <- crit_range[2]
                  filtered_group_df <- filtered_group_df %>% dplyr::filter(dplyr::between(.data[[crit_var]], min_val, max_val))
                  if(verbose && sink_active) message(sprintf("    - Applied %s: %s [%.1f, %.1f]", group_name, crit_var, min_val, max_val))
                }
                ids_passed_this_group <- unique(filtered_group_df[[id_col]])
              } else { ids_passed_this_group <- group_ids_initial } # Keep all if no criteria defined/available
              
              ids_removed_this_group <- setdiff(group_ids_initial, ids_passed_this_group); n_ids_removed_this_group <- length(ids_removed_this_group)
              ids_passed_inclusion_list[[group_name]] <- ids_passed_this_group; dataset_preproc_log[[dset_name]][[paste0("n_ids_removed_inclusion_", group_name)]] <- n_ids_removed_this_group
              if (verbose) message(sprintf("    - Result (Inclusion %s): Removed %d IDs. Kept %d.", group_name, n_ids_removed_this_group, length(ids_passed_this_group)))
            } # End loop groups
            
            ids_passed_any_group_criteria <- unique(unlist(ids_passed_inclusion_list))
            ids_not_in_defined_groups <- setdiff(ids_after_min_entries, ids_processed_inclusion)
            if(length(ids_not_in_defined_groups) > 0 && verbose) message(sprintf("  - Note: %d individuals kept (not subject to group inclusion filters).", length(ids_not_in_defined_groups)))
            final_ids_to_keep <- union(ids_passed_any_group_criteria, ids_not_in_defined_groups)
          } # End if can_do_group_filter
        } # End if can_do_filtering (min entries didn't remove all)
        
        # --- Log Final ID count ---
        n_ids_final <- length(final_ids_to_keep); n_ids_removed_total_subjects <- n_ids_start - n_ids_final
        dataset_preproc_log[[dset_name]]$ids_final <- n_ids_final; dataset_preproc_log[[dset_name]]$ids_removed_total <- n_ids_removed_total_subjects
        if(verbose) message(sprintf("  - Final Filter Summary: Kept %d / %d unique IDs (Removed %d total).", n_ids_final, n_ids_start, n_ids_removed_total_subjects))
        
        # --- Apply Final Row Filter ---
        if (n_ids_final == 0) {
          warning(sprintf("Dataset '%s': No individuals remaining. Emptying data frames.", dset_name), call. = FALSE)
          for(df_name in names(data_processed[[dset_name]])) { if(is.data.frame(data_processed[[dset_name]][[df_name]])) data_processed[[dset_name]][[df_name]] <- data_processed[[dset_name]][[df_name]][0,] }
          dataset_preproc_log[[dset_name]]$rows_final <- 0; dataset_preproc_log[[dset_name]]$rows_removed_total <- n_rows_start
        } else {
          original_crit_df_for_indices <- data_original[[dset_name]][[crit_source_file]] # Use original row order
          original_row_count_source <- nrow(original_crit_df_for_indices)
          rows_to_keep_indices <- which(original_crit_df_for_indices[[id_col]] %in% final_ids_to_keep)
          n_rows_final <- length(rows_to_keep_indices); dataset_preproc_log[[dset_name]]$rows_final <- n_rows_final; dataset_preproc_log[[dset_name]]$rows_removed_total <- original_row_count_source - n_rows_final
          if(verbose) message(sprintf("  - Applying final row filter: Keeping %d / %d rows.", n_rows_final, original_row_count_source))
          
          # Apply indices to all DFs
          for (df_name in names(data_processed[[dset_name]])) {
            current_df_original <- data_original[[dset_name]][[df_name]]
            if(is.data.frame(current_df_original)) {
              rows_before <- nrow(current_df_original)
              if (rows_before == original_row_count_source) {
                data_processed[[dset_name]][[df_name]] <- current_df_original[rows_to_keep_indices, , drop = FALSE]
                rows_after <- nrow(data_processed[[dset_name]][[df_name]])
                dataset_preproc_log[[dset_name]][[paste0("rows_removed_", df_name)]] <- rows_before - rows_after
                if(verbose) message(sprintf("    - Filtered '%s': Kept %d / %d rows.", df_name, rows_after, rows_before))
              } else {
                warning(sprintf("Dataset '%s': Row count mismatch '%s' (%d) vs source '%s' (%d). Skipping row filter.", dset_name, df_name, rows_before, crit_source_file, original_row_count_source), call. = FALSE)
                dataset_preproc_log[[dset_name]][[paste0("rows_removed_", df_name)]] <- 0 # No rows removed
              }
            } # else keep non-df elements
          } # End loop applying row filter
        } # End if n_ids_final > 0
        
        
        # --- Step 4: Calculate and Add UV data (if possible) ---
        # Check if both SUV and VOL data exist and are valid after filtering
        if ("data_suv_bi" %in% names(data_processed[[dset_name]]) &&
            "data_vol_bi" %in% names(data_processed[[dset_name]]) && # Check if data_vol_bi exists
            is.data.frame(data_processed[[dset_name]]$data_suv_bi) &&
            is.data.frame(data_processed[[dset_name]]$data_vol_bi) &&
            nrow(data_processed[[dset_name]]$data_suv_bi) == nrow(data_processed[[dset_name]]$data_vol_bi) && # Ensure same rows after filtering
            nrow(data_processed[[dset_name]]$data_suv_bi) > 0) { # Ensure not empty
          
          if(verbose) message(sprintf("  - Calculating UV data (SUV * Volume)..."))
          data_suv_filt <- data_processed[[dset_name]]$data_suv_bi
          data_vol_filt <- data_processed[[dset_name]]$data_vol_bi
          id_col_present_suv <- id_col %in% names(data_suv_filt)
          id_col_present_vol <- id_col %in% names(data_vol_filt)
          
          # Identify numeric columns (excluding ID) in both dataframes
          numeric_suv_cols <- names(data_suv_filt)[sapply(data_suv_filt, is.numeric)]
          numeric_vol_cols <- names(data_vol_filt)[sapply(data_vol_filt, is.numeric)]
          if(id_col_present_suv) numeric_suv_cols <- setdiff(numeric_suv_cols, id_col)
          if(id_col_present_vol) numeric_vol_cols <- setdiff(numeric_vol_cols, id_col)
          
          # Find common numeric columns to multiply
          common_numeric_cols <- intersect(numeric_suv_cols, numeric_vol_cols)
          
          if (length(common_numeric_cols) > 0) {
            # Perform element-wise multiplication only on common numeric columns
            # Initialize data_uv with non-numeric columns if needed (like ID)
            if(id_col_present_suv && id_col_present_vol) {
              data_uv <- data.frame(placeholder_id = data_suv_filt[[id_col]])
              names(data_uv)[1] <- id_col
            } else {
              data_uv <- data.frame(matrix(ncol = 0, nrow = nrow(data_suv_filt))) # Empty df with correct rows
            }
            
            # Multiply common columns
            data_uv[common_numeric_cols] <- data_suv_filt[, common_numeric_cols, drop=FALSE] * data_vol_filt[, common_numeric_cols, drop=FALSE]
            
            # Add the calculated data_uv to the list for this dataset
            data_processed[[dset_name]]$data_uv_bi <- data_uv
            if(verbose) message(sprintf("    -> Added 'data_uv_bi' data frame with %d columns.", ncol(data_uv)))
            
          } else {
            warning(sprintf("Dataset '%s': No common numeric columns found between 'data_suv_bi' and 'data_vol_bi' after filtering. Cannot calculate 'data_uv'.", dset_name), call. = FALSE)
            # Optionally add an empty data_uv frame? Or leave it NULL? Let's leave it NULL.
            data_processed[[dset_name]]$data_uv <- NULL
          }
        } else {
          if(verbose) message(sprintf("  - Skipping UV data calculation: 'data_suv_bi' or 'data_vol_bi' missing, empty, non-dataframe, or have mismatched rows after filtering."))
          # Ensure data_uv element doesn't exist if calculation skipped
          data_processed[[dset_name]]$data_uv <- NULL
        }
        
        # --- Step 5: Ensure Time Variable Exists ---
        if(verbose) message(sprintf("  - Processing time variable..."))
        time_calc_success <- FALSE
        
        # Get config params for time
        date_col_name <- preproc_params$date_column %||% "ScanDate"
        rel_time_col_name <- preproc_params$relative_time_column %||% "years.from.baseline"
        time_col_name <- preproc_params$time_output_column %||% "time" # Final desired name
        center_time <- isTRUE(preproc_params$center_time_variable %||% TRUE) # Default to TRUE
        date_format_arg <- preproc_params$date_format # NULL if not specified
        
        crit_df_ref <- data_processed[[dset_name]][[crit_source_file]] # Reference the df in the list
        time_added_or_found <- FALSE
        calculated_time_vector <- NULL
        
        # --- Attempt 1: Calculate from Date Column ---
        if (date_col_name %in% names(crit_df_ref)) {
          if(verbose) message(sprintf("    - Found date column '%s'. Attempting calculation.", date_col_name))
          temp_dates <- crit_df_ref[[date_col_name]]
          converted_dates <- tryCatch({
            if (is.null(date_format_arg)) as.Date(temp_dates)
            else as.Date(temp_dates, format = date_format_arg)
          }, error = function(e) { NULL }) # Return NULL on error
          
          if (!is.null(converted_dates) && !all(is.na(converted_dates))) {
            # Use mean date per ID for baseline/centering
            mean_baseline_dates <- ave(converted_dates, crit_df_ref[[id_col]], FUN = function(x) mean(x, na.rm = TRUE))
            time_diff_days <- as.numeric(converted_dates - mean_baseline_dates)
            days_in_year <- 365.2425
            calculated_time_vector <- time_diff_days / days_in_year
            time_added_or_found <- TRUE
            if(verbose) message(sprintf("      -> Calculated time from '%s'.", date_col_name))
          } else {
            warning(sprintf("Dataset '%s': Found date column '%s' but failed to convert to Date objects. Will check for '%s'.",
                            dset_name, date_col_name, rel_time_col_name), call.=FALSE)
          }
        } else {
          if(verbose) message(sprintf("    - Date column '%s' not found. Checking for '%s'.", date_col_name, rel_time_col_name))
        }
        
        # --- Attempt 2: Use Relative Time Column (if date method failed/skipped) ---
        if (!time_added_or_found) {
          if (rel_time_col_name %in% names(crit_df_ref)) {
            if(verbose) message(sprintf("    - Found relative time column '%s'. Using it.", rel_time_col_name))
            temp_rel_time <- crit_df_ref[[rel_time_col_name]]
            if (!is.numeric(temp_rel_time)) {
              warning(sprintf("Dataset '%s': Relative time column '%s' found but is not numeric. Cannot use it.", dset_name, rel_time_col_name), call.=FALSE)
            } else if (all(is.na(temp_rel_time))) {
              warning(sprintf("Dataset '%s': Relative time column '%s' contains only NAs. Cannot use it.", dset_name, rel_time_col_name), call.=FALSE)
            } else {
              calculated_time_vector <- temp_rel_time # Use the column directly first
              time_added_or_found <- TRUE
              if(verbose) message(sprintf("      -> Using values from '%s'.", rel_time_col_name))
            }
          } else {
            if(verbose) message(sprintf("    - Relative time column '%s' not found.", rel_time_col_name))
          }
        }
        
        # --- Final Check and Processing ---
        if (!time_added_or_found || is.null(calculated_time_vector)) {
          # Stop if neither source column was found or usable
          stop(sprintf("Dataset '%s': Required time information could not be found or calculated using specified columns ('%s' or '%s'). Please check data and configuration.",
                       dset_name, date_col_name, rel_time_col_name))
        } else {
          # --- Optional Centering (Recommended) ---
          if(center_time) {
            # Center around the mean *per subject* even if input was relative time
            if(!id_col %in% names(crit_df_ref)){
              warning(sprintf("Dataset '%s': Cannot center time variable because ID column '%s' is missing from '%s'. Using uncentered time.", dset_name, id_col, crit_source_file), call.=FALSE)
            } else {
              mean_time_per_id <- ave(calculated_time_vector, crit_df_ref[[id_col]], FUN = function(x) mean(x, na.rm=TRUE))
              # Replace NA means (e.g., for single-visit subjects after filtering) with 0? Or keep NA? Let's keep NA for time if mean is NA.
              # mean_time_per_id[is.na(mean_time_per_id)] <- 0 # Optional: center single visits at 0
              calculated_time_vector <- calculated_time_vector - mean_time_per_id
              if(verbose) message("      -> Centered time variable around mean per subject.")
              time_calc_success <- TRUE
            }
          } else {
            if(verbose) message("      -> Using time variable without centering per subject.")
          }
          
          # --- Assign final time column ---
          # Add/overwrite the column directly in the data_processed list
          data_processed[[dset_name]][[crit_source_file]][[time_col_name]] <- calculated_time_vector
          if(verbose) message(sprintf("    -> Final time variable stored in column '%s'.", time_col_name))
          
          # Update crit_df reference if needed for subsequent steps within this loop iteration
          # (Though filtering is done, this ensures consistency if crit_df is used later)
          crit_df <- data_processed[[dset_name]][[crit_source_file]]
        }
        # --- End Time Variable Processing ---
        
        # Optional: Stop processing dataset if time calculation failed and time is essential?
        if (!time_calc_success) { stop("time variable creation failed.") }
        
        
      } # End if can_do_filtering
      
      # --- Store validation report for this dataset ---
      validation_report[[dset_name]] <- dataset_report
      # Store preprocessing log for this dataset
      preprocessing_log[[dset_name]] <- dataset_preproc_log
      
    } # End loop through datasets
  }, warning = warning_handler) # End of withCallingHandlers Block
  
  # --- Final message ---
  if (sink_active && verbose) message("\n--- Finished Data Preprocessing ---")
  else if(verbose) cat("\n--- Finished Data Preprocessing ---\n")
  
  # --- Return Results ---
  return(list(
    data = data_processed,
    config = config,
    validation_report = validation_report,
    preprocessing_log = preprocessing_log,
    log_file = actual_log_path
  ))
}