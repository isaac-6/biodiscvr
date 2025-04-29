# Place this in R/preprocess_datasets.R

#' Preprocess Loaded Datasets Based on Configuration Criteria
#'
#' Filters datasets based on minimum entries per ID and group-specific
#' inclusion criteria (e.g., AGE, MMSE, GCDR ranges) defined in the
#' configuration file.
#'
#' @param checked_data List. The output from `check_and_prepare_data`, which
#'   must contain `$data` (the list of loaded and potentially renamed datasets)
#'   and `$config` (the loaded configuration).
#' @param verbose Logical. If TRUE, print messages about the filtering process
#'   and counts of removed individuals/rows. Defaults to TRUE.
#'
#' @return A list structure identical to `checked_data$data`, but with data frames
#'   filtered according to the preprocessing rules specified in the configuration.
#'   Individuals not meeting the minimum entry count or group-specific inclusion
#'   criteria are removed.
#'
#' @details
#' This function performs the following preprocessing steps for each dataset:
#' 1.  Reads preprocessing parameters from `checked_data$config$preprocessing`
#'     (e.g., `min_entries_per_id`, `id_column`, `criteria_source_file`).
#' 2.  Reads group-specific inclusion criteria (e.g., `AGE`, `MMSE`, `GCDR` ranges)
#'     from `checked_data$config$power_params` under keys like `CU` and `CI`.
#' 3.  Filters individuals based on `min_entries_per_id` using the data in the
#'     `criteria_source_file` (typically 'data').
#' 4.  Applies inclusion criteria range filters separately for CU (DX=0) and
#'     CI (DX=1) groups, using the data in `criteria_source_file`.
#'     - Issues warnings if required criteria columns (e.g., AGE, MMSE, GCDR)
#'       are missing in a dataset's `criteria_source_file`, skipping the
#'       filter for that specific criterion in that dataset.
#' 5.  Determines the final set of individuals (IDs) who pass both minimum
#'     entries and their respective group's inclusion criteria.
#' 6.  Filters *all* data frames within the dataset (e.g., 'data', 'data_suv_bi')
#'     to retain only data corresponding to the final set of included individuals.
#'
#' @export
#' @importFrom dplyr group_by summarise filter n select all_of pull between distinct
#' @importFrom rlang .data `%||%`
preprocess_datasets <- function(checked_data, verbose = TRUE) {
  
  # --- Input Validation ---
  if (!is.list(checked_data) || !all(c("data", "config") %in% names(checked_data))) {
    stop("'checked_data' must be a list containing 'data' and 'config' (output of check_and_prepare_data).")
  }
  config <- checked_data$config
  if (!is.list(config)) stop("'checked_data$config' must be a list.")
  
  # --- Extract Configuration ---
  # Use rlang::`%||%` for safe defaults if parameters are missing
  preproc_params <- config$preprocessing %||% list()
  id_col <- preproc_params$id_column %||% "RID"
  min_entries <- preproc_params$min_entries_per_id %||% 1 # Default to 1 if not specified
  crit_source_file <- preproc_params$criteria_source_file %||% "data"
  
  power_params_config <- config$inclusion_criteria %||% list()
  group_defs <- list( # Define groups and their criteria locations
    CU = list(DX_val = 0L, criteria = power_params_config$CU %||% list()),
    CI = list(DX_val = 1L, criteria = power_params_config$CI %||% list())
    # Add more groups here if defined in config
  )
  
  # --- Initialize Output Data Structure ---
  data_original <- checked_data$data
  data_processed <- data_original # Start with a copy to modify
  
  if (verbose) message("--- Starting Dataset Preprocessing ---")
  
  # --- Loop Through Each Dataset ---
  for (dset_name in names(data_processed)) {
    if (verbose) message("\nProcessing dataset: ", dset_name)
    current_dataset_list <- data_processed[[dset_name]]
    
    # Check if the file designated for criteria exists
    if (!crit_source_file %in% names(current_dataset_list) || is.null(current_dataset_list[[crit_source_file]])) {
      warning(sprintf("Dataset '%s': Criteria source file '%s' not found or is NULL. Skipping preprocessing for this dataset.",
                      dset_name, crit_source_file), call. = FALSE)
      next # Skip to next dataset
    }
    crit_df <- current_dataset_list[[crit_source_file]]
    
    # Check for essential columns in the criteria data frame
    if (!id_col %in% names(crit_df)) {
      warning(sprintf("Dataset '%s': ID column '%s' not found in criteria source file '%s'. Skipping preprocessing.",
                      dset_name, id_col, crit_source_file), call.=FALSE)
      next
    }
    if (!"DX" %in% names(crit_df)) {
      warning(sprintf("Dataset '%s': DX column not found in criteria source file '%s'. Skipping group-specific inclusion criteria filtering.",
                      dset_name, crit_source_file), call.=FALSE)
      # Allow min_entries filter to proceed, but skip group criteria
      can_do_group_filter <- FALSE
    } else {
      can_do_group_filter <- TRUE
    }
    
    
    n_rows_start <- nrow(crit_df)
    n_ids_start <- length(unique(crit_df[[id_col]]))
    if (verbose) message(sprintf("  - Initial state: %d rows, %d unique IDs.", n_rows_start, n_ids_start))
    
    # --- Step 1: Filter by Minimum Entries ---
    if (min_entries > 1) {
      entry_counts <- crit_df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
        dplyr::summarise(n_entries = dplyr::n(), .groups = "drop")
      
      ids_meeting_min_entries <- entry_counts |>
        dplyr::filter(.data$n_entries >= min_entries) |>
        dplyr::pull(!!rlang::sym(id_col)) # Use pull with !!sym() for robustness
      
      n_ids_after_min_entries <- length(ids_meeting_min_entries)
      n_ids_removed_min_entries <- n_ids_start - n_ids_after_min_entries
      
      if (verbose) {
        message(sprintf("  - Filter (min_entries >= %d): Kept %d / %d IDs (Removed %d).",
                        min_entries, n_ids_after_min_entries, n_ids_start, n_ids_removed_min_entries))
      }
      
      if (n_ids_after_min_entries == 0) {
        warning(sprintf("Dataset '%s': No individuals remaining after min_entries filter. Skipping further processing.", dset_name), call. = FALSE)
        # Empty out all data frames for this dataset
        for(df_name in names(current_dataset_list)) {
          if(is.data.frame(current_dataset_list[[df_name]])) {
            data_processed[[dset_name]][[df_name]] <- current_dataset_list[[df_name]][0, , drop = FALSE]
          }
        }
        next # Skip to next dataset
      }
      # Filter the criteria dataframe immediately to simplify next steps
      crit_df <- crit_df |>
        dplyr::filter(.data[[id_col]] %in% ids_meeting_min_entries)
      
    } else {
      ids_meeting_min_entries <- unique(crit_df[[id_col]]) # All IDs included if min_entries <= 1
      n_ids_after_min_entries <- n_ids_start
      if (verbose) message("  - Filter (min_entries <= 1): Step skipped, all IDs kept.")
    }
    
    
    # --- Step 2: Filter by Group-Specific Inclusion Criteria ---
    ids_passed_inclusion <- c() # Vector to store IDs passing group criteria
    
    if(can_do_group_filter) {
      initial_ids_for_inclusion <- unique(crit_df[[id_col]]) # IDs remaining after min_entries filter
      ids_processed_inclusion <- character(0) # Track IDs evaluated
      ids_passed_inclusion_list <- list() # Store passing IDs per group
      
      for(group_name in names(group_defs)) {
        group_info <- group_defs[[group_name]]
        dx_value <- group_info$DX_val
        criteria <- group_info$criteria
        required_criteria_vars <- names(criteria) # e.g., c("AGE", "MMSE", "GCDR")
        
        # Filter data for the current group being evaluated
        group_crit_df <- crit_df |> dplyr::filter(.data$DX == dx_value)
        group_ids_initial <- unique(group_crit_df[[id_col]])
        ids_processed_inclusion <- union(ids_processed_inclusion, group_ids_initial) # Mark these IDs as processed
        
        n_group_ids_initial <- length(group_ids_initial)
        if (n_group_ids_initial == 0) {
          if(verbose) message(sprintf("  - Filter (Inclusion %s, DX=%d): No individuals in this group.", group_name, dx_value))
          ids_passed_inclusion_list[[group_name]] <- character(0)
          next # Skip to next group
        }
        
        if (verbose) message(sprintf("  - Filter (Inclusion %s, DX=%d): Applying criteria to %d individuals...", group_name, dx_value, n_group_ids_initial))
        
        # Check for missing criteria columns and build filter expression
        valid_filters <- list()
        missing_crit_cols <- setdiff(required_criteria_vars, names(group_crit_df))
        if (length(missing_crit_cols) > 0) {
          warning(sprintf("Dataset '%s': Group '%s' is missing inclusion criteria columns in '%s': %s. Filtering only on available criteria.",
                          dset_name, group_name, crit_source_file, paste(missing_crit_cols, collapse=", ")), call. = FALSE)
        }
        
        available_criteria_vars <- intersect(required_criteria_vars, names(group_crit_df))
        
        if (length(available_criteria_vars) == 0) {
          if(verbose) message("    - No criteria columns available for this group. Keeping all individuals in this group.")
          ids_passed_this_group <- group_ids_initial # Keep all if no criteria cols exist
        } else {
          # Filter the group data frame based on available criteria ranges
          filtered_group_df <- group_crit_df
          for (crit_var in available_criteria_vars) {
            crit_range <- criteria[[crit_var]]
            if (!is.numeric(crit_range) || length(crit_range) != 2) {
              warning(sprintf("Dataset '%s': Invalid range definition for criterion '%s' in group '%s'. Skipping this filter.",
                              dset_name, crit_var, group_name), call.=FALSE)
              next # Skip this specific criterion filter
            }
            min_val <- crit_range[1]
            max_val <- crit_range[2]
            
            # Apply filter using dplyr::between or base R
            filtered_group_df <- filtered_group_df |>
              dplyr::filter(dplyr::between(.data[[crit_var]], min_val, max_val))
            
            if (verbose) message(sprintf("    - Applied %s filter: %s >= %.2f & %s <= %.2f", group_name, crit_var, min_val, crit_var, max_val))
            
          }
          # Get the IDs that remained after all applicable filters for this group
          ids_passed_this_group <- unique(filtered_group_df[[id_col]])
        }
        
        ids_passed_inclusion_list[[group_name]] <- ids_passed_this_group
        n_group_ids_passed <- length(ids_passed_this_group)
        n_group_ids_removed <- n_group_ids_initial - n_group_ids_passed
        if (verbose) {
          message(sprintf("    - Result %s: Kept %d / %d individuals (Removed %d).",
                          group_name, n_group_ids_passed, n_group_ids_initial, n_group_ids_removed))
        }
        
      } # End loop over groups (CU, CI)
      
      # Combine IDs from all groups that passed their respective criteria
      ids_passed_inclusion <- unique(unlist(ids_passed_inclusion_list))
      
      # What about individuals who were present after min_entries but had DX values not in group_defs?
      # We should arguably keep them if they passed min_entries, as they weren't subject to inclusion criteria.
      # Find IDs that were eligible but not processed (e.g., DX = NA or other value)
      ids_not_processed_by_group <- setdiff(initial_ids_for_inclusion, ids_processed_inclusion)
      if(length(ids_not_processed_by_group) > 0 && verbose) {
        message(sprintf("  - Note: %d individuals were not assigned to CU/CI groups (e.g., missing DX) and were not subject to inclusion criteria.",
                        length(ids_not_processed_by_group)))
      }
      # Final set includes those passing group criteria + those not processed by group criteria
      final_ids_to_keep <- union(ids_passed_inclusion, ids_not_processed_by_group)
      
    } else { # Cannot do group filter (DX missing)
      final_ids_to_keep <- ids_meeting_min_entries # Keep all who passed min_entries
      if(verbose) message("  - Skipping group-specific inclusion criteria (DX column missing).")
    }
    
    
    # --- Step 3: Determine Row Indices to Keep and Apply to All Data Frames ---
    n_ids_final <- length(final_ids_to_keep)
    n_ids_removed_total <- n_ids_start - n_ids_final
    if (verbose) message(sprintf("  - Final filter: Keeping %d / %d total unique IDs (Total removed: %d).",
                                 n_ids_final, n_ids_start, n_ids_removed_total))
    
    if (n_ids_final == 0) {
      warning(sprintf("Dataset '%s': No individuals remaining after all preprocessing steps. Emptying data frames.", dset_name), call. = FALSE)
      # Empty out all data frames for this dataset
      for(df_name in names(current_dataset_list)) {
        if(is.data.frame(current_dataset_list[[df_name]])) {
          data_processed[[dset_name]][[df_name]] <- current_dataset_list[[df_name]][0, , drop = FALSE]
        }
      }
      next # Skip to next dataset
    }
    
    # Get the original criteria data frame to find row indices
    original_crit_df <- data_original[[dset_name]][[crit_source_file]]
    original_row_count_source <- nrow(original_crit_df)
    
    # Find the indices of rows in the *original* source dataframe that belong to the final IDs
    # This assumes the filtering was based on individuals, not specific rows within individuals yet.
    # If inclusion criteria filtering *removed specific rows* from crit_df earlier,
    # we need to re-filter the original based on final_ids_to_keep to get the right indices.
    # Let's recalculate the rows to keep based on the final ID list applied to the original source df.
    rows_to_keep_indices <- which(original_crit_df[[id_col]] %in% final_ids_to_keep)
    n_rows_final <- length(rows_to_keep_indices)
    
    if (verbose) message(sprintf("  - Applying row filter: Keeping %d / %d rows based on final IDs. Applying to all data frames...",
                                 n_rows_final, original_row_count_source))
    
    
    # Iterate through all data frames within the current dataset and apply row subsetting
    for (df_name in names(current_dataset_list)) {
      current_df_original <- data_original[[dset_name]][[df_name]] # Get the original df
      
      if(is.data.frame(current_df_original)) {
        # Check if the row count matches the original source dataframe
        if (nrow(current_df_original) == original_row_count_source) {
          # Apply the filter using the row indices identified from the source df
          filtered_df <- current_df_original[rows_to_keep_indices, , drop = FALSE]
          rows_removed <- original_row_count_source - nrow(filtered_df) # Should match n_rows_start - n_rows_final
          
          # Update the data frame in the list being built
          data_processed[[dset_name]][[df_name]] <- filtered_df
          
          if (verbose && rows_removed > 0) {
            message(sprintf("    - Filtered '%s': Kept %d / %d rows.", df_name, nrow(filtered_df), original_row_count_source ))
          } else if (verbose) {
            message(sprintf("    - Filtered '%s': Kept all %d rows (matched source row count).", df_name, nrow(filtered_df)))
          }
        } else {
          # Row counts don't match - cannot safely apply index filter
          warning(sprintf("Dataset '%s': Data frame '%s' has %d rows, but criteria source file '%s' has %d rows. Cannot apply row-index filter reliably. Keeping '%s' unfiltered.",
                          dset_name, df_name, nrow(current_df_original), crit_source_file, original_row_count_source, df_name), call. = FALSE)
          # Keep the original dataframe in the processed list if row counts don't match
          data_processed[[dset_name]][[df_name]] <- current_df_original
        }
      } else {
        # Keep non-data-frame elements as they are
        data_processed[[dset_name]][[df_name]] <- current_df_original
      }
    } # End loop filtering data frames within dataset
    
  } # End loop through datasets
  
  data_processed$config <- config
  
  if (verbose) message("\n--- Finished Dataset Preprocessing ---")
  
  # Return only the $data part, now filtered
  return(data_processed)
}
