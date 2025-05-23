#' Run Iterative Regional Ablation on Discovered Biomarkers
#'
#' Reads a CSV file of discovered biomarkers (typically from `run_experiments`),
#' and for each biomarker, performs an iterative single-region ablation using
#' `region_ablation` to minimize aggregated Sample Size Estimate (SSE).
#' Results of the ablation process are saved to a new CSV file.
#'
#' @param discovered_biomarkers_csv_path Character string. Path to the input CSV
#'   file containing discovered biomarkers. Expected columns include at least:
#'   `regs_numerator`, `regs_denominator`, `var_composition`, and identifiers
#'   like `discovery_dataset`, `experiment_tag`, `group_evaluated`.
#' @param prepared_data_list List. The output from `preprocess_data()`,
#'   containing filtered data for multiple datasets.
#' @param config List. The loaded main configuration object from
#'   `preprocess_data$config`.
#' @param output_ablation_results_csv_path Character string. Full path to save the
#'   ablation results as a *new* CSV file. If the file exists, it will be
#'   overwritten. Directory created if needed.
#' @param datasets_for_ablation_eval Character vector. Names of the datasets within
#'   `prepared_data_list` to use for evaluating SSE *during the ablation process*.
#'   Defaults to using all datasets present in `prepared_data_list`.
#' @param groups_for_ablation_eval Character vector. Groups ("CU", "CI") to evaluate
#'   within each dataset *during the ablation process*. Defaults to `c("CU", "CI")`.
#' @param features_for_ablation Character vector or NULL. A specific list of feature
#'   names to consider for ablation. If NULL (default), common features across
#'   `datasets_for_ablation_eval` will be identified and used. The regions from
#'   the input CSV must be a subset of these.
#' @param verbose Logical. Print progress messages? Defaults to TRUE.
#'
#' @return A data frame containing the aggregated results of the ablation process
#'   for all processed biomarkers, which is also saved to
#'   `output_ablation_results_csv_path`. Each row represents one initial biomarker
#'   and includes its original definition, initial and final aggregated SSE,
#'   the refined regions, and a summary of the ablation process (e.g., number
#'   of regions removed). Returns `NULL` if the input CSV cannot be read or no
#'   biomarkers are processed.
#'
#' @details
#' Workflow:
#' 1. Reads the `discovered_biomarkers_csv_path`.
#' 2. Identifies common features across `datasets_for_ablation_eval` if
#'    `features_for_ablation` is not provided.
#' 3. Iterates through each row (discovered biomarker) in the input CSV.
#' 4. For each biomarker:
#'    a. Validates that its defined numerator and denominator regions are present
#'       within the set of features available for ablation.
#'    b. Calls `region_ablation` with the biomarker's definition,
#'       the specified evaluation datasets/groups, and other parameters.
#'    c. Collects key results from the ablation: initial and final aggregated SSE,
#'       initial and final region counts, and number of regions removed.
#' 5. Combines all results into a single data frame.
#' 6. Saves the aggregated ablation results table to `output_ablation_results_csv_path`.
#'
#' This function is intended to assess the stability and identify core components
#' of promising biomarkers found in earlier discovery stages.
#'
#' @export
#' @importFrom utils read.csv
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr bind_rows select all_of filter n_distinct between
#' @importFrom stringr str_split
#' @importFrom rlang `%||%` .data sym is_scalar_character is_scalar_logical is_scalar_integerish
#' @importFrom stats setNames median quantile IQR sd na.omit aggregate as.formula
run_ablation <- function(
    discovered_biomarkers_csv_path,
    prepared_data_list,
    config,
    output_ablation_results_csv_path,
    datasets_for_ablation_eval = NULL,
    groups_for_ablation_eval = NULL,
    features_for_ablation = NULL,
    verbose = TRUE) {
  
  # --- Input Validation ---
  stopifnot(
    rlang::is_scalar_character(discovered_biomarkers_csv_path),
    file.exists(discovered_biomarkers_csv_path),
    is.list(prepared_data_list), length(prepared_data_list) > 0,
    is.list(config),
    rlang::is_scalar_character(output_ablation_results_csv_path),
    is.null(datasets_for_ablation_eval) || (is.character(datasets_for_ablation_eval) && length(datasets_for_ablation_eval) > 0),
    is.null(groups_for_ablation_eval) || (is.character(groups_for_ablation_eval) && all(groups_for_ablation_eval %in% c("CU", "CI"))),
    is.null(features_for_ablation) || is.character(features_for_ablation),
    rlang::is_scalar_logical(verbose)
  )
  missing_dsets <- setdiff(datasets_for_ablation_eval, names(prepared_data_list))
  if (length(missing_dsets) > 0) stop("Datasets for ablation not found in prepared_data_list: ", paste(missing_dsets, collapse=", "))
  id_col <- config$preprocessing$id_column %||% "RID"
  
  
  # --- Load Discovered Biomarkers ---
  if(verbose) message("--- Loading Discovered Biomarkers for Ablation ---")
  discovery_df <- try(readr::read_csv(discovered_biomarkers_csv_path, show_col_types = FALSE), silent = TRUE)
  if (inherits(discovery_df, "try-error")) {
    stop("Failed to read discovered biomarkers CSV: ", discovered_biomarkers_csv_path, "\nError: ", conditionMessage(attr(discovery_df,"condition")))
  }
  if (nrow(discovery_df) == 0) {
    warning("Discovered biomarkers CSV is empty: ", discovered_biomarkers_csv_path, call. = FALSE); return(NULL)
  }
  message(sprintf("Loaded %d discovered biomarkers to process for ablation.", nrow(discovery_df)))
  
  # Check for essential columns from discovery results
  required_discovery_cols <- c("regs_numerator", "regs_denominator", "var_composition") # Minimum needed
  # Add other columns you want to carry over to the output, e.g., experiment_tag, discovery_dataset, group_evaluated, fitness_value
  carry_over_cols <- c("experiment_tag", "discovery_dataset", "group_evaluated", "fitness_value", "bilateral") # Example
  carry_over_cols <- intersect(names(discovery_df), c(carry_over_cols, "evaluation_dataset", "evaluation_group")) # in case evaluation csv is given
  all_needed_cols <- unique(c(required_discovery_cols, carry_over_cols))
  
  missing_discovery_cols <- setdiff(all_needed_cols, names(discovery_df))
  if (length(missing_discovery_cols) > 0) {
    stop("Input discovery CSV is missing required/carry-over columns: ", paste(missing_discovery_cols, collapse=", "))
  }
  
  
  # --- Initialize List to Store All Ablation Results ---
  all_ablation_summary_results <- list()
  result_counter <- 0
  
  # --- Loop Through Each Discovered Biomarker ---
  if(verbose) message("\n--- Starting Ablation Process for Each Discovered Biomarker ---")
  for (i in seq_len(nrow(discovery_df))) {
    
    biomarker_row <- discovery_df[i, , drop = FALSE]
    # Create a unique ID for this original biomarker for clarity
    original_biomarker_id_str <- paste(
      biomarker_row$discovery_dataset %||% "UnknownDset",
      biomarker_row$experiment_tag %||% "NoTag",
      biomarker_row$group_evaluated %||% "NoGroup",
      "Row", i, sep="_"
    )
    if(verbose) message(sprintf("\nProcessing Biomarker %d/%d (ID: %s)...", i, nrow(discovery_df), original_biomarker_id_str))
    
    # Extract definition
    regs_num_str <- biomarker_row$regs_numerator[[1]] %||% ""
    regs_den_str <- biomarker_row$regs_denominator[[1]] %||% ""
    var_comp <- biomarker_row$var_composition[[1]]
    
    if (!nzchar(regs_num_str) || !nzchar(regs_den_str) || is.na(var_comp)) {
      warning(sprintf("Skipping biomarker ID '%s': Missing/invalid definition.", original_biomarker_id_str), call. = FALSE)
      next
    }
    initial_num_regions <- stringr::str_split(regs_num_str, ", ")[[1]]
    initial_den_regions <- stringr::str_split(regs_den_str, ", ")[[1]]
    initial_num_regions <- initial_num_regions[nzchar(initial_num_regions)]
    initial_den_regions <- initial_den_regions[nzchar(initial_den_regions)]
    
    if (length(initial_num_regions) == 0 || length(initial_den_regions) == 0) {
      warning(sprintf("Skipping biomarker ID '%s': Parsed definition resulted in empty numerator or denominator.", original_biomarker_id_str), call. = FALSE)
      next
    }
    
    
    
    # group to evaluate. If NULL, takes the group from the data row.
    if (is.null(groups_for_ablation_eval)) {
      if("evaluation_group" %in% names(biomarker_row)){
        # this is the csv from run_evaluation
        eval_groups <- biomarker_row$evaluation_group
      } else {
        # this is the csv from run_experiments
        eval_groups <- biomarker_row$group_evaluated
      }
      
    } else {
      eval_groups <- groups_for_ablation_eval
    }
    
    # cohorts/datasets to evaluate. If NULL, takes the one(s) used for discover from the data row.
    if (is.null(datasets_for_ablation_eval)) {
      if("evaluation_dataset" %in% names(biomarker_row)){
        # this is the csv from run_evaluation
        eval_dsets = strsplit(biomarker_row$evaluation_dataset, "\\|")[[1]]
      } else {
        # this is the csv from run_experiments
        eval_dsets = strsplit(biomarker_row$discovery_dataset, "\\|")[[1]]
      }
      
    } else {
      eval_dsets <- datasets_for_ablation_eval
    }
    
    # --- Determine Features for Ablation (if not provided) ---
    if (is.null(features_for_ablation)) {
      # message("Identifying common features across datasets specified for ablation evaluation...")
      common_features_list <- list()
      for (dset_name in eval_dsets) {
        dset_data_list <- prepared_data_list[[dset_name]]
        if (!"data_suv_bi" %in% names(dset_data_list) || !is.data.frame(dset_data_list$data_suv_bi)) next
        current_dset_features <- setdiff(names(dset_data_list$data_suv_bi), id_col)
        if (length(current_dset_features) > 0) common_features_list[[dset_name]] <- current_dset_features
      }
      if (length(common_features_list) == 0) stop("No datasets for ablation had valid SUV features.")
      features_for_ablation <- Reduce(intersect, common_features_list)
      if (length(features_for_ablation) == 0) stop("No common features found across datasets for ablation.")
      # message(sprintf("Using %d common features for ablation.", length(features_for_ablation)))
    }
    
    # Validate initial regions against features_for_ablation
    missing_initial_num <- setdiff(initial_num_regions, features_for_ablation)
    missing_initial_den <- setdiff(initial_den_regions, features_for_ablation)
    if (length(missing_initial_num) > 0 || length(missing_initial_den) > 0) {
      # warning(sprintf("Skipping biomarker ID '%s': Initial regions not all found in 'features_for_ablation'. Missing Num: [%s], Missing Den: [%s]",
      #                 original_biomarker_id_str, paste(missing_initial_num, collapse=","), paste(missing_initial_den, collapse=",")), call. = FALSE)
      # next
      warning(sprintf("Biomarker ID '%s': Initial regions not all found in 'features_for_ablation'. Using intersection. \nMissing Num: [%s], Missing Den: [%s]",
                      original_biomarker_id_str, paste(missing_initial_num, collapse=","), paste(missing_initial_den, collapse=",")), call. = FALSE)
      
      initial_num_regions <- intersect(initial_num_regions, features_for_ablation)
      initial_den_regions <- intersect(initial_den_regions, features_for_ablation)
    }
    
    
    # --- Call region_ablation ---
    if(verbose) message("  Running iterative ablation...")
    ablation_output <- try(region_ablation(
      initial_numerator_regs = initial_num_regions,
      initial_denominator_regs = initial_den_regions,
      var_composition = var_comp,
      prepared_data_list = prepared_data_list, # Pass the full list
      config = config,
      datasets_to_evaluate = eval_dsets, # Datasets for SSE in ablation
      groups_to_evaluate = eval_groups,     # Groups for SSE in ablation
      features = features_for_ablation, # The common set for this ablation run
      verbose = FALSE # Keep internal ablation less verbose, summarize here
    ), silent = TRUE)
    
    if (inherits(ablation_output, "try-error") || is.null(ablation_output)) {
      warning(sprintf("Ablation failed for biomarker ID '%s'. Error: %s",
                      original_biomarker_id_str,
                      if(inherits(ablation_output, "try-error")) conditionMessage(attr(ablation_output,"condition")) else "NULL result"),
              call. = FALSE)
      # Store a row indicating ablation failure
      summary_row <- data.frame(
        biomarker_source_id = original_biomarker_id_str,
        dplyr::select(biomarker_row, dplyr::any_of(carry_over_cols)), # Carry over original info
        initial_num_regions_str = regs_num_str,
        initial_den_regions_str = regs_den_str,
        initial_num_count = length(initial_num_regions),
        initial_den_count = length(initial_den_regions),
        ablation_status = "Failed",
        final_num_regions_str = NA_character_,
        final_den_regions_str = NA_character_,
        final_num_count = NA_integer_,
        final_den_count = NA_integer_,
        regions_removed_count = NA_integer_,
        initial_aggregated_sse = NA_real_,
        final_aggregated_sse = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      if(verbose) message("  Ablation completed.")
      # Construct summary row from ablation_output
      summary_row <- data.frame(
        biomarker_source_id = original_biomarker_id_str,
        dplyr::select(biomarker_row, dplyr::any_of(carry_over_cols)), # Carry over original info
        initial_num_regions_str = regs_num_str,
        initial_den_regions_str = regs_den_str,
        initial_num_count = length(initial_num_regions),
        initial_den_count = length(initial_den_regions),
        ablation_status = "Success",
        final_num_regions_str = paste(ablation_output$final_numerator_regs, collapse=", "),
        final_den_regions_str = paste(ablation_output$final_denominator_regs, collapse=", "),
        final_num_count = length(ablation_output$final_numerator_regs),
        final_den_count = length(ablation_output$final_denominator_regs),
        regions_removed_count = nrow(ablation_output$ablation_log),
        initial_aggregated_sse = ablation_output$initial_aggregated_sse,
        final_aggregated_sse = ablation_output$final_aggregated_sse,
        stringsAsFactors = FALSE
      )
      # Optionally, could also save ablation_output$ablation_log (detailed steps)
      # or ablation_output$per_dataset_sse_final to separate files if needed.
    }
    
    result_counter <- result_counter + 1
    all_ablation_summary_results[[result_counter]] <- summary_row
    
  } # End loop over discovered biomarkers
  
  # --- Combine All Ablation Summary Results ---
  if (length(all_ablation_summary_results) == 0) {
    warning("No ablation summary results were generated.", call. = FALSE)
    final_ablation_df <- data.frame()
  } else {
    final_ablation_df <- dplyr::bind_rows(all_ablation_summary_results)
    message(sprintf("\n--- Ablation Processing Complete: Generated %d summary rows. ---", nrow(final_ablation_df)))
  }
  
  # --- Save Final Ablation Table ---
  if (!is.null(output_ablation_results_csv_path)) {
    tryCatch({
      output_dir <- dirname(output_ablation_results_csv_path)
      if (!dir.exists(output_dir)) {
        message("Creating output directory for ablation results CSV: ", output_dir)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      if(dir.exists(output_dir)){
        readr::write_csv(final_ablation_df, output_ablation_results_csv_path, na = "NA")
        message("Ablation summary results saved to: ", output_ablation_results_csv_path)
      } else {
        warning("Failed to create output directory '", output_dir, "'. Cannot write ablation CSV.", call.=FALSE)
      }
    }, error = function(e) {
      warning(sprintf("Failed to write ablation results to '%s': %s", output_ablation_results_csv_path, e$message), call. = FALSE)
    })
  }
  
  invisible(final_ablation_df)
}