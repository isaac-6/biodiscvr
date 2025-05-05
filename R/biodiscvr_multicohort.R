#' Run Multi-Cohort Biomarker Discovery (using a Genetic Algorithm by default)
#'
#' Performs GA optimization for biomarker definitions, evaluating fitness based
#' on aggregated performance across multiple specified datasets. The GA aims
#' to maximize the aggregated fitness.
#'
#' @param preprocessed_data List. The output from `preprocess_datasets`, containing
#'   the potentially filtered `$data` list for multiple cohorts and `$config`.
#' @param datasets_to_run Character vector. Names of the datasets within
#'   `preprocessed_data` to include in the multi-cohort fitness evaluation.
#' @param group Character string. The group ("CU" or "CI") to use for evaluating
#'   fitness via the internal `.calculate_fitness` function within each cohort.
#' @param config List. The loaded configuration object from `check_and_prepare_data`.
#' @param features Character vector. Names of columns in `data_suv_bi` frames
#'   to be considered as features. If NULL (default), uses common numeric columns
#'   (excluding ID) found across *all* specified `datasets_to_run`.
#' @param fixed_numerator_regs fixed numerator regions.
#' @param fixed_denominator_regs fixed denominator regions.
#' @param min_bounds Numeric vector or single number. Lower bounds for GA parameters.
#'   Length must match the number of features if specified as a vector. Defaults to 0.1.
#' @param max_bounds Numeric vector or single number. Upper bounds for GA parameters.
#'   Length must match the number of features if specified as a vector. Defaults to 2.9.
#' @param var_composition Numeric. Parameter passed to `.calculate_fitness` indicating
#'   how features are combined (e.g., 0 = volume-weighed, 1 = arithmetic mean).
#' @param reference_fitness Numeric vector or NULL. A vector of weights, one for
#'   each dataset in `datasets_to_run`, used to calculate the weighted sum for
#'   the aggregated fitness. If NULL (default), a warning is issued and weights
#'   default to 1 for all datasets (simple sum/average). Must be the same length
#'   as `datasets_to_run`.
#' @param bilateral Logical. Metadata flag. Defaults to TRUE.
#' @param experiment_tag Character string or NULL. Metadata tag.
#' @param output_csv_name Character string or NULL. Path to CSV file for appending
#'   the main results (best solution, aggregated fitness). Directory created if needed.
#' @param save_plot Boolean, to choose if a plot of the search is saved (fitness vs generation)
#' @param output_dir Path to the output folder (to save csv, and plots)
#' @param ga_seed Numeric or NULL. Seed for GA.
#' @param ... Additional arguments passed to `GA::ga`.
#'
#' @return A list containing:
#' \describe{
#'   \item{ga_result_object}{The raw object returned by `GA::ga`.}
#'   \item{best_solution}{A list containing the decoded best biomarker definition
#'     (e.g., `numerator_regions`, `denominator_regions`, `chromosome_vector`).}
#'   \item{aggregated_fitness}{The final aggregated fitness value (higher is better)
#'     for the best solution.}
#'   \item{result_row}{The data frame row written to `output_csv_name`.}
#'   \item{per_dataset_metrics}{A data frame containing the detailed Rep, SepAB, SSE metrics
#'     for the best solution evaluated on each individual dataset included in the run.}
#' }
#' Returns `NULL` if GA fails or setup issues occur.
#'
#' @details
#' Optimizes a biomarker definition across multiple cohorts simultaneously.
#' 1. **Setup:** Identifies common features across specified datasets.
#' 2. **Multi-Fitness Wrapper:** Defines a fitness function for `GA::ga`.
#'    - This function iterates through `datasets_to_run`.
#'    - For each dataset, it calls the internal `.calculate_fitness` function (higher score is better).
#'    - It calculates a weighted sum of the fitness scores using `reference_fitness` weights.
#'    - It returns this aggregated score directly to `GA::ga` (which maximizes).
#' 3. **GA Run:** Executes `GA::ga` using the multi-fitness wrapper. Assumes `real-valued` GA type.
#' 4. **Result Processing:** Decodes the best overall chromosome found using the same logic as `biodiscvr_single`.
#' 5. **Re-evaluation (Per Dataset):** Evaluates the *best overall* chromosome on
#'    *each individual dataset* using `.feval_group` to get detailed Rep, SepAB,
#'    SSE metrics specific to that dataset for the optimal biomarker.
#' 6. **Logging:** Appends the main result (best solution, aggregated fitness, datasets used)
#'    to `output_csv_name` if provided.
#' 7. **Return:** Provides a list with the GA object, best solution details, final aggregated
#'    fitness, the result row logged, and the detailed per-dataset metrics.
#'
#' Ensure the internal `.calculate_fitness` and `.feval_group` functions exist within the package.
#' The logic for decoding the chromosome and calculating the biomarker value (likely via
#' an internal `.calculate_cvr` function as seen in `biodiscvr_single`) *must* be consistent.
#'
#' @export
#' @importFrom GA ga
#' @importFrom lme4 lmerControl
#' @importFrom longpower lmmpower
#' @importFrom dplyr left_join select all_of filter bind_rows
#' @importFrom stats na.omit as.formula setNames sd residuals coef quantile median weighted.mean
#' @importFrom rlang `%||%` .data sym is_scalar_character 
#' @importFrom methods is
#' @importFrom utils sessionInfo installed.packages packageName
biodiscvr_multicohort <- function(preprocessed_data,
                                  datasets_to_run,
                                  group,
                                  config,
                                  features = NULL,
                                  fixed_numerator_regs = NULL,
                                  fixed_denominator_regs = NULL,
                                  min_bounds = NULL,
                                  max_bounds = NULL,
                                  var_composition = 1,
                                  reference_fitness = NULL, 
                                  bilateral = TRUE,
                                  experiment_tag = NULL,
                                  output_csv_name = NULL,
                                  save_plot = FALSE,
                                  output_dir = NULL,
                                  ga_seed = NULL,
                                  ...) {
  
  # --- Input Validation ---
  if (!is.list(preprocessed_data)) {
    stop("'preprocessed_data' must be a list containing a folder per dataset.")
  }
  if (!rlang::is_character(datasets_to_run) || length(datasets_to_run) == 0) { # Use rlang for type check
    stop("'datasets_to_run' must be a non-empty character vector.")
  }
  missing_dsets <- setdiff(datasets_to_run, names(preprocessed_data))
  if (length(missing_dsets) > 0) {
    stop("Datasets not found in preprocessed_data: ", paste(missing_dsets, collapse=", "))
  }
  # Validate reference_fitness
  if (!is.null(reference_fitness)) {
    if(!is.numeric(reference_fitness) || length(reference_fitness) != length(datasets_to_run)) {
      stop("'reference_fitness' must be a numeric vector with the same length as 'datasets_to_run'.")
    }
    if(any(is.na(reference_fitness)) || any(reference_fitness == 0)) {
      stop("'reference_fitness' values cannot be NA or zero")
    }
  } else {
    warning("Argument 'reference_fitness' not provided. Using equal weights (1) for all datasets in fitness aggregation.", call. = FALSE)
    reference_fitness <- rep(1, length(datasets_to_run)) # Default: vector of 1
    names(reference_fitness) = datasets_to_run
  }
  
  # Other validations (copy/adapt from biodiscvr_single)
  stopifnot(
    is.character(group), length(group) == 1, group %in% c("CU", "CI"),
    is.list(config),
    "model_equations" %in% names(config), "power_params" %in% names(config),
    "preprocessing" %in% names(config), "genetic_algorithm" %in% names(config),
    is.null(features) || is.character(features),
    is.null(fixed_numerator_regs) || is.character(fixed_numerator_regs),         # Validate fixed regs
    is.null(fixed_denominator_regs) || is.character(fixed_denominator_regs),
    is.null(min_bounds) || is.numeric(min_bounds),
    is.null(max_bounds) || is.numeric(max_bounds),
    !is.null(var_composition), is.numeric(var_composition), length(var_composition) == 1,
    is.logical(bilateral), length(bilateral) == 1,
    is.null(experiment_tag) || rlang::is_scalar_character(experiment_tag), # Use rlang check
    is.null(output_csv_name) || rlang::is_scalar_character(output_csv_name),
    is.null(ga_seed) || (is.numeric(ga_seed) && length(ga_seed)==1),
    is.null(output_dir) || (is.character(output_dir) && length(output_dir) == 1), # Validate plot dir
    is.logical(save_plot), length(save_plot) == 1 # Validate save_plot flag
  )
  # Maybe ensure required packages are available
  # ... (checks for GA, lme4, longpower) ...
  
  if (!is.null(fixed_numerator_regs) && !is.null(fixed_denominator_regs)) {
    stop("No need to run biodiscvr if numerator and denominator are both provided.", call. = FALSE)
  }
  
  # --- Configuration & Setup ---
  id_col <- config$preprocessing$id_column %||% "RID"
  all_power_params_internal <- config$power_params
  ga_config <- config$genetic_algorithm %||% list()
  # Convert formulas once
  eq_all_internal <- try(stats::as.formula(config$model_equations$eq_all_string), silent = TRUE)
  eq_group_internal <- try(stats::as.formula(config$model_equations$eq_group_string), silent = TRUE)
  if (inherits(eq_all_internal, "try-error")) stop("Failed to parse 'eq_all_string'")
  if (inherits(eq_group_internal, "try-error")) stop("Failed to parse 'eq_group_string'")
  # lmer Control
  lmer_control_internal <- lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE, check.conv.singular = "ignore",
                                             optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
  
  
  
  
  
  
  # --- *** NEW: Pre-GA Data Check Across Cohorts *** ---
  message("--- Performing Pre-GA Data Sufficiency Check Across Cohorts ---")
  min_rows_threshold <- 30 # Use consistent thresholds (or get from config?)
  min_group_members_threshold <- 20
  target_group_dx_val <- if(group == "CI") 1L else 0L
  valid_datasets_for_mc <- character(0) # Keep track of datasets passing the check
  
  for (dset_name in datasets_to_run) {
    current_data_list <- preprocessed_data[[dset_name]]
    
    # Check if essential 'data' table exists and has required columns
    if (!"data" %in% names(current_data_list) || !is.data.frame(current_data_list$data) ||
        !all(c(id_col, "DX", "AB") %in% names(current_data_list$data))) {
      warning(sprintf("Dataset '%s': Skipping pre-check due to missing 'data' table or required columns (ID, DX, AB).", dset_name), call. = FALSE)
      next # Cannot check this dataset
    }
    
    data_clinical <- current_data_list$data
    
    # Filter for the target group and AB status needed for evaluation
    initial_group_data_for_eval <- data_clinical |>
      dplyr::filter(.data$DX == target_group_dx_val, .data$AB == TRUE, !is.na(.data$AB))
    
    n_initial_group_rows <- nrow(initial_group_data_for_eval)
    n_initial_group_ids <- dplyr::n_distinct(initial_group_data_for_eval[[id_col]])
    
    # Check if sufficient data exists for this dataset/group combination
    if (n_initial_group_ids >= min_group_members_threshold && n_initial_group_rows >= min_rows_threshold) {
      message(sprintf("   - Dataset '%s': Sufficient data found for Group '%s' (%d IDs, %d rows).",
                      dset_name, group, n_initial_group_ids, n_initial_group_rows))
      valid_datasets_for_mc <- c(valid_datasets_for_mc, dset_name) # Add to valid list
    } else {
      warning(sprintf("Dataset '%s', Group '%s': Insufficient initial data for multi-cohort fitness evaluation (Found %d IDs / %d rows, need >= %d IDs / %d rows with DX=%d, AB=TRUE). Excluding from multi-cohort run.",
                      dset_name, group, n_initial_group_ids, n_initial_group_rows, min_group_members_threshold, min_rows_threshold, target_group_dx_val), call. = FALSE)
    }
  } # End loop checking datasets
  
  # --- Check if enough datasets remain ---
  min_required_cohorts <- 2 # Define minimum needed for multi-cohort analysis
  if (length(valid_datasets_for_mc) < min_required_cohorts) {
    warning(sprintf("Insufficient number of datasets (%d) passed the pre-GA data check for Group '%s' (minimum required: %d). Skipping multi-cohort GA.",
                    length(valid_datasets_for_mc), group, min_required_cohorts), call. = FALSE)
    
    # --- Construct partial result row indicating skip ---
    # (Similar structure to the one in biodiscvr_single's skip logic)
    result_row <- data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      experiment_tag = experiment_tag %||% NA_character_,
      # discovery_dataset field doesn't apply here, use datasets_included
      datasets_included = paste(datasets_to_run, collapse = "|"), # Original requested
      group_evaluated = group,
      bilateral = bilateral,
      fitness_aggregation = "skipped", # Indicate skip
      fitness_value_aggregated = NA_real_,
      regs_numerator = "SKIPPED_INSUFFICIENT_COHORTS",
      regs_denominator = "SKIPPED_INSUFFICIENT_COHORTS",
      var_composition = var_composition,
      reference_vector = "SKIPPED_INSUFFICIENT_COHORTS",
      # Add NA placeholders for metrics if your CSV expects them
      Rep = NA_real_, SepAB = NA_real_, SSE = NA_real_,
      # Add NA GA params
      ga_type = "skipped", ga_popSize = NA, ga_maxiter = NA, ga_seed_used = ga_seed %||% NA, ga_runtime_sec = 0,
      stringsAsFactors = FALSE
    )
    
    # --- Append skip info to CSV (Optional) ---
    if (!is.null(output_csv_name)) {
      # ... (logic to ensure output_dir exists and call .append_to_csv) ...
      output_dir <- dirname(output_csv_name) # Assuming output_csv_name is the full path
      if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) }
      if(dir.exists(output_dir)){ .append_to_csv(result_row, output_csv_name) }
      else { warning("Failed to create output directory '", output_dir, "'. Cannot write skipped MC CSV.", call.=FALSE) }
    }
    
    # --- Return indicating skip ---
    # Return NULL as the main function failed to proceed
    return(NULL)
  } else {
    # --- Update datasets_to_run for the rest of the function ---
    # Only proceed with the datasets that passed the check
    message(sprintf("Proceeding with multi-cohort analysis for Group '%s' using %d valid datasets: %s",
                    group, length(valid_datasets_for_mc), paste(valid_datasets_for_mc, collapse=", ")))
    datasets_to_run <- valid_datasets_for_mc # IMPORTANT: Overwrite the list for subsequent steps
    
    # Also need to filter the reference_fitness vector if it was provided
    if (!is.null(reference_fitness)) {
      original_names <- names(reference_fitness) # Assumes it was named by dataset
      if (is.null(original_names)) stop("Provided 'reference_fitness' must be a named vector.")
      keep_indices <- original_names %in% datasets_to_run
      reference_fitness <- reference_fitness[keep_indices]
      # Re-normalize weights if needed (optional, depends on aggregation method)
      # if(sum(reference_fitness) > 0) reference_fitness <- reference_fitness / sum(reference_fitness)
      # else reference_fitness <- rep(1/length(reference_fitness), length(reference_fitness))
      message(sprintf("   - Filtered reference_fitness vector to match %d valid datasets.", length(reference_fitness)))
    }
  }
  # --- *** END: Pre-GA Data Check Across Cohorts *** ---
  
  

  
  # --- Optimized Version ---
  
  message("Identifying common features across specified datasets...")
  id_col <- config$preprocessing$id_column %||% "RID"
  
  # Initialize common features tracking
  common_features <- NULL
  processed_at_least_one <- FALSE
  
  for (dset_name in datasets_to_run) {
    dset_data_list <- preprocessed_data[[dset_name]]
    
    # Ensure dataset contains expected data
    if (!"data_suv_bi" %in% names(dset_data_list) || is.null(dset_data_list$data_suv_bi)) {
      warning(sprintf("Dataset '%s' is missing 'data_suv_bi'. Skipping.", dset_name), call. = FALSE)
      next
    }
    
    # Extract valid features
    current_features_valid <- setdiff(names(dset_data_list$data_suv_bi), id_col)
    
    if (length(current_features_valid) == 0) {
      warning(sprintf("Dataset '%s' has no valid feature columns in 'data_suv_bi'. Skipping.", dset_name), call. = FALSE)
      next
    }
    
    # Establish common features
    if (!processed_at_least_one) {
      common_features <- current_features_valid
      processed_at_least_one <- TRUE
    } else {
      common_features <- intersect(common_features, current_features_valid)
    }
    
    # Early exit if intersection is empty
    if (length(common_features) == 0) {
      stop(sprintf("No common features remain after processing '%s'. Cannot proceed.", dset_name))
    }
  }
  
  # Ensure at least one dataset was processed
  if (!processed_at_least_one || length(common_features) == 0) {
    stop("No common features found across all datasets. Analysis cannot proceed.")
  }
  
  message(sprintf("Found %d common features across %d datasets.", length(common_features), length(datasets_to_run)))
  
  # --- Validate Fixed Regions ---
  validate_fixed_regions <- function(fixed_regs, common_features, type) {
    if (is.null(fixed_regs)) return(NULL)
    
    valid_regs <- intersect(fixed_regs, common_features)
    excluded_regs <- setdiff(fixed_regs, valid_regs)
    
    if (length(excluded_regs) > 0) {
      warning(sprintf("Excluded %s not found in all datasets: %s", type, paste(excluded_regs, collapse = ", ")), call. = FALSE)
    }
    
    if (length(valid_regs) == 0) {
      stop(sprintf("No valid %s common across all datasets. Cannot proceed.", type))
    }
    
    message(sprintf("Using %d valid %s (common across all datasets).", length(valid_regs), type))
    return(valid_regs)
  }
  
  fixed_numerator_regs <- validate_fixed_regions(fixed_numerator_regs, common_features, "fixed numerator regions")
  fixed_denominator_regs <- validate_fixed_regions(fixed_denominator_regs, common_features, "fixed denominator regions")
  
  # --- Determine Features for GA ---
  # Adjust features if a user override is provided
  if (!is.null(features)) {
    features <- setdiff(features, c(fixed_numerator_regs, fixed_denominator_regs))
    features <- intersect(features, common_features)
    
    if (length(features) == 0) stop("User-specified features do not overlap with common features.")
    
    n_features <- length(features)
    
    message(sprintf("Using %d user-specified common features.", length(features)))
  }
  
  
  # --- GA Bounds & Variables (n_vars = n_features) ---
  n_vars <- n_features # Correction based on feedback
  # Assume real-valued only
  min_bounds_default <- 0.1; max_bounds_default <- 2.9
  min_bounds <- min_bounds %||% min_bounds_default
  max_bounds <- max_bounds %||% max_bounds_default
  if (length(min_bounds) == 1) min_bounds <- rep(min_bounds, n_vars)
  if (length(max_bounds) == 1) max_bounds <- rep(max_bounds, n_vars)
  if (length(min_bounds) != n_vars || length(max_bounds) != n_vars) {
    stop(sprintf("GA bounds dimensions (%d, %d) do not match n_features (%d)", length(min_bounds), length(max_bounds), n_features))
  }
  
  # --- Constrain literature-based regions ---
  # This can avoid having two directions of well-performing biomarkers,
  # although the fact that a negative t is currently penalised makes it redundant.
  # This future-proofs other fitness metrics, like focusing on SSE alone.
  aux_t <- c("entorhinal", "amygdala", "superiortemporal")
  max_bounds[features %in% aux_t] <- 1.9
  
  aux_r <- c("brainstem", "hemiwm", "whole_cerebellum", "cerebellum_cortex", "inferior_cerebgm")
  min_bounds[features %in% aux_r] <- 1.1
  
  rm(aux_t, aux_r)
  
  
  # --- Define Multi-Cohort Fitness Wrapper ---
  # Capture variables needed from this environment
  datasets_to_run_internal <- datasets_to_run
  preprocessed_data_internal <- preprocessed_data
  var_composition_internal <- var_composition
  group_internal <- group
  # reference_fitness_internal <- reference_fitness
  
  .multi_fitness_wrapper <- function(chromosome) {
    individual_fitness_values <- rep(NA_real_, length(datasets_to_run_internal))
    names(individual_fitness_values) <- datasets_to_run_internal
    
    for (i in seq_along(datasets_to_run_internal)) {
      dset_name_internal <- datasets_to_run_internal[i]
      dataset_data_internal <- preprocessed_data_internal[[dset_name_internal]]
      
      fitness_value <- tryCatch({
        .calculate_fitness( # Call the package's internal function
          chromosome = chromosome,
          features = features,
          fixed_numerator_regs = fixed_numerator_regs,
          fixed_denominator_regs = fixed_denominator_regs,
          var_composition = var_composition_internal,
          cohort_data = dataset_data_internal,
          group = group_internal,
          eq_all = eq_all_internal,
          eq_group = eq_group_internal,
          all_power_params = all_power_params_internal,
          lmer_control = lmer_control_internal
        )
      }, error = function(e) { NA_real_ }) # Return NA on error
      
      individual_fitness_values[i] <- fitness_value
    } # End dataset loop
    
    # Aggregate Fitness
    valid_indices <- !is.na(individual_fitness_values) & is.finite(individual_fitness_values)
    if (sum(valid_indices) < length(datasets_to_run_internal)) {
      return(-Inf) # No valid fitness calculated
    } else {
      # Aggregate
      aggregated_fitness <- .multi_fitness(individual_fitness_values, reference_fitness)
    }
    
    # Return value directly (GA::ga maximizes)
    if(is.na(aggregated_fitness) || !is.finite(aggregated_fitness)){
      return(-Inf) # Return worst possible value if aggregation fails
    } else {
      return(aggregated_fitness)
    }
  } # End .multi_fitness_wrapper definition
  
  
  # --- Run the GA ---
  message(sprintf("--- Running Multi-Cohort GA (Datasets: %s; Group: %s) ---",
                  paste(datasets_to_run, collapse=", "), group))
  # --- Prepare GA args ---
  ga_popSize <- ga_config$popSize %||% 50
  ga_maxiter <- ga_config$maxiter %||% 100
  # ... (set ga_run, ga_pmutation, ga_pcrossover, ga_elitism from ga_config) ...
  ga_run <- ga_config$run %||% (ga_maxiter + 10)
  ga_pmutation <- ga_config$pmutation %||% 0.5
  ga_pcrossover <- ga_config$pcrossover %||% 0.8
  elitism_raw <- ga_config$elitism_prop %||% 2
  ga_parallel <- ga_config$parallel %||% FALSE
  ga_monitor <- ga_config$monitor %||% FALSE
  ga_elitism <- max(1, round(elitism_raw))
  if(elitism_raw > 0 && elitism_raw < 1) ga_elitism <- max(1, round(elitism_raw * ga_popSize))
  effective_seed <- ga_seed %||% sample.int(.Machine$integer.max, 1)
  
  ga_args <- list(
    type = "real-valued", # Only support real-valued
    fitness = .multi_fitness_wrapper, # Use the multi-cohort wrapper
    lower = min_bounds,
    upper = max_bounds,
    popSize = ga_popSize,
    maxiter = ga_maxiter,
    run = ga_run,
    pmutation = ga_pmutation,
    pcrossover = ga_pcrossover,
    elitism = ga_elitism,
    monitor = ga_monitor,
    parallel = ga_parallel,
    seed = effective_seed
  )
  # Add specific real-valued operators if in config
  if(!is.null(ga_config$mutation)) ga_args$mutation <- ga_config$mutation
  if(!is.null(ga_config$crossover)) ga_args$crossover <- ga_config$crossover
  # Add extra arguments from ...
  extra_args <- list(...)
  ga_args <- c(ga_args, extra_args[!names(extra_args) %in% names(ga_args)])
  
  # --- Execute GA ---
  ga_result_obj <- try(do.call(GA::ga, ga_args), silent = TRUE)
  
  # --- Process GA Results ---
  if (inherits(ga_result_obj, "try-error")) {
    warning(sprintf("Multi-Cohort GA failed. Error: %s", conditionMessage(attr(ga_result_obj, "condition"))), call.=FALSE)
    return(NULL)
  }
  message("--- Multi-Cohort GA completed ---")
  
  
  
  
  
  
  # --- 8b. Save GA Convergence Plot (Optional) ---
  if (save_plot && !is.null(output_dir) && !inherits(ga_result_obj, "try-error")) {
    
    # --- Create Directory if needed ---
    # (Attempt creation even if just saving one plot, might be first time)
    if (!dir.exists(output_dir)) {
      message("Creating output directory for plot: ", output_dir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) # Use showWarnings=FALSE for cleaner output if dir exists
    }
    
    # --- Proceed only if directory exists (or was successfully created) ---
    if(dir.exists(output_dir)){
      # --- Construct Filename ---
      plot_file_base <- paste(
        paste(datasets_to_run, collapse="."),
        group,
        experiment_tag %||% "untagged", # Handle NULL experiment tag
        "discovery_progression",
        sep = "_"
      )
      # Replace potentially problematic characters
      plot_file_base <- gsub("[^a-zA-Z0-9_.-]", "_", plot_file_base)
      plot_filename <- paste0(plot_file_base, ".tiff")
      full_plot_path <- file.path(output_dir, plot_filename)
      
      # --- Attempt to save plot (will overwrite if exists) ---
      message("   - Saving GA plot to: ", full_plot_path)
      tryCatch({
        # Use grDevices::tiff for saving
        grDevices::tiff(filename = full_plot_path,
                        width = 7, height = 5, units = "in", # Adjust size
                        res = 300, # Adjust resolution
                        compression = "lzw")
        
        # Create the plot (add a title)
        plot_title <- sprintf("GA Convergence: %s (%s%s)",
                              paste(datasets_to_run, collapse="."),
                              group,
                              if(!is.null(experiment_tag)) paste0(" - ", experiment_tag) else "" )
        plot(ga_result_obj, main = plot_title)
        
        # Close the device to write the file
        grDevices::dev.off()
        
      }, error = function(e) {
        warning(sprintf("Dataset '%s', Group '%s': Failed to save GA plot to '%s'. Error: %s",
                        paste(datasets_to_run, collapse="."), group, full_plot_path, conditionMessage(e)), call. = FALSE)
        # Ensure device is closed if error occurred after opening but before dev.off()
        if (grDevices::dev.cur() != 1L) try(grDevices::dev.off(), silent=TRUE) # Use try() for safety
      }) # End tryCatch
    } else {
      warning("Output plot directory '", output_dir, "' does not exist and could not be created. Cannot save plot.", call.=FALSE)
    } # End if dir.exists
  } else if (save_plot && is.null(output_dir)) {
    warning("`save_plot` is TRUE, but `output_dir` was not provided. Plot not saved.", call. = FALSE)
  } # End if save_plot
  
  
  
  
  best_chromosome <- ga_result_obj@solution[1, ]
  # Fitness value from GA object is already the maximized aggregated value
  final_aggregated_fitness <- ga_result_obj@fitnessValue[1]
  
  
  # --- Decode Best Chromosome (Logic from biodiscvr_single) ---
  best_regs_numerator <- character(0)
  best_regs_denominator <- character(0)
  best_regs_numerator <- features[best_chromosome < 1]
  best_regs_denominator <- features[best_chromosome > 2]
  
  # Override numerator if fixed_numerator_regs is provided
  if (!is.null(fixed_numerator_regs)) {
    best_regs_numerator <- fixed_numerator_regs
  }
  
  # Override denominator if fixed_denominator_regs is provided
  if (!is.null(fixed_denominator_regs)) {
    best_regs_denominator <- fixed_denominator_regs
  }
  
  if (length(best_regs_numerator) == 0 || length(best_regs_denominator) == 0) {
    warning("Best solution from multi-cohort GA resulted in empty numerator/denominator. Cannot evaluate further.", call.=FALSE)
    # Construct a partial result row indicating failure after GA
    result_row <- data.frame( # Create a row with NAs for metrics
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      experiment_tag = experiment_tag %||% NA_character_,
      datasets_included = paste(datasets_to_run, collapse = "|"),
      group_evaluated = group, bilateral = bilateral,
      fitness_aggregation = "weighted_sum", # Indicate aggregation type
      fitness_value_aggregated = final_aggregated_fitness,
      regs_numerator = "DECODE_FAILED", regs_denominator = "DECODE_FAILED",
      var_composition = var_composition,
      # ... Add GA params ...
      Rep = NA_real_, SepAB = NA_real_, SSE = NA_real_, # Add NA placeholders for metrics
      stringsAsFactors = FALSE
    )
    # Optionally log this failure row
    if (!is.null(output_csv_name)) { .append_to_csv(result_row, output_csv_name) }
    # Return list indicating partial success
    return(list(ga_result_object=ga_result_obj, best_solution=list(chromosome=best_chromosome), aggregated_fitness=final_aggregated_fitness, result_row=result_row, per_dataset_metrics=NULL))
  }
  
  
  # --- Re-evaluate Best Solution on Each Dataset ---
  message("--- Evaluating best solution on individual datasets ---")
  per_dataset_results_list <- list() # Initialize list to store results rows
  
  for (dset_name in datasets_to_run) {
    message(" - Evaluating on: ", dset_name)
    current_dset_data_list <- preprocessed_data[[dset_name]] # Get the list for this dataset
    
    # Initialize metrics and status for this dataset
    metrics <- stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
    eval_status <- "Evaluation Pending" # Initial status
    
    # --- Calculate biomarker value using BEST chromosome ---
    # Use the same function/logic as biodiscvr_single's re-evaluation
    # Assumes existence of internal .calculate_cvr or similar
    final_value <- tryCatch({
      # Example call, replace with your actual calculation function/logic
      # This MUST be consistent with how value is calculated in fitness functions
      .calculate_cvr(
        chromosome = best_chromosome, # The best one found by multi-cohort GA
        dataset_cohort_data = current_dset_data_list, # Pass this dataset's list
        features = features, # Common features used in GA
        var_composition = var_composition,
        # Pass fixed regs used in THIS multi-cohort run (might be NULL)
        fixed_numerator_regs = fixed_numerator_regs,
        fixed_denominator_regs = fixed_denominator_regs,
        verbose = FALSE # Typically keep internal calls non-verbose
      )
    }, error = function(e) {
      warning(sprintf("Dataset '%s': Error during final value calculation for best solution: %s",
                      dset_name, conditionMessage(e)), call.=FALSE)
      return(NULL) # Return NULL on error
    })
    
    if (is.null(final_value)) {
      warning(sprintf("Dataset '%s': Failed final value calculation for best solution. Metrics cannot be calculated.", dset_name), call.=FALSE)
      eval_status <- "Biomarker Calculation Failed"
      # metrics remain NA
    } else {
      # --- Prepare final data for evaluation ---
      # Add the calculated 'value' column to the 'data' part for .feval_group
      # Ensure essential columns exist in the clinical data part
      required_clinical_cols <- c(id_col, "time", "DX", "AB") # Align with .feval_group needs
      if (!"data" %in% names(current_dset_data_list) ||
          !all(required_clinical_cols %in% names(current_dset_data_list$data))) {
        warning(sprintf("Dataset '%s': Missing 'data' table or required clinical columns. Cannot calculate metrics.", dset_name), call.=FALSE)
        eval_status <- "Missing Clinical Data"
        # metrics remain NA
      } else {
        final_data_for_eval <- current_dset_data_list$data
        # Add/overwrite value column safely
        if(nrow(final_data_for_eval) == length(final_value)) {
          final_data_for_eval$value <- final_value
        } else {
          warning(sprintf("Dataset '%s': Length mismatch between final value (%d) and data rows (%d). Cannot merge for metrics.",
                          dset_name, length(final_value), nrow(final_data_for_eval)), call.=FALSE)
          eval_status <- "Value/Data Length Mismatch"
          # metrics remain NA - skip to storing results
          final_data_for_eval <- NULL # Prevent further processing
        }
        
        if (!is.null(final_data_for_eval)) {
          # Remove rows with NA in essential columns AFTER adding value
          essential_eval_cols <- c("value", "time", id_col, "DX", "AB")
          final_data_for_eval <- stats::na.omit(final_data_for_eval[, essential_eval_cols, drop = FALSE])
          
          # --- Get final metrics using .feval_group ---
          min_rows_threshold <- 30 # Example threshold
          min_group_members_threshold <- 20 # Example threshold
          current_group_dx_val <- if(group == "CI") 1L else 0L
          current_group_member_count <- sum(final_data_for_eval$DX == current_group_dx_val, na.rm = TRUE)
          
          if(nrow(final_data_for_eval) < min_rows_threshold || current_group_member_count < min_group_members_threshold ) {
            warning(sprintf("Dataset '%s', Group '%s': Not enough final data rows (%d) or group members (%d) to calculate metrics reliably.",
                            dset_name, group, nrow(final_data_for_eval), current_group_member_count), call.=FALSE)
            eval_status <- "Insufficient Data for Metrics"
            # metrics remain NA
          } else {
            # Attempt to calculate metrics
            metrics_result <- try(.feval_group(
              data = final_data_for_eval,
              group = group,
              eq_all = eq_all_internal, # Use internal formula objects
              eq_group = eq_group_internal,
              all_power_params = all_power_params_internal,
              lmer_control = lmer_control_internal
            ), silent = TRUE)
            
            if (inherits(metrics_result, "try-error")) {
              warning(sprintf("Dataset '%s', Group '%s': .feval_group failed during re-evaluation. Error: %s",
                              dset_name, group, conditionMessage(attr(metrics_result, "condition"))), call.=FALSE)
              eval_status <- "Metrics Calculation Error"
              # metrics remain NA
            } else if (anyNA(metrics_result)) {
              warning(sprintf("Dataset '%s', Group '%s': .feval_group returned NA values during re-evaluation.",
                              dset_name, group), call.=FALSE)
              eval_status <- "Metrics Calculation NA"
              metrics <- metrics_result # Store the NAs
            } else {
              # Success!
              metrics <- metrics_result
              eval_status <- "Success"
            }
          } # End check for sufficient data
        } # End if final_data_for_eval is valid
      } # End check for clinical data validity
    } # End check for final_value validity
    
    # --- Store per-dataset metrics including status ---
    per_dataset_results_list[[dset_name]] <- data.frame(
      dataset = dset_name,
      group_evaluated = group,
      status = eval_status, # Add the status column
      Rep = metrics["Rep"],
      SepAB = metrics["SepAB"],
      SSE = metrics["SSE"],
      stringsAsFactors = FALSE
    )
    
  } # End loop re-evaluating datasets
  
  # Convert the list to a single data frame
  per_dataset_metrics_df <- dplyr::bind_rows(per_dataset_results_list)
  
  # --- Helper function for safe rounding and extraction ---
  extract_round_metric <- function(metric_list, metric_name, order, digits) {
    vals <- sapply(order, function(name) {
      res <- metric_list[[name]][metric_name]
      if (is.null(res) || is.na(res)) return(NA_character_) # Handle missing results
      format(round(res, digits), nsmall = digits)
    })
    paste(vals, collapse = "_")
  }
  
  # --- Concatenate Metrics ---
  evaluation_order <- datasets_to_run
  Rep_All <- extract_round_metric(per_dataset_results_list, "Rep", evaluation_order, 3)
  SepAB_All <- extract_round_metric(per_dataset_results_list, "SepAB", evaluation_order, 2)
  SSE_All <- extract_round_metric(per_dataset_results_list, "SSE", evaluation_order, 1)
  
  
  # --- Construct Main Result Row ---
  result_row <- data.frame(
    # === Metadata ===
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    experiment_tag = experiment_tag %||% NA_character_,
    discovery_dataset = paste(datasets_to_run, collapse = "|"),
    group_evaluated = group,          # Group used for fitness/metrics
    bilateral = bilateral,
    var_composition = var_composition, # Logged input parameter
    
    # === Fitness Info ===
    # fitness_metric_used = "custom", # Less informative now we know it's .calculate_fitness
    # maximized_fitness = TRUE,      # Implied by design
    fitness_value = final_aggregated_fitness,
    
    # === Evaluation Metrics (Prefixed) ===
    # This adds columns like Rep_CI, SepAB_CI, SSE_CI etc. dynamically
    # metrics_for_df,
    
    # as.list(final_metrics),
    # === Concatenated Evaluation Metrics ===
    Rep_All = Rep_All,
    SepAB_All = SepAB_All,
    SSE_All = SSE_All,
    
    # === Biomarker Definition ===
    regs_numerator = paste(best_regs_numerator, collapse = ", "), # Use '|' as separator
    regs_denominator = paste(best_regs_denominator, collapse = ", "),
    reference_vector = 1, # Added reference vector. This function only has a single cohort.
    
    # === GA Info ===
    # ga_type = ga_result_obj@type,
    ga_popSize = ga_result_obj@popSize,
    ga_maxiter = ga_result_obj@maxiter,
    ga_run = ga_result_obj@run,
    ga_iter = ga_result_obj@iter,
    ga_seed_used = effective_seed,
    
    # Ensure strings aren't factors
    stringsAsFactors = FALSE
  )
  
  # --- Append Main Results to CSV (Optional) ---
  output_csv_path <- paste0(
    ifelse(substring(output_dir, nchar(output_dir)) == "/", output_dir, paste0(output_dir, "/")),
    ifelse(grepl("\\.csv$", output_csv_name, ignore.case = TRUE), output_csv_name, paste0(output_csv_name, ".csv"))
  )
  if (!is.null(output_csv_path)) {
    output_dir <- dirname(output_csv_path)
    if (!dir.exists(output_dir)) {
      message("(Re)Creating output directory for CSV: ", output_dir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if(dir.exists(output_dir)){
      append_success <- .append_to_csv(result_row, output_csv_path)
      if(append_success) message("   - Results appended to: ", output_csv_path)
    } else {
      warning("Failed to create output directory '", output_dir, "'. Cannot write CSV.", call.=FALSE)
    }
  }
  
  
  # --- Return Structured Results ---
  return(list(
    # result_row contains the main aggregated fitness, regions (string), GA params etc.
    # This is the primary record for the CSV log.
    result_row = result_row,
    best_regs_numerator = best_regs_numerator,
    best_regs_denominator = best_regs_denominator
    
    # per_dataset_metrics provides the breakdown of performance on individual cohorts
    # for the best multi-cohort solution. Essential for interpretation.
    # per_dataset_metrics = per_dataset_metrics_df # Returned as data frame
    
    # Optionally include the decoded best solution if needed upstream, otherwise omit
    # best_solution = list(
    #   numerator_regions = best_regs_numerator, # Vector
    #   denominator_regions = best_regs_denominator, # Vector
    #   chromosome_vector = best_chromosome
    # ),
    
    # Optionally include raw GA object if deep diagnostics are needed, otherwise omit
    # ga_result_object = ga_result_obj
    
    # aggregated_fitness is already captured within result_row
  ))
}
