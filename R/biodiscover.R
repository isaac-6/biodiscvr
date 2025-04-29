#' Run Single-Dataset Biomarker Discovery using Genetic Algorithm
#'
#' Performs GA optimization to find biomarker definitions (e.g., feature weights
#' or region selections) for a single dataset, maximizing a custom fitness
#' function based on internal evaluation metrics.
#'
#' @param dataset_data List containing the prepared 'data' and 'data_suv_bi'
#'   data frames for the target dataset (e.g., `prepared_data_list$ADNI`).
#' @param dataset_name Character string. The name of the dataset being processed
#'   (e.g., "ADNI"). Used for metadata logging.
#' @param group Character string. The group ("CU" or "CI") whose specific data
#'   characteristics will be used for the primary fitness evaluation by the
#'   internal `.calculate_fitness` function.
#' @param config List. The loaded configuration object, typically from
#'   `check_and_prepare_data$config`. Must contain sections like
#'   `model_equations`, `power_params`, `preprocessing`, and `genetic_algorithm`.
#' @param features Character vector. Names of columns in `dataset_data$data_suv_bi`
#'   to be considered as features by the GA. If NULL (default), uses all numeric
#'   columns excluding the ID column specified in `config$preprocessing$id_column`.
#' @param fixed_numerator_regs fixed numerator regions.
#' @param fixed_denominator_regs fixed denominator regions.
#' @param min_bounds Numeric vector or single number. Lower bounds for GA parameters.
#'   Matched to the number of variables in the GA encoding (e.g., number of features
#'   if using real-valued weights). Defaults to 0.1 if NULL.
#' @param max_bounds Numeric vector or single number. Upper bounds for GA parameters.
#'   Matched to the number of variables in the GA encoding. Defaults to 2.9 if NULL.
#' @param var_composition Numeric. A parameter passed to the internal
#'   `.calculate_fitness` function, likely indicating how regions/features
#'   are combined (e.g., 0 = volume-weighed, 1 = arithmetic mean, etc.).
#' @param bilateral Logical. Metadata flag indicating if the analysis considers
#'   bilateral regions (logged in the output). Defaults to TRUE.
#' @param experiment_tag Character string or NULL. An optional tag to identify this
#'   specific run/experiment batch in the output CSV.
#' @param output_csv_name Character string or NULL. If provided, the full path
#'   to a CSV file where the single result row for this run will be appended.
#'   Handles header creation automatically. Directory will be created if needed.
#' @param save_plot Boolean, to choose if a plot of the search is saved (fitness vs generation)
#' @param output_dir Path to the output folder (to save csv, and plots)
#' @param ga_seed Numeric or NULL. Seed for the GA's random number generator
#'   for reproducibility. If NULL, a random seed is used.
#' @param ... Additional arguments passed directly to the underlying GA function
#'   (e.g., `GA::ga`).
#'
#' @return A single-row data frame containing the results for the best solution
#'   found by the GA, matching the structure intended for the output CSV.
#'   Columns include metadata, fitness value, final evaluation metrics (Rep, SepAB,
#'   SSE for the evaluated `group`), the best biomarker definition (e.g.,
#'   regions), and GA parameters. Returns `NULL` if the GA fails, essential
#'   data is missing, or critical steps cannot be completed.
#'
#' @details
#' This function orchestrates the GA optimization for a single dataset:
#' 1.  **Initialization:** Sets up parameters, identifies features, prepares bounds.
#' 2.  **Fitness Function:** Defines an internal wrapper around the package's
#'     `.calculate_fitness` function. This wrapper handles:
#'       a. Decoding the GA chromosome (solution vector) into a biomarker definition.
#'          **(PLACEHOLDER - requires specific implementation)**
#'       b. Calculating the biomarker `value` based on the definition and the
#'          dataset's SUV data. **(PLACEHOLDER - requires specific implementation)**
#'       c. Preparing the necessary data subset for evaluation.
#'       d. Calling `.calculate_fitness` with all required arguments.
#'       e. Returning the negative of the fitness value (since `.calculate_fitness`
#'          returns "higher is better" and `GA::ga` minimizes).
#' 3.  **GA Execution:** Runs the genetic algorithm (using `GA::ga` by default)
#'     with the wrapper fitness function and parameters derived from the `config`
#'     and function arguments.
#' 4.  **Result Processing:** Extracts the best solution found by the GA.
#' 5.  **Re-evaluation:** Recalculates the biomarker value using the best solution
#'     and calls the internal `.feval_group` function to obtain the final set of
#'     evaluation metrics (Rep, SepAB, SSE) for the specified `group`.
#' 6.  **Output:** Constructs a single-row data frame containing comprehensive
#'     results and metadata.
#' 7.  **CSV Logging (Optional):** If `output_csv_path` is provided, appends the
#'     result row to the specified CSV file using the internal `.append_to_csv` helper.
#'
#' Ensure the internal `.calculate_fitness` function exists within the package
#' and accepts the arguments passed by the wrapper function defined herein.
#' The PLACEHOLDER sections for decoding the GA solution and calculating the
#' biomarker value *must* be implemented according to your specific methodology.
#'
#' @export
#' @importFrom GA ga
#' @importFrom lme4 lmerControl
#' @importFrom longpower lmmpower
#' @importFrom dplyr left_join select all_of filter 
#' @importFrom stats na.omit as.formula setNames sd residuals coef quantile 
#' @importFrom rlang `%||%` .data sym 
#' @importFrom methods is 
biodiscvr_single <- function(dataset_data,
                             dataset_name,
                             group,
                             config,
                             features = NULL,
                             fixed_numerator_regs = NULL,
                             fixed_denominator_regs = NULL,
                             min_bounds = NULL,
                             max_bounds = NULL,
                             var_composition, 
                             bilateral = TRUE,
                             experiment_tag = NULL,
                             output_csv_name = NULL, # Optional CSV logging
                             save_plot = FALSE,
                             output_dir = NULL,
                             ga_seed = NULL,
                             ...) {
  
  # --- Input Validation & Defaults ---
  stopifnot(
    is.list(dataset_data), all(c("data", "data_suv_bi") %in% names(dataset_data)),
    is.data.frame(dataset_data$data), is.data.frame(dataset_data$data_suv_bi),
    is.character(dataset_name), length(dataset_name) == 1,
    is.character(group), length(group) == 1, group %in% c("CU", "CI"), # Extend if needed
    is.list(config),
    # Check essential config sections needed later
    "model_equations" %in% names(config), "power_params" %in% names(config),
    "preprocessing" %in% names(config), "genetic_algorithm" %in% names(config),
    is.null(features) || is.character(features),
    is.null(min_bounds) || is.numeric(min_bounds),
    is.null(max_bounds) || is.numeric(max_bounds),
    !is.null(var_composition), is.numeric(var_composition), length(var_composition) == 1, # Validate var_composition
    is.logical(bilateral), length(bilateral) == 1,
    is.null(experiment_tag) || (is.character(experiment_tag) && length(experiment_tag)==1),
    is.null(output_csv_name) || (is.character(output_csv_name) && length(output_csv_name)==1),
    is.null(ga_seed) || (is.numeric(ga_seed) && length(ga_seed)==1),
    is.null(output_dir) || (is.character(output_dir) && length(output_dir) == 1), # Validate plot dir
    is.null(fixed_numerator_regs) || is.character(fixed_numerator_regs),         # Validate fixed regs
    is.null(fixed_denominator_regs) || is.character(fixed_denominator_regs),
    is.logical(save_plot), length(save_plot) == 1 # Validate save_plot flag
  )
  
  if (!is.null(fixed_numerator_regs) && !is.null(fixed_denominator_regs)) {
    stop("No need to run biodiscvr if numerator and denominator are both provided.", call. = FALSE)
  }
  
  # --- Ensure required packages are available ---
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop("Package 'GA' needed for biodiscvr_single function. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' needed for biodiscvr_single function. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("longpower", quietly = TRUE)) {
    stop("Package 'longpower' needed for biodiscvr_single function. Please install it.", call. = FALSE)
  }
  
  
  # --- Configuration Extraction ---
  id_col <- config$preprocessing$id_column %||% "RID"
  all_power_params <- config$power_params
  ga_config <- config$genetic_algorithm %||% list() # Safely get GA config
  
  # --- Formula Conversion ---
  eq_all_str <- config$model_equations$eq_all_string
  eq_group_str <- config$model_equations$eq_group_string
  if(is.null(eq_all_str) || is.null(eq_group_str)) stop("Missing 'eq_all_string' or 'eq_group_string' in config$model_equations")
  eq_all <- try(stats::as.formula(eq_all_str), silent = TRUE)
  eq_group <- try(stats::as.formula(eq_group_str), silent = TRUE)
  if (inherits(eq_all, "try-error")) stop("Failed to parse 'eq_all_string': ", eq_all_str)
  if (inherits(eq_group, "try-error")) stop("Failed to parse 'eq_group_string': ", eq_group_str)
  
  # --- lmer Control ---
  lmer_control <- lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                                    check.conv.singular = "ignore",
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
  
  
  # --- Validate and Filter Fixed Regions (if provided) ---
  # Work with copies so we don't modify the original arguments directly in this scope yet
  # We will reassign to the original names after validation.
  valid_fixed_num_regs <- fixed_numerator_regs
  valid_fixed_den_regs <- fixed_denominator_regs
  all_features <- names(dataset_data$data_suv_bi)
  
  if (!is.null(fixed_numerator_regs)) {
    initial_fixed_num <- fixed_numerator_regs
    valid_fixed_num_regs <- intersect(initial_fixed_num, all_features)
    missing_fixed_num <- setdiff(initial_fixed_num, valid_fixed_num_regs)
    
    if (length(missing_fixed_num) > 0) {
      warning(sprintf("Dataset '%s': Provided fixed_numerator_regs not found in available all_features and were excluded: %s",
                      dataset_name, paste(missing_fixed_num, collapse = ", ")), call. = FALSE)
    }
    if (length(valid_fixed_num_regs) == 0) {
      warning(sprintf("Dataset '%s': None of the provided fixed_numerator_regs were found in available all_features Cannot proceed.",
                      dataset_name), call. = FALSE)
      return(NULL) # Fail if NO valid fixed numerator regions remain
    }
    message(sprintf("   - Using %d valid fixed numerator regions.", length(valid_fixed_num_regs)))
  }
  
  if (!is.null(fixed_denominator_regs)) {
    initial_fixed_den <- fixed_denominator_regs
    valid_fixed_den_regs <- intersect(initial_fixed_den, all_features)
    missing_fixed_den <- setdiff(initial_fixed_den, valid_fixed_den_regs)
    
    if (length(missing_fixed_den) > 0) {
      warning(sprintf("Dataset '%s': Provided fixed_denominator_regs not found in available all_features and were excluded: %s",
                      dataset_name, paste(missing_fixed_den, collapse = ", ")), call. = FALSE)
    }
    if (length(valid_fixed_den_regs) == 0) {
      warning(sprintf("Dataset '%s': None of the provided fixed_denominator_regs were found in available all_features Cannot proceed.",
                      dataset_name), call. = FALSE)
      return(NULL) # Fail if NO valid fixed denominator regions remain
    }
    message(sprintf("   - Using %d valid fixed denominator regions.", length(valid_fixed_den_regs)))
  }
  
  # --- Reassign the validated & filtered regions back to the main variables ---
  # These variables will be used subsequently (passed to fitness or used directly)
  fixed_numerator_regs <- valid_fixed_num_regs
  fixed_denominator_regs <- valid_fixed_den_regs
  # --- end of fixed regions prep
  
  
  
  # --- Feature Selection ---
  data_suv <- dataset_data$data_suv_bi
  # if (!id_col %in% names(data_suv)) {
  #   warning(sprintf("Dataset '%s': ID column '%s' not found in data_suv_bi. Cannot run GA.", dataset_name, id_col), call. = FALSE)
  #   return(NULL)
  # }
  if (is.null(features)) {
    potential_features <- setdiff(names(data_suv), id_col)
    is_numeric_col <- sapply(data_suv[, potential_features, drop = FALSE], is.numeric)
    features <- potential_features[is_numeric_col]
    
    if(!is.null(fixed_numerator_regs)) {
      features <- features[!(features %in% fixed_numerator_regs)]
    } else if(!is.null(fixed_denominator_regs)) {
      features <- features[!(features %in% fixed_denominator_regs)]
    }
    
    if (length(features) == 0) {
      warning(sprintf("Dataset '%s': No numeric feature columns found in data_suv_bi. Cannot run GA.", dataset_name), call. = FALSE)
      return(NULL)
    }
    message(sprintf("Dataset '%s': Using %d automatically identified numeric features.", dataset_name, length(features)))
  } else {
    
    if(!is.null(fixed_numerator_regs)) {
      features <- features[!(features %in% fixed_numerator_regs)]
    } else if(!is.null(fixed_denominator_regs)) {
      features <- features[!(features %in% fixed_denominator_regs)]
    }
    
    missing_features <- setdiff(features, names(data_suv))
    if (length(missing_features) > 0) {
      warning(sprintf("Dataset '%s': Specified features not found: %s. Cannot run GA.", dataset_name, paste(missing_features, collapse=", ")), call. = FALSE)
      return(NULL)
    }
    if (length(features) == 0) {
      warning(sprintf("Dataset '%s': No features provided or remaining. Cannot run GA.", dataset_name), call. = FALSE)
      return(NULL)
    }
  }
  n_features <- length(features) # CVR will search among these
  
  # --- GA Bounds & Variables ---
  # *** This depends heavily on your GA type and encoding in .calculate_fitness ***
  # Example assumes "real-valued" with weights for num/den per feature
  ga_type <- ga_config$type %||% "real-valued"
  n_vars <- n_features
  # Add logic for other types if necessary
  
  if(!is.null(fixed_numerator_regs)) {
    min_bounds <- min_bounds %||% 1.1
    max_bounds <- max_bounds %||% 2.9
  } else if(!is.null(fixed_denominator_regs)) {
    min_bounds <- min_bounds %||% 0.1
    max_bounds <- max_bounds %||% 1.9
  } else {
    min_bounds <- min_bounds %||% 0.1
    max_bounds <- max_bounds %||% 2.9
  }
  
  if (ga_type == "real-valued") {
    if (length(min_bounds) == 1) min_bounds <- rep(min_bounds, n_vars)
    if (length(max_bounds) == 1) max_bounds <- rep(max_bounds, n_vars)
    if (length(min_bounds) != n_vars || length(max_bounds) != n_vars) {
      stop(sprintf("GA bounds dimensions (%d, %d) do not match expected n_vars (%d)", length(min_bounds), length(max_bounds), n_vars))
    }
  }
  
  # --- Prepare other shared data ---
  data_clinical <- dataset_data$data
  required_clinical_cols <- c(id_col, "time", "DX", "AB") # Ensure these align with formulas/logic
  if (!all(required_clinical_cols %in% names(data_clinical))) {
    warning(sprintf("Dataset '%s': Clinical data missing required columns: %s. Cannot run GA.",
                    dataset_name, paste(setdiff(required_clinical_cols, names(data_clinical)), collapse=", ")), call. = FALSE)
    return(NULL)
  }
  
  # ============================================================
  # == WRAPPER FITNESS FUNCTION for GA::ga ==
  # ============================================================
  ga_fitness_wrapper <- function(chromosome) {
    # Ensure .calculate_fitness (internal function) is accessible
    # It needs chromosome, features, var_composition, cohort_data, group,
    # eq_all, eq_group, all_power_params, lmer_control
    fitness_value <- tryCatch({
      .calculate_fitness( # Call the internal function
        chromosome = chromosome,
        features = features,
        var_composition = var_composition,
        fixed_numerator_regs,
        fixed_denominator_regs,
        cohort_data = dataset_data, # Pass the single dataset list
        group = group,
        eq_all = eq_all,
        eq_group = eq_group,
        all_power_params = all_power_params,
        lmer_control = lmer_control
      )
    }, error = function(e) {
      warning("Error inside .calculate_fitness: ", e$message, call. = FALSE)
      return(NA_real_) # Return NA on error within the fitness calc
    })
    
    if (is.na(fitness_value) || !is.finite(fitness_value)) {
      return(-Inf) # Return high penalty for maximisation
    } else {
      return(fitness_value)
    }
  }
  # ============================================================
  # == END WRAPPER FITNESS FUNCTION ==
  # ============================================================
  
  
  # --- 7. Run the Genetic Algorithm ---
  message(sprintf("   - Running GA: Dataset=%s, Group=%s, Fitness=Custom(.calculate_fitness), Maximize=TRUE",
                  dataset_name, group))
  
  # Map config params to GA::ga args
  ga_popSize <- ga_config$popSize %||% 32
  ga_maxiter <- ga_config$maxiter %||% 300
  ga_run <- ga_config$run %||% (ga_maxiter)
  ga_pmutation <- ga_config$pmutation %||% 0.9
  ga_pcrossover <- ga_config$pcrossover %||% 0.8
  elitism_raw <- ga_config$elitism_prop %||% 2
  ga_elitism <- max(1, round(elitism_raw))
  if(elitism_raw > 0 && elitism_raw < 1) ga_elitism <- max(1, round(elitism_raw * ga_popSize))
  effective_seed <- ga_seed %||% sample.int(.Machine$integer.max, 1) # Capture the seed being used
  
  ga_args <- list(
    fitness = ga_fitness_wrapper,
    popSize = ga_popSize,
    maxiter = ga_maxiter,
    run = ga_run,
    pmutation = ga_pmutation,
    pcrossover = ga_pcrossover,
    elitism = ga_elitism,
    monitor = T,
    seed = effective_seed
  )
  
  if (ga_type == "real-valued") {
    ga_args$type <- "real-valued"
    ga_args$lower <- min_bounds
    ga_args$upper <- max_bounds
    if(!is.null(ga_config$mutation)) ga_args$mutation <- ga_config$mutation
    if(!is.null(ga_config$crossover)) ga_args$crossover <- ga_config$crossover
  } else {
    stop("Unsupported GA type specified in config: ", ga_type)
  }
  
  # Add extra arguments passed via ...
  extra_args <- list(...)
  ga_args <- c(ga_args, extra_args[!names(extra_args) %in% names(ga_args)])
  
  # --- Execute GA ---
  ga_result_obj <- try(do.call(GA::ga, ga_args), silent = TRUE)
  
  
  
  
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
        dataset_name,
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
                              dataset_name,
                              group,
                              if(!is.null(experiment_tag)) paste0(" - ", experiment_tag) else "" )
        plot(ga_result_obj, main = plot_title)
        
        # Close the device to write the file
        grDevices::dev.off()
        
      }, error = function(e) {
        warning(sprintf("Dataset '%s', Group '%s': Failed to save GA plot to '%s'. Error: %s",
                        dataset_name, group, full_plot_path, conditionMessage(e)), call. = FALSE)
        # Ensure device is closed if error occurred after opening but before dev.off()
        if (grDevices::dev.cur() != 1L) try(grDevices::dev.off(), silent=TRUE) # Use try() for safety
      }) # End tryCatch
    } else {
      warning("Output plot directory '", output_dir, "' does not exist and could not be created. Cannot save plot.", call.=FALSE)
    } # End if dir.exists
  } else if (save_plot && is.null(output_dir)) {
    warning("`save_plot` is TRUE, but `output_dir` was not provided. Plot not saved.", call. = FALSE)
  } # End if save_plot
  
  
  
  
  
  # --- 8. Process and Log Results ---
  if (inherits(ga_result_obj, "try-error")) {
    warning(sprintf("Dataset '%s', Group '%s': GA failed. Error: %s",
                    dataset_name, group, conditionMessage(attr(ga_result_obj, "condition"))), call. = FALSE)
    return(NULL)
  }
  message("   - GA completed.")
  
  # Extract best solution and fitness
  best_solution_vector <- ga_result_obj@solution[1, ]
  actual_best_fitness <- (ga_result_obj@fitnessValue[1])
  
  # --- Re-evaluate the best solution to get Rep, SepAB, SSE ---
  # Decode best_solution_vector into region names. MUST be consistent
  # with the decoding method anticipated by .calculate_fitness.
  if(ga_type == "real-valued"){
    best_regs_numerator <- features[best_solution_vector < 1]
    best_regs_denominator <- features[best_solution_vector > 2]
  } else {
    stop("Decoding not implemented for GA type: ", ga_type)
  }
  
  # Override numerator if fixed_numerator_regs is provided
  if (!is.null(fixed_numerator_regs)) {
    best_regs_numerator <- fixed_numerator_regs
  }
  
  # Override denominator if fixed_denominator_regs is provided
  if (!is.null(fixed_denominator_regs)) {
    best_regs_denominator <- fixed_denominator_regs
  }
  
  if(length(best_regs_numerator)==0 || length(best_regs_denominator)==0){
    warning(sprintf("Dataset '%s', Group '%s': Best GA solution -> empty numerator/denominator.", dataset_name, group), call. = FALSE)
    # You might still want to log this failed result
    # Create a data frame with NAs for metrics but including other info
    final_metrics <- stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
    # Skip further processing, go straight to constructing the output row with NAs
  } else {
    # Calculate final biomarker value with best solution
    # (Use the same logic as your internal .calculate_fitness would use)
    final_value <- tryCatch({
      val <- .calculate_cvr(best_solution_vector, features = features, dataset_cohort_data = dataset_data, var_composition = var_composition,
                            fixed_numerator_regs = fixed_numerator_regs,
                            fixed_denominator_regs = fixed_denominator_regs)
    }, error = function(e) { NULL })
    if (is.null(final_value)){
      warning("Failed final value calculation for best solution.", call.=FALSE)
      return(NULL)
    }
    
    
    # --- Prepare final data for evaluation ---
    dataset_data$data$value <- final_value
    
    # --- Get final metrics using .feval_group ---
    if(nrow(dataset_data$data) < 10 || sum(dataset_data$data$DX == (if(group=="CI") 1L else 0L)) < 5 ) { # Check group size too
      warning(sprintf("Dataset '%s', Group '%s': Not enough final data rows (%d) or group members to calculate metrics reliably.", dataset_name, group, nrow(dataset_data$data)), call.=FALSE)
      final_metrics <- stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
    } else {
      final_metrics <- .feval_group(
        data = dataset_data$data, group = group, eq_all = eq_all, eq_group = eq_group,
        all_power_params = all_power_params, lmer_control = lmer_control
      )
    }
  } # End else block for non-empty regions
  
  
  # --- 9. Construct Result Row ---
  
  # # --- Create prefixed metric names ---
  # # Takes the names from final_metrics (Rep, SepAB, SSE) and adds the group prefix
  # metric_names_prefixed <- paste(names(final_metrics), group, sep = "_") # e.g., "Rep_CI", "SepAB_CI", "SSE_CI"
  # # Create a named list/vector suitable for data.frame construction
  # metrics_for_df <- as.list(stats::setNames(final_metrics, metric_names_prefixed))
  
  
  # --- Define composition (Use input var_composition as requested) ---
  # No calculation needed here, just log the input parameter
  composition_value <- var_composition
  
  
  result_row <- data.frame(
    # === Metadata ===
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    experiment_tag = experiment_tag %||% NA_character_,
    discovery_dataset = dataset_name, # Renamed from 'dataset'
    group_evaluated = group,          # Group used for fitness/metrics
    bilateral = bilateral,
    var_composition = composition_value, # Logged input parameter
    
    # === Fitness Info ===
    # fitness_metric_used = "custom", # Less informative now we know it's .calculate_fitness
    # maximized_fitness = TRUE,      # Implied by design
    fitness_value = actual_best_fitness, # The direct output of .calculate_fitness
    
    # === Evaluation Metrics (Prefixed) ===
    # This adds columns like Rep_CI, SepAB_CI, SSE_CI etc. dynamically
    # metrics_for_df,
    
    as.list(final_metrics),
    
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
  
  # --- 10. Append to CSV (Optional) ---
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
  
  # --- 11. Return Result Row ---
  return(result_row)
}


