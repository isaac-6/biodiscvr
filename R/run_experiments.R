#' Run a Defined Set of Single- and Multi-Cohort Discovery Experiments
#'
#' Orchestrates multiple biomarker discovery runs using specified single-cohort
#' (`biodiscvr_single`) and multi-cohort (`biodiscvr_multicohort`) functions.
#' Executes a sequence of experiments defined in a configuration file, including
#' baseline discovery, runs with predefined fixed regions, optional second-iteration
#' runs based on initial fixed-region results, and finally multi-cohort optimization
#' using single-cohort results as reference. All primary results are appended to
#' a single CSV file.
#'
#' @param prepared_data_list List. The output from `preprocess_datasets()`,
#'   containing filtered data for multiple datasets (e.g., `list(ADNI=list(data=..., data_suv_bi=...), MAYO=...)`).
#' @param config List. The loaded main configuration object (output from
#'   `check_and_prepare_data$config`). Must contain sections like
#'   `preprocessing`, `model_equations`, `power_params`, `genetic_algorithm`.
#' @param groups groups to run the experiments for. Defaults to c("CU", "CI")
#' @param experiments_config_path Character string. Path to the YAML file defining
#'   the experiments to run. See Details for expected structure.
#' @param output_csv_name Character string. Full path to the CSV file where
#'   results from ALL experimental runs (single and multi-cohort) will be appended.
#'   The file will be created with headers if it doesn't exist. Directory created if needed.
#' @param output_dir Character string or NULL. Path to the directory where
#'   GA convergence plots from *single-cohort* runs should be saved.
#'   If NULL, plots are not saved. Directory created if needed.
#' @param save_plots Logical. If `TRUE` (and `output_dir` is provided),
#'   save the GA convergence plot for each single-cohort run. Defaults to `TRUE`.
#'   (Note: Multi-cohort plot saving might need separate handling if desired).
#' @param datasets_to_run Character vector. Names of the datasets within
#'   `prepared_data_list` to include in the experiments (e.g., `c("ADNI", "MAYO")`).
#'   Defaults to using all datasets present in `prepared_data_list`.
#' @param experiment_master_tag Character string. An overall tag prepended to
#'   individual run tags for easier filtering of results. Defaults to the
#'   current date (YYYYMMDD).
#' @param base_ga_seed Numeric or NULL. A base seed for the GA runs. If provided,
#'   a unique seed derived from this base (`base_ga_seed + run_counter`) will be
#'   used for each GA run for reproducibility. If NULL, random seeds are used.
#'
#' @return Invisibly returns a list containing the path to the main output CSV file
#'   and potentially other summary information or paths generated during the run.
#'   The primary output is the side effect of appending results to the CSV file.
#'
#' @details
#'   **Workflow:**
#'   1.  Loads experiment definitions from `experiments_config_path`.
#'   2.  Identifies features common across all specified `datasets_to_run`.
#'   3.  Validates any fixed regions defined in experiments against common features.
#'   4.  **Single-Cohort Phase:** Iterates through each experiment definition:
#'       *   Runs `biodiscvr_single` for CU and CI groups on each dataset.
#'       *   If an experiment uses fixed regions and `run_second_iteration` is TRUE,
#'         it immediately runs a second iteration fixing the regions found in the
#'         first iteration (fixing numerator if denominator was initially fixed, and vice-versa).
#'       *   Appends results of each `biodiscvr_single` run to `output_csv_name`.
#'       *   Stores the single-cohort fitness values achieved in the first iteration runs.
#'   5.  **Multi-Cohort Phase:** Iterates through each experiment definition again:
#'       *   Retrieves the stored single-cohort fitness values for the current experiment/group
#'         to be used as the `reference_fitness` vector.
#'       *   Calls `biodiscvr_multicohort` to perform optimization across all
#'         `datasets_to_run` simultaneously, using the reference fitness.
#'       *   Appends the main result row from `biodiscvr_multicohort` to the *same*
#'         `output_csv_name`.
#'       *   (Note: Per-dataset metrics from multi-cohort are returned by the function
#'         but not automatically logged to the main CSV in this implementation).
#'
#'   **Experiments Configuration (`experiments.yaml`):**
#'   Should be a YAML file with a top-level key `experiments`. Each key under
#'   `experiments` defines a unique run configuration:
#'   ```yaml
#'   experiments:
#'     ExperimentName1:
#'       description: "A brief description"
#'       fixed_numerator_regs: NULL # or list ["RegionA", "RegionB"]
#'       fixed_denominator_regs: NULL # or list ["RegionC"]
#'       run_second_iteration: FALSE # or TRUE (only relevant if one side is fixed)
#'     ExperimentName2:
#'       # ... other definitions ...
#'   ```
#'
#' @export
#' @importFrom utils packageName packageVersion sessionInfo head str
#' @importFrom yaml read_yaml
#' @importFrom rlang `%||%` is_character
#' @importFrom dplyr bind_rows
#' @importFrom stats setNames
run_experiments <- function(prepared_data_list,
                            config,
                            groups = c("CU", "CI"),
                            experiments_config_path,
                            output_csv_name,
                            output_dir = NULL,
                            save_plots = TRUE,
                            datasets_to_run = names(prepared_data_list),
                            experiment_master_tag = format(Sys.time(), "%Y%m%d"),
                            base_ga_seed = 42) {
  
  # --- Initial Validations & Setup ---
  stopifnot(
    is.list(prepared_data_list), length(prepared_data_list) > 0,
    is.list(config),
    rlang::is_scalar_character(experiments_config_path),
    rlang::is_scalar_character(output_csv_name),
    is.null(output_dir) || rlang::is_scalar_character(output_dir),
    rlang::is_scalar_logical(save_plots),
    is.character(datasets_to_run), length(datasets_to_run) > 0,
    rlang::is_scalar_character(experiment_master_tag),
    is.null(base_ga_seed) || rlang::is_scalar_double(base_ga_seed)
  )
  # Check datasets
  missing_dsets <- setdiff(datasets_to_run, names(prepared_data_list))
  if (length(missing_dsets) > 0) stop("Datasets not found: ", paste(missing_dsets, collapse=", "))
  # Check config
  req_conf_sections <- c("preprocessing", "model_equations", "power_params", "genetic_algorithm")
  if (!all(req_conf_sections %in% names(config))) stop("Config missing required sections.")
  
  # --- Prepare Common Elements ---
  id_col <- config$preprocessing$id_column %||% "RID"
  all_power_params <- config$power_params
  var_composition <- config$genetic_algorithm$var_composition %||% 1
  if (! var_composition %in% 0:3) stop("Invalid var_composition in config.")
  # Formulas
  eq_all_str <- config$model_equations$eq_all_string
  eq_group_str <- config$model_equations$eq_group_string
  if(is.null(eq_all_str) || is.null(eq_group_str)) stop("Missing formula strings.")
  eq_all <- try(stats::as.formula(eq_all_str)); if(inherits(eq_all, "try-error")) stop("Invalid eq_all_str")
  eq_group <- try(stats::as.formula(eq_group_str)); if(inherits(eq_group, "try-error")) stop("Invalid eq_group_str")
  # LMER Control
  lmer_control <- lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                                    check.conv.singular = "ignore",
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
  
  # --- Load Experiments Configuration ---
  if (!file.exists(experiments_config_path)) stop("Experiments config file not found: ", experiments_config_path)
  experiments_yaml <- try(yaml::read_yaml(experiments_config_path))
  if (inherits(experiments_yaml, "try-error") || !is.list(experiments_yaml) || !"experiments" %in% names(experiments_yaml)) {
    stop("Failed to load or parse experiments from: ", experiments_config_path)
  }
  experiments_list <- experiments_yaml$experiments
  if(!is.list(experiments_list) || length(experiments_list) == 0) stop("No experiments defined in config file.")
  message(sprintf("Loaded %d experiment definitions.", length(experiments_list)))
  
  # --- Find Common Features & Validate Predefined Sets ---
  message("--- Preparing for Experiments ---")
  common_features <- NULL
  processed_at_least_one <- FALSE
  for (dset_name in datasets_to_run) {
    dset_data_list <- prepared_data_list[[dset_name]]
    if (!"data_suv_bi" %in% names(dset_data_list) || !is.data.frame(dset_data_list$data_suv_bi)) next
    current_features_all <- names(dset_data_list$data_suv_bi)
    current_features_valid <- setdiff(current_features_all, id_col)
    if (length(current_features_valid) == 0) next
    if (!processed_at_least_one) { common_features <- current_features_valid; processed_at_least_one <- TRUE }
    else { common_features <- intersect(common_features, current_features_valid) }
    if (length(common_features) == 0) break
  }
  if (!processed_at_least_one || length(common_features) == 0) {
    stop("Could not identify any common features across the specified datasets. Cannot proceed.")
  }
  message(sprintf("Identified %d common features for analysis.", length(common_features)))
  
  # Validate predefined fixed sets defined within experiments_list
  validated_experiments <- list()
  for (exp_name in names(experiments_list)) {
    exp_def <- experiments_list[[exp_name]]
    valid_exp_def <- exp_def # Start with original definition
    valid_exp_def$fixed_numerator_regs_validated <- NULL
    valid_exp_def$fixed_denominator_regs_validated <- NULL
    valid_exp <- TRUE # Assume valid initially
    
    if (!is.null(exp_def$fixed_numerator_regs)) {
      valid_num <- intersect(exp_def$fixed_numerator_regs, common_features)
      if (length(valid_num) != length(exp_def$fixed_numerator_regs)) warning(sprintf("Exp '%s': Some fixed num regions excluded.", exp_name), call.=FALSE)
      if (length(valid_num) == 0) { warning(sprintf("Exp '%s': No valid fixed num regions remain. Skipping.", exp_name), call.=FALSE); valid_exp <- FALSE }
      else { valid_exp_def$fixed_numerator_regs_validated <- valid_num }
    }
    if (valid_exp && !is.null(exp_def$fixed_denominator_regs)) {
      valid_den <- intersect(exp_def$fixed_denominator_regs, common_features)
      if (length(valid_den) != length(exp_def$fixed_denominator_regs)) warning(sprintf("Exp '%s': Some fixed den regions excluded.", exp_name), call.=FALSE)
      if (length(valid_den) == 0) { warning(sprintf("Exp '%s': No valid fixed den regions remain. Skipping.", exp_name), call.=FALSE); valid_exp <- FALSE }
      else { valid_exp_def$fixed_denominator_regs_validated <- valid_den }
    }
    
    if(valid_exp) {
      validated_experiments[[exp_name]] <- valid_exp_def
    }
  }
  if(length(validated_experiments) == 0) stop("No valid experiments remaining after validating fixed regions.")
  message(sprintf("Validated %d experiments.", length(validated_experiments)))
  
  
  # --- Store Single-Cohort Results for Multi-Cohort Phase ---
  single_cohort_fitness_store <- list() # Structure: $ExpName$GroupName$DatasetName -> fitness
  
  # --- Experiment Sequence ---
  message("\n--- Starting Single-Cohort Experiment Runs ---")
  run_counter <- 0
  
  # == Iterate Through Validated Experiments ==
  for (exp_name in names(validated_experiments)) {
    exp_def <- validated_experiments[[exp_name]]
    message(sprintf("\n*** Running Experiment: %s (%s) ***", exp_name, exp_def$description %||% "No description"))
    
    # Use the *validated* fixed regions for this experiment
    fixed_num_iter1 <- exp_def$fixed_numerator_regs_validated # Can be NULL
    fixed_den_iter1 <- exp_def$fixed_denominator_regs_validated # Can be NULL
    run_second_iter <- isTRUE(exp_def$run_second_iteration)
    
    # Store results from the first iteration if second iteration is needed
    first_iter_regions <- list() # Store list(num=..., den=...) keyed by Dataset_Group
    
    # --- Run First Iteration (Baseline or Predefined Fixed) ---
    for (group in groups) {
      single_cohort_fitness_store[[exp_name]] <- single_cohort_fitness_store[[exp_name]] %||% list()
      single_cohort_fitness_store[[exp_name]][[group]] <- list()
      
      for (dset_name in datasets_to_run) {
        run_counter <- run_counter + 1
        current_seed <- if(!is.null(base_ga_seed)) base_ga_seed else 42
        current_tag <- paste(experiment_master_tag, exp_name, group, "Iter1", sep = "_")
        
        message(sprintf("Running Iter1: Dataset=%s, Group=%s, Tag=%s", dset_name, group, current_tag))
        run_result_list <- try(biodiscvr_single(
          dataset_data = prepared_data_list[[dset_name]],
          dataset_name = dset_name, 
          group = group, 
          config = config,
          features = common_features, # Use common features
          var_composition = var_composition,
          fixed_numerator_regs = fixed_num_iter1, # Use validated fixed for Iter1
          fixed_denominator_regs = fixed_den_iter1,
          experiment_tag = current_tag,
          output_csv_name = output_csv_name,
          output_dir = output_dir, 
          save_plot = save_plots,
          ga_seed = current_seed,
          # Pass bounds from config if needed by biodiscvr_single
          min_bounds = config$genetic_algorithm$min_bounds,
          max_bounds = config$genetic_algorithm$max_bounds
        ), silent = TRUE)
        
        
        # --- Store Fitness & Regions ---
        if (!inherits(run_result_list, "try-error") && is.list(run_result_list) && !is.null(run_result_list$result_row)) {
          fitness_val <- run_result_list$result_row$fitness_value
          single_cohort_fitness_store[[exp_name]][[group]][[dset_name]] <- ifelse(is.finite(fitness_val), fitness_val, NA) # Store NA if Inf/-Inf
          
          if (run_second_iter && (!is.null(fixed_num_iter1) || !is.null(fixed_den_iter1))) { # Only store if fixed & iter2 needed
            result_key <- paste(dset_name, group, sep="_")
            first_iter_regions[[result_key]] <- list(
              num = run_result_list$best_regs_numerator,
              den = run_result_list$best_regs_denominator
            )
          }
        } else {
          warning(sprintf("Run failed/returned invalid result: Dataset=%s, Group=%s, Tag=%s.", dset_name, group, current_tag), call.=FALSE)
          single_cohort_fitness_store[[exp_name]][[group]][[dset_name]] <- NA
        }
      } # end dataset loop (Iter1)
    } # end group loop (Iter1)
    
    
    # --- Run Second Iteration (if requested and applicable) ---
    if (run_second_iter && (!is.null(fixed_num_iter1) || !is.null(fixed_den_iter1)) && is.null(fixed_num_iter1) != is.null(fixed_den_iter1) ) { # Only run if exactly ONE side was fixed
      message(sprintf("\n*** Running Second Iteration for Experiment: %s ***", exp_name))
      
      for (group in groups) {
        for (dset_name in datasets_to_run) {
          result_key <- paste(dset_name, group, sep="_")
          first_iter_res <- first_iter_regions[[result_key]]
          
          if (is.null(first_iter_res)) {
            message(sprintf("Skipping Iter2 for %s/%s: No valid result from Iter1.", dset_name, group)); next
          }
          
          new_fixed_num <- NULL; new_fixed_den <- NULL; iter2_tag_suffix <- ""
          if (!is.null(fixed_den_iter1)) { # Original run fixed Denominator
            new_fixed_num <- first_iter_res$num; iter2_tag_suffix <- "FixNumFound"
            if(is.null(new_fixed_num) || length(new_fixed_num)==0) {message(sprintf("Skipping Iter2 for %s/%s: No numerator found in Iter1.",dset_name,group)); next}
          } else { # Original run fixed Numerator
            new_fixed_den <- first_iter_res$den; iter2_tag_suffix <- "FixDenFound"
            if(is.null(new_fixed_den) || length(new_fixed_den)==0) {message(sprintf("Skipping Iter2 for %s/%s: No denominator found in Iter1.",dset_name,group)); next}
          }
          
          run_counter <- run_counter + 1
          current_seed <- if(!is.null(base_ga_seed)) base_ga_seed else NULL
          current_tag <- paste(experiment_master_tag, exp_name, group, iter2_tag_suffix, sep = "_")
          
          message(sprintf("Running Iter2: Dataset=%s, Group=%s, Tag=%s", dset_name, group, current_tag))
          run_result_list_iter2 <- try(biodiscvr_single(
            dataset_data = prepared_data_list[[dset_name]],
            dataset_name = dset_name, group = group, config = config,
            features = common_features,
            var_composition = var_composition,
            fixed_numerator_regs = new_fixed_num, # Use newly fixed regions
            fixed_denominator_regs = new_fixed_den,
            experiment_tag = current_tag,
            output_csv_name = output_csv_name,
            output_dir = output_dir, 
            save_plot = save_plots,
            ga_seed = current_seed,
            min_bounds = config$genetic_algorithm$min_bounds,
            max_bounds = config$genetic_algorithm$max_bounds
          ), silent = TRUE)
          
          if (inherits(run_result_list_iter2, "try-error") || is.null(run_result_list_iter2)) {
            warning(sprintf("Run failed/returned NULL: Dataset=%s, Group=%s, Tag=%s.", dset_name, group, current_tag), call.=FALSE)
          }
          # No need to store results from Iter2 for multi-cohort phase
        } # end dataset loop (Iter2)
      } # end group loop (Iter2)
    } else if (run_second_iter) {
      message(sprintf("Skipping Second Iteration for Exp '%s': Not applicable (requires exactly one side fixed in Iter1).", exp_name))
    }# end if run_second_iter
    
  } # end loop through experiments
  
  
  # ==================================================
  # == Multi-Cohort Phase ==
  # ==================================================
  message("\n--- Starting Multi-Cohort Experiment Runs ---")
  
  # Check if the multi-cohort function exists (replace with actual name if different)
  if (!exists("biodiscvr_multicohort") || !is.function(biodiscvr_multicohort)){
    warning("Function 'biodiscvr_multicohort' is not defined or available. Skipping multi-cohort phase.", call.=FALSE)
  } else {
    for (exp_name in names(validated_experiments)) { # Use validated experiments
      exp_def <- validated_experiments[[exp_name]]
      message(sprintf("\n*** Running Multi-Cohort for Experiment: %s ***", exp_name))
      
      # --- Prepare fixed regions for this experiment ---
      fixed_num <- exp_def$fixed_numerator_regs_validated # Use validated
      fixed_den <- exp_def$fixed_denominator_regs_validated # Use validated
      
      # --- Run multi-cohort optimization for EACH group ---
      for (group in groups) {
        
        # --- Fetch the reference fitness vector ---
        ref_fitness_list <- single_cohort_fitness_store[[exp_name]][[group]]
        # Filter for datasets actually being run and remove NAs
        ref_fitness_list <- ref_fitness_list[names(ref_fitness_list) %in% datasets_to_run]
        ref_fitness_vector <- unlist(ref_fitness_list)
        ref_fitness_vector <- ref_fitness_vector[!is.na(ref_fitness_vector)]
        
        if(length(ref_fitness_vector) == 0) {
          warning(sprintf("Skipping Multi-Cohort run for Exp '%s', Group '%s': No valid single-cohort reference fitness values found for target datasets.", exp_name, group), call.=FALSE)
          next
        }
        # Ensure names match datasets_to_run included in the vector
        datasets_for_mc_ref <- names(ref_fitness_vector)
        
        
        run_counter <- run_counter + 1
        current_seed <- if(!is.null(base_ga_seed)) base_ga_seed else NULL
        current_tag <- paste(experiment_master_tag, exp_name, group, "MultiCohort", sep = "_")
        
        message(sprintf("Running MultiCohort: Exp=%s, Group=%s, Tag=%s (on %d datasets)", exp_name, group, current_tag, length(datasets_for_mc_ref)))
        
        # --- Call the multi-cohort function ---
        mc_result_list <- try(biodiscvr_multicohort(
          preprocessed_data = prepared_data_list, # Pass the main list
          datasets_to_run = datasets_for_mc_ref, # Use datasets with valid reference
          group = group,
          config = config,
          features = common_features, # Pass common features
          fixed_numerator_regs = fixed_num,
          fixed_denominator_regs = fixed_den,
          var_composition = var_composition,
          reference_fitness = ref_fitness_vector, # Pass the named vector
          bilateral = TRUE, # Get from config?
          experiment_tag = current_tag,
          output_csv_name = output_csv_name,
          output_dir = output_dir, 
          save_plot = save_plots, 
          ga_seed = current_seed,
          # Pass bounds from config if needed by biodiscvr_multicohort
          min_bounds = config$genetic_algorithm$min_bounds,
          max_bounds = config$genetic_algorithm$max_bounds
          # Add other necessary args from config$genetic_algorithm
        ), silent = TRUE)
        
        if (inherits(mc_result_list, "try-error")) {
          warning(sprintf("Multi-Cohort run failed: Exp=%s, Group=%s, Tag=%s. Error: %s", exp_name, group, current_tag, conditionMessage(attr(mc_result_list,"condition"))), call.=FALSE)
        } else {
          # Optionally process/log per_dataset_metrics from mc_result_list if needed
          message(sprintf("Multi-Cohort run completed for Exp '%s', Group '%s'.", exp_name, group))
        }
      } # end group loop for multi-cohort
    } # end experiment loop for multi-cohort
  } # end else block (if biodiscvr_multicohort exists)
  
  
  message("\n--- All Experiments Completed ---")
  message("Results logged to: ", output_csv_name)
  if(save_plots && !is.null(output_dir)) message("Single-cohort plots saved to: ", output_dir)
  
  invisible(list(
    output_csv_name = output_csv_name,
    output_dir = output_dir,
    single_cohort_fitness = single_cohort_fitness_store # Return stored fitness
  ))
}