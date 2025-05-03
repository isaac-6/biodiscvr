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
#' @param features Character vector. Names of columns in `dataset_data$data_suv_bi`
#'   to be considered as features by the GA. If NULL (default), uses all numeric
#'   columns excluding the ID column specified in `config$preprocessing$id_column`.
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
#'   a unique seed derived from this base (`base_ga_seed`) will be
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
#' @importFrom progress progress_bar
run_experiments <- function(prepared_data_list, 
                            config, 
                            groups = c("CU", "CI"), 
                            experiments_config_path, 
                            features = NULL, 
                            output_csv_name, 
                            output_dir = NULL, 
                            save_plots = TRUE, 
                            datasets_to_run = names(prepared_data_list), 
                            experiment_master_tag = format(Sys.time(), "%Y%m%d"), 
                            base_ga_seed = 42) {
  
  # --- Initial Validations ---
  stopifnot(is.list(prepared_data_list), length(prepared_data_list) > 0, is.list(config),
            rlang::is_scalar_character(experiments_config_path),
            is.null(features) | is.character(features),
            rlang::is_scalar_character(output_csv_name),
            is.null(output_dir) | rlang::is_scalar_character(output_dir),
            rlang::is_scalar_logical(save_plots),
            is.character(datasets_to_run), length(datasets_to_run) > 0,
            rlang::is_scalar_character(experiment_master_tag),
            is.null(base_ga_seed) | rlang::is_scalar_double(base_ga_seed))
  
  # Validate datasets
  missing_dsets <- setdiff(datasets_to_run, names(prepared_data_list))
  if (length(missing_dsets) > 0) stop("Datasets not found: ", paste(missing_dsets, collapse=", "))
  datasets_to_run <- intersect(datasets_to_run, names(prepared_data_list))
  
  # Validate config
  req_conf_sections <- c("preprocessing", "model_equations", "power_params", "genetic_algorithm")
  if (!all(req_conf_sections %in% names(config))) stop("Config missing required sections.")
  
  # --- Setup Common Variables ---
  id_col <- config$preprocessing$id_column %||% "RID"
  ignore_cols <- config$preprocessing$ignore_columns %||% "idx"
  
  var_composition <- config$genetic_algorithm$var_composition %||% 1
  if (!var_composition %in% 0:3) stop("Invalid var_composition in config.")
  
  # Load experiments from YAML
  if (!file.exists(experiments_config_path)) stop("Experiments config file not found: ", experiments_config_path)
  experiments_list <- yaml::read_yaml(experiments_config_path)$experiments
  if (!is.list(experiments_list) || length(experiments_list) == 0) stop("No experiments defined in config file.")
  
  message(sprintf("Loaded %d experiment definitions.", length(experiments_list)))
  
  # --- Identify Common Features Across Datasets ---
  message("Identifying common features across datasets...")
  common_features <- Reduce(intersect, lapply(datasets_to_run, function(dset) {
    data_suv <- prepared_data_list[[dset]]$data_suv_bi
    if (is.null(data_suv)) return(NULL)
    setdiff(names(data_suv), c(id_col, ignore_cols))
  }))
  
  if (!is.null(features)) {
    common_features <- intersect(common_features, features)
  }
  
  if (length(common_features) == 0) stop("No common features identified across datasets.")
  
  message(sprintf("Found %d common features.", length(common_features)))
  
  
  # --- progress bar ---
  # Multicohort take as many time as the number of cohorts
  # Experiments run twice if there is a second iteration (for fixed regions)
  num_2iter <- sum(sapply(experiments_list, function(x) x$run_second_iteration))
  aux <- ifelse(length(datasets_to_run) > 1, length(datasets_to_run) +1, 1)
  total_iterations <- aux * (length(experiments_list)+num_2iter) * length(groups)
  
  pb <- progress_bar$new(
    format = "Progress [:bar] :percent (:current/:total) | ETA: :eta \n",
    total = total_iterations,
    clear = FALSE,
    width = 50
  )
  
  
  
  # --- Single-Cohort Experiment Runs ---
  single_cohort_fitness_store <- list()
  run_counter <- 0
  
  for (exp_name in names(experiments_list)) {
    exp_def <- experiments_list[[exp_name]]
    message(sprintf("Running Single-Cohort Experiment: %s", exp_name))

    fixed_num_iter1 <- exp_def$fixed_numerator_regs
    fixed_den_iter1 <- exp_def$fixed_denominator_regs
    run_second_iter <- isTRUE(exp_def$run_second_iteration)

    first_iter_regions <- list()

    for (group in groups) {
      single_cohort_fitness_store[[exp_name]][[group]] <- list()

      for (dset_name in datasets_to_run) {
        run_counter <- run_counter + 1
        current_tag <- paste(experiment_master_tag, exp_name, group, "Iter1", sep="_")

        message(sprintf("Running Iter1: %s, %s", dset_name, group))

        run_result <- try(biodiscvr_single(
          dataset_data = prepared_data_list[[dset_name]],
          dataset_name = dset_name,
          group = group,
          config = config,
          features = common_features,
          var_composition = var_composition,
          fixed_numerator_regs = fixed_num_iter1,
          fixed_denominator_regs = fixed_den_iter1,
          experiment_tag = current_tag,
          output_csv_name = output_csv_name,
          output_dir = output_dir,
          save_plot = save_plots,
          ga_seed = base_ga_seed
        ), silent = TRUE)

        if (!inherits(run_result, "try-error") && !is.null(run_result$result_row)) {
          fitness_val <- run_result$result_row$fitness_value
          single_cohort_fitness_store[[exp_name]][[group]][[dset_name]] <- ifelse(is.finite(fitness_val), fitness_val, NA)

          if (run_second_iter) {
            first_iter_regions[[paste(dset_name, group, sep="_")]] <- list(
              num = run_result$best_regs_numerator,
              den = run_result$best_regs_denominator
            )
          }
        }
        
        pb$tick()
      }
    }

    # --- Second Iteration (if required) ---
    if (run_second_iter) {
      message(sprintf("Running Second Iteration for %s", exp_name))

      for (group in groups) {
        for (dset_name in datasets_to_run) {
          result_key <- paste(dset_name, group, sep="_")
          first_iter_res <- first_iter_regions[[result_key]]
          if (is.null(first_iter_res)) next

          new_fixed_num <- if (!is.null(fixed_den_iter1)) first_iter_res$num else NULL
          new_fixed_den <- if (!is.null(fixed_num_iter1)) first_iter_res$den else NULL
          iter2_tag <- paste(experiment_master_tag, exp_name, group, "Iter2", sep="_")

          run_result_iter2 <- try(biodiscvr_single(
            dataset_data = prepared_data_list[[dset_name]],
            dataset_name = dset_name,
            group = group,
            config = config,
            features = common_features,
            fixed_numerator_regs = new_fixed_num,
            fixed_denominator_regs = new_fixed_den,
            experiment_tag = iter2_tag,
            output_csv_name = output_csv_name,
            output_dir = output_dir,
            save_plot = save_plots,
            ga_seed = base_ga_seed
          ), silent = TRUE)
          
          pb$tick()
        }
      }
    }
  }
  
  # --- Multi-Cohort Phase (with Second Iteration) ---
  message("Starting Multi-Cohort Experiment Runs...")
  
  for (exp_name in names(experiments_list)) {
    exp_def <- experiments_list[[exp_name]]
    
    for (group in groups) {
      ref_fitness_vector <- unlist(single_cohort_fitness_store[[exp_name]][[group]])
      ref_fitness_vector <- ref_fitness_vector[!is.na(ref_fitness_vector)] # Remove NAs
      
      if (length(ref_fitness_vector) < 2) {
        warning(sprintf("Skipping Multi-Cohort run for '%s' - Group '%s': Not enough valid fitness values.", exp_name, group), call. = FALSE)
        next
      }
      
      # --- First Multi-Cohort Iteration ---
      first_iter_tag <- paste(experiment_master_tag, exp_name, group, "MultiCohort_Iter1", sep = "_")
      
      message(sprintf("Running First MultiCohort Iteration: %s (%d datasets)", first_iter_tag, length(ref_fitness_vector)))
      
      first_iter_result <- try(biodiscvr_multicohort(
        preprocessed_data = prepared_data_list,
        datasets_to_run = datasets_to_run,
        group = group,
        config = config,
        features = common_features,
        fixed_numerator_regs = exp_def$fixed_numerator_regs,
        fixed_denominator_regs = exp_def$fixed_denominator_regs,
        reference_fitness = ref_fitness_vector,
        experiment_tag = first_iter_tag,
        output_csv_name = output_csv_name,
        output_dir = output_dir,
        save_plot = save_plots,
        ga_seed = base_ga_seed
      ), silent = TRUE)
      
      if (inherits(first_iter_result, "try-error") || is.null(first_iter_result)) {
        warning(sprintf("Multi-Cohort Iter1 failed for '%s' - Group '%s'.", exp_name, group), call. = FALSE)
        next
      }
      
      
      # it takes more time to run the multicohort, but we count it as one iteration
      pb$tick()
      
      # Store best numerator and denominator from first multi-cohort iteration
      best_num_regions <- first_iter_result$best_regs_numerator
      best_den_regions <- first_iter_result$best_regs_denominator
      
      # --- Determine Fixed Regions for Second Iteration ---
      new_fixed_num <- if (!is.null(exp_def$fixed_denominator_regs)) best_num_regions else NULL
      new_fixed_den <- if (!is.null(exp_def$fixed_numerator_regs)) best_den_regions else NULL
      
      message(sprintf("DEBUG: Running second multi-cohort iteration for %s - Group %s (Numerator Fixed: %s, Denominator Fixed: %s)", 
                      exp_name, group, 
                      !is.null(new_fixed_num), 
                      !is.null(new_fixed_den)))
      
      if (!is.null(new_fixed_num) || !is.null(new_fixed_den)) {
        second_iter_tag <- paste(experiment_master_tag, exp_name, group, "MultiCohort_Iter2", sep = "_")
        
        message(sprintf("Running Second MultiCohort Iteration: %s", second_iter_tag))
        
        second_iter_result <- try(biodiscvr_multicohort(
          preprocessed_data = prepared_data_list,
          datasets_to_run = datasets_to_run,
          group = group,
          config = config,
          features = common_features,
          fixed_numerator_regs = new_fixed_num,
          fixed_denominator_regs = new_fixed_den,
          reference_fitness = ref_fitness_vector,
          experiment_tag = second_iter_tag,
          output_csv_name = output_csv_name,
          output_dir = output_dir,
          save_plot = save_plots,
          ga_seed = base_ga_seed
        ), silent = TRUE)
        
        if (inherits(second_iter_result, "try-error") || is.null(second_iter_result)) {
          warning(sprintf("Multi-Cohort Iter2 failed for '%s' - Group '%s'.", exp_name, group), call. = FALSE)
        }
        
        # it takes more time to run the multicohort, but we count it as one iteration
        pb$tick()
      }
    }
  }
  
  message("Multi-Cohort experiments completed.")
  
  message("All experiments completed.")
  message("Results logged to: ", output_csv_name)
  
  invisible(list(
    output_csv_name = output_csv_name,
    output_dir = output_dir
    # single_cohort_fitness = single_cohort_fitness_store # Return stored fitness
  ))
}