#' Evaluate Predefined Literature Biomarkers
#'
#' Calculates and evaluates the performance (Rep, SepAB, SSE) of a predefined
#' set of literature-based biomarkers across specified datasets and groups.
#'
#' @param lit_biomarker_file_path Path to literature biomarkers file. A named list where each element defines a
#'   literature biomarker. Each element must be a list containing at least:
#'   `$name` (character), `$numerator_regions` (character vector),
#'   `$denominator_regions` (character vector), and `$var_composition` (numeric, 0-3, or -1 if calculation is not traditional SUVR).
#'   Optionally include `$description` and `$source_citation`.
#' @param prepared_data_list List. The output from `preprocess_data()`,
#'   containing filtered data for multiple datasets.
#' @param config List. The loaded main configuration object.
#' @param datasets_to_evaluate Character vector. Names of the datasets within
#'   `prepared_data_list` to evaluate the biomarkers on. Defaults to using all.
#' @param groups_to_evaluate Character vector. Groups ("CU", "CI") to evaluate.
#'   Defaults to `c("CU", "CI")`.
#' @param calculate_ci Logical. Calculate 95% CIs? Defaults to `FALSE`.
#' @param nsim Integer. Number of bootstrap simulations if `calculate_ci` is `TRUE`.
#'   Defaults to 100.
#' @param output_csv_path Character string or NULL. Path to save the evaluation
#'   results as a *new* CSV file (overwrites if exists).
#' @param id_col Character string. Name of the unique ID column. Defaults to "RID".
#' @param verbose Logical. Print progress messages? Defaults to TRUE.
#'
#' @return A data frame containing the evaluation results. Each row corresponds to
#'   one literature biomarker evaluated on one specific dataset and group.
#'   Returns `NULL` if essential data is missing or no biomarkers can be evaluated.
#'
#' @details
#' Iterates through each literature biomarker definition and each specified
#' dataset/group. For each combination, it checks if the required regions exist
#' in the dataset, calculates the composite biomarker value using a simplified
#' internal helper (`.calculate_value_from_regions`), prepares the data, checks
#' for sufficient samples, and calls the appropriate evaluation function
#' (`.feval_group` or `.feval_group_95CI`). Results are aggregated and optionally saved.
#'
#' @export
#' @importFrom dplyr bind_rows select all_of left_join filter n_distinct between
#' @importFrom rlang `%||%` .data sym is_scalar_character is_scalar_logical is_scalar_integerish
#' @importFrom stats setNames median quantile IQR sd na.omit aggregate as.formula
#' @importFrom readr write_csv locale read_csv 
#' @importFrom lme4 lmerControl
#' @importFrom longpower lmmpower
#' @importFrom methods is
evaluate_literature_biomarkers <- function(lit_biomarker_file_path = system.file("files", "literature_biomarkers.yaml", package = "biodiscvr"),
                                           prepared_data_list,
                                           config,
                                           datasets_to_evaluate = names(prepared_data_list),
                                           groups_to_evaluate = c("CU", "CI"),
                                           calculate_ci = FALSE,
                                           nsim = 1000,
                                           output_csv_path = NULL,
                                           id_col = NULL,
                                           verbose = TRUE) {
  
  # --- Input Validation ---
  stopifnot(
    is.null(lit_biomarker_file_path) || rlang::is_scalar_character(lit_biomarker_file_path),
    is.list(prepared_data_list), length(prepared_data_list) > 0,
    is.list(config),
    is.character(datasets_to_evaluate), length(datasets_to_evaluate) > 0,
    is.character(groups_to_evaluate), all(groups_to_evaluate %in% c("CU", "CI")),
    rlang::is_scalar_logical(calculate_ci),
    rlang::is_scalar_integerish(nsim), nsim > 1,
    is.null(output_csv_path) || rlang::is_scalar_character(output_csv_path),
    is.null(id_col) || rlang::is_scalar_character(id_col),
    rlang::is_scalar_logical(verbose)
  )
  # Check datasets exist
  missing_dsets <- setdiff(datasets_to_evaluate, names(prepared_data_list))
  if (length(missing_dsets) > 0) stop("Datasets not found: ", paste(missing_dsets, collapse=", "))
  # Get ID col
  id_col <- id_col %||% config$preprocessing$id_column %||% "RID"
  
  # open literature biomarkers file
  if (!file.exists(lit_biomarker_file_path)) {
    stop("Literature biomarker definition file not found: ", lit_biomarker_file_path)
  }
  loaded_yaml <- yaml::read_yaml(lit_biomarker_file_path)
  literature_defs <- loaded_yaml$biomarkers # Access the list under the 'biomarkers' key
  
  # Validate structure of literature_defs
  for(bname in names(literature_defs)) {
    bdef <- literature_defs[[bname]]
    req_keys <- c("name", "numerator_regions", "denominator_regions", "var_composition")
    has_regions <- all(c("numerator_regions", "denominator_regions", "var_composition") %in% names(bdef))
    has_method <- "calculation_method" %in% names(bdef)
    if (!has_regions && !has_method) stop("Biomarker '", bname, "' definition needs num/den regions OR calculation_method.")
    if (has_method && bdef$calculation_method == "DDS_infcgm" && !"reference_region" %in% names(bdef)) {
      stop("Biomarker '", bname, "' with 'DDS_infcgm' method requires 'reference_region' definition.")
    }
    if(!rlang::is_scalar_character(bdef$name)) stop("Biomarker '", bname, "' needs a valid character name.")
    if(bdef$var_composition != -1) {
      if(!is.list(bdef) || !all(req_keys %in% names(bdef))) stop("Invalid definition for literature biomarker '", bname, "'. Must be list with name, numerator_regions, denominator_regions, var_composition.")
      if(!is.character(bdef$numerator_regions) || length(bdef$numerator_regions)==0) stop("Biomarker '", bname, "' needs non-empty character vector for numerator_regions.")
      if(!is.character(bdef$denominator_regions) || length(bdef$denominator_regions)==0) stop("Biomarker '", bname, "' needs non-empty character vector for denominator_regions.")
      if(!rlang::is_scalar_integerish(bdef$var_composition) || !bdef$var_composition %in% 0:3) stop("Biomarker '", bname, "' needs var_composition between 0 and 3.")
    }

  }

  # --- Prepare Common Elements ---
  all_power_params <- config$power_params
  eq_all_str <- config$model_equations$eq_all_string
  eq_group_str <- config$model_equations$eq_group_string
  if(is.null(eq_all_str) || is.null(eq_group_str)) stop("Missing formula strings.")
  eq_all <- try(stats::as.formula(eq_all_str)); if(inherits(eq_all, "try-error")) stop("Invalid eq_all_str")
  eq_group <- try(stats::as.formula(eq_group_str)); if(inherits(eq_group, "try-error")) stop("Invalid eq_group_str")
  lmer_control <- lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                                    check.conv.singular = "ignore",
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))

  # --- Initialize Results List ---
  all_evaluation_results <- list()
  result_counter <- 0
  
  # --- Loop 1: Through Datasets ---
  message("\n--- Starting Literature Biomarker Evaluation ---")
  for (eval_dset_name in datasets_to_evaluate) {
    if(verbose) message(sprintf("\nProcessing Dataset: %s", eval_dset_name))
    current_eval_data_list <- prepared_data_list[[eval_dset_name]]
    
    # Prepare data for this dataset once
    if (!"data" %in% names(current_eval_data_list) || !is.data.frame(current_eval_data_list$data) ||
        !"data_suv_bi" %in% names(current_eval_data_list) || !is.data.frame(current_eval_data_list$data_suv_bi)) {
      warning(sprintf("Skipping dataset '%s': Missing 'data' or 'data_suv_bi'.", eval_dset_name), call. = FALSE); next
    }
    
    # Add composite reference region
    comp_regs <- c("brainstem", "hemiwm", "whole_cerebellum")
    if(all((comp_regs) %in% names(current_eval_data_list$data_suv_bi))) {
      # Mean
      current_eval_data_list$data_suv_bi$composite <- rowMeans(current_eval_data_list$data_suv_bi[,comp_regs], na.rm = F)
      current_eval_data_list$data_suv_bi$composite <- rowSums(current_eval_data_list$data_suv_bi[,comp_regs], na.rm = F)
    } else {
      warning(sprintf("Can not compute the composite reference region: '%s' missing", comp_regs[!(comp_regs %in% names(current_eval_data_list$data_suv_bi))]))
    }
    
    
    data_clinical <- current_eval_data_list$data
    data_suv <- current_eval_data_list$data_suv_bi
    available_features <- setdiff(names(data_suv), id_col)
    required_clinical_cols <- c(id_col, "time", "DX", "AB")
    if (length(available_features) == 0 || !all(required_clinical_cols %in% names(data_clinical))) {
      warning(sprintf("Skipping dataset '%s': Missing features or required clinical columns.", eval_dset_name), call. = FALSE); next
    }
    

    
    # --- Loop 2: Through Groups ---
    for (eval_group in groups_to_evaluate) {
      if(verbose) message(sprintf("  Evaluating Group: %s", eval_group))
      target_group_dx_val <- if(eval_group == "CI") 1L else 0L
      
      # --- Loop 3: Through Literature Biomarkers ---
      for (biomarker_name in names(literature_defs)) {
        biomarker_def <- literature_defs[[biomarker_name]]
        # if(verbose) message(sprintf("    Biomarker: %s", biomarker_name))
        
        calculated_value <- NULL # Initialize
        eval_status <- "Pending Calculation"
        
        # --- Calculate Biomarker Value based on definition type ---
        if (!is.null(biomarker_def$calculation_method) &&
            biomarker_def$calculation_method == "DDS_infcgm") {
          # --- Braak Staging Calculation ---
          if(verbose) message(sprintf("      Calculating %s via Braak Staging...", biomarker_name))
          # Check if all stage regions exist in this dataset
          all_stage_regions <- unique(unlist(biomarker_def$stage_definitions))
          ref_reg <- biomarker_def$reference_region
          missing_regs <- setdiff(c(all_stage_regions, ref_reg), available_features)
          if(length(missing_regs) > 0) {
            warning(sprintf("Biomarker '%s' on %s: Missing regions required for staging: %s", biomarker_name, eval_dset_name, paste(missing_regs, collapse=",")), call.=FALSE)
            eval_status <- "Missing Staging Regions"
          } else {
            calculated_value <- tryCatch({
              .calculate_DDS(
                dataset_data_list = current_eval_data_list, # Pass list with data & data_suv
                stage_definitions = biomarker_def$stage_definitions,
                reference_region = ref_reg,
                id_col = id_col,
                verbose = FALSE # Keep internal call quiet
              )
              
            }, error = function(e) { NULL })
            if(is.null(calculated_value)) eval_status <- "Staging Calculation Error"
          }
          
        } else if (!is.null(biomarker_def$numerator_regions)) {
          # --- Standard Num/Den Calculation ---
          if(verbose) message(sprintf("      Calculating %s via Num/Den Ratio...", biomarker_name))
          num_regions <- biomarker_def$numerator_regions
          den_regions <- biomarker_def$denominator_regions
          var_comp <- biomarker_def$var_composition
          
          missing_num <- setdiff(num_regions, available_features)
          missing_den <- setdiff(den_regions, available_features)
          if (length(missing_num) > 0 || length(missing_den) > 0) {
            eval_status <- "Missing Regions"
            warning(sprintf("Biomarker '%s' on %s: Missing regions. Num:[%s], Den:[%s]", biomarker_name, eval_dset_name, paste(missing_num, collapse=","), paste(missing_den, collapse=",")), call.=FALSE)
          } else {
            calculated_value <- tryCatch({
              .tr2suvr( # Use the helper for this
                dataset_cohort_data = current_eval_data_list, # Pass list with data & data_suv
                var_composition = var_comp,
                fixed_numerator_regs = num_regions,
                fixed_denominator_regs = den_regions,
                verbose = F
                # Pass UV/Vol if needed
              )
            }, error = function(e) { NULL })
            if(is.null(calculated_value)) eval_status <- "Value Calculation Error"
          }
        } else {
          # Invalid definition
          warning(sprintf("Biomarker '%s' has invalid definition.", biomarker_name), call.=FALSE)
          eval_status <- "Invalid Definition"
        }
        
        # --- Check calculated value and proceed with evaluation ---
        metrics <- NULL # Initialize metrics for this biomarker/dataset/group
        if (eval_status == "Pending Calculation") { # Only proceed if calculation seemed ok
          if (all(is.na(calculated_value))) {
            eval_status <- "Value Calculation NA"
          } else if (any(calculated_value <= 0, na.rm = TRUE)) {
            eval_status <- "Value Calculation NonPositive"
          } else {
            # --- Prepare data and evaluate metrics ---
            eval_data <- data_clinical
            if (nrow(eval_data) != length(calculated_value)) {
              warning(sprintf("Biomarker '%s' on %s/%s: Row mismatch value/clinical.", biomarker_name, eval_dset_name, eval_group), call. = FALSE)
              eval_status <- "Data Merge Failed"
              if (calculate_ci) metrics <- stats::setNames(rep(NA_real_, 9), c("Rep", "Rep.lower", "Rep.upper", "SepAB", "SepAB.lower", "SepAB.upper", "SSE", "SSE.lower", "SSE.upper")) else stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
            } else {
              eval_data$value <- log(calculated_value) # take log to fit LMER
              essential_eval_cols <- c("value", "time", id_col, "DX", "AB")
              eval_data <- stats::na.omit(eval_data[, essential_eval_cols, drop = FALSE])
              
              # Check sufficiency
              min_rows_threshold <- 10; min_group_members_threshold <- 5
              group_subset <- eval_data |> dplyr::filter(.data$DX == target_group_dx_val, .data$AB == TRUE)
              n_group_rows <- nrow(group_subset); n_group_ids <- dplyr::n_distinct(group_subset[[id_col]])
              
              if(n_group_rows < min_rows_threshold || n_group_ids < min_group_members_threshold ) {
                warning(sprintf("Biomarker '%s' on %s/%s: Insufficient data rows (%d) or unique IDs (%d).", biomarker_name, eval_dset_name, eval_group, n_group_rows, n_group_ids), call.=FALSE)
                eval_status <- "Insufficient Data"
                if (calculate_ci) metrics <- stats::setNames(rep(NA_real_, 9), c("Rep", "Rep.lower", "Rep.upper", "SepAB", "SepAB.lower", "SepAB.upper", "SSE", "SSE.lower", "SSE.upper")) else stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
              } else {
                # Call evaluation function
                eval_func <- if(calculate_ci) .feval_group_95CI else .feval_group
                eval_args <- list(data = eval_data, group = eval_group, eq_all = eq_all, eq_group = eq_group, all_power_params = all_power_params, lmer_control = lmer_control)
                if(calculate_ci) eval_args$nsim <- nsim
                metrics_result <- try(do.call(eval_func, eval_args), silent = TRUE)
                
                if (inherits(metrics_result, "try-error")) {
                  warning(sprintf("Biomarker '%s' on %s/%s: Eval func failed: %s", biomarker_name, eval_dset_name, eval_group, conditionMessage(attr(metrics_result,"condition"))), call.=FALSE)
                  eval_status <- "Metrics Calculation Error"
                  if (calculate_ci) metrics <- stats::setNames(rep(NA_real_, 9), c("Rep", "Rep.lower", "Rep.upper", "SepAB", "SepAB.lower", "SepAB.upper", "SSE", "SSE.lower", "SSE.upper")) else stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
                } else if (anyNA(metrics_result)) {
                  warning(sprintf("Biomarker '%s' on %s/%s: Eval func returned NA.", biomarker_name, eval_dset_name, eval_group), call.=FALSE)
                  eval_status <- "Metrics Calculation NA"
                  metrics <- metrics_result
                } else {
                  eval_status <- "Success"
                  metrics <- metrics_result
                }
              } # End sufficiency check else
            } # End data merge check else
          }
        } # End if status was Pending Calculation
        
        # --- Ensure metrics structure exists even if failed ---
        if (is.null(metrics)) {
          if (calculate_ci) metrics <- stats::setNames(rep(NA_real_, 9), c("Rep", "Rep.lower", "Rep.upper", "SepAB", "SepAB.lower", "SepAB.upper", "SSE", "SSE.lower", "SSE.upper"))
          else metrics <- stats::setNames(rep(NA_real_, 3), c("Rep", "SepAB", "SSE"))
        }
        metrics_df <- as.data.frame(as.list(metrics))
        
        # --- Store Result Row ---
        result_counter <- result_counter + 1
        all_evaluation_results[[result_counter]] <- data.frame(
          biomarker_name = biomarker_def$name,
          evaluation_dataset = eval_dset_name,
          evaluation_group = eval_group,
          evaluation_status = eval_status,
          var_composition = biomarker_def$var_composition %||% NA_integer_, # Use NA if method based
          num_regions = paste(biomarker_def$numerator_regions %||% "-", collapse=", "),
          den_regions = paste(biomarker_def$denominator_regions %||% "-", collapse=", "),
          calculation_method = biomarker_def$calculation_method %||% "num_den_ratio",
          source_citation = biomarker_def$source_citation %||% NA_character_,
          stringsAsFactors = FALSE
        ) |> dplyr::bind_cols(metrics_df)
        
      } # End loop over biomarkers
    } # End loop over groups
  } # End loop over datasets
  
  # --- Combine All Results ---
  if (length(all_evaluation_results) == 0) {
    warning("No literature evaluation results were generated.", call. = FALSE)
    final_results_df <- data.frame()
  } else {
    final_results_df <- dplyr::bind_rows(all_evaluation_results)
    message(sprintf("\n--- Literature Evaluation Complete: Generated %d result rows. ---", nrow(final_results_df)))
  }
  
  # --- Save Final Table (Optional) ---
  if (!is.null(output_csv_path)) {
    # ... (logic to create dir and write CSV using readr::write_csv) ...
    tryCatch({
      output_dir <- dirname(output_csv_path)
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      if(dir.exists(output_dir)){
        readr::write_csv(final_results_df, output_csv_path, na = "NA")
        message("Literature evaluation results saved to: ", output_csv_path)
      } else { warning("Failed to create output directory '", output_dir, "'. Cannot write CSV.", call.=FALSE) }
    }, error = function(e) { warning(sprintf("Failed to write literature evaluation results to '%s': %s", output_csv_path, e$message), call. = FALSE) })
  }
  
  invisible(final_results_df)
}