#' Iteratively Refine Biomarker by Single-Region Ablation to Minimize SSE
#'
#' Takes an initial biomarker definition (numerator and denominator regions) and
#' iteratively removes the single region whose removal most decreases the
#' aggregated Sample Size Estimate (SSE) across specified datasets and groups.
#' The process stops when no single region removal further decreases the SSE.
#'
#' @param initial_numerator_regs Character vector of initial numerator regions.
#' @param initial_denominator_regs Character vector of initial denominator regions.
#' @param var_composition Numeric. Parameter indicating how features are combined
#'   (e.g., 0 = volume-weighed, 1 = arithmetic mean). Passed to internal CVR calculation.
#' @param prepared_data_list List. The output from `preprocess_datasets()`,
#'   containing filtered data for multiple datasets.
#' @param config List. The loaded main configuration object from `check_and_prepare_data$config`.
#'   Must contain `preprocessing$id_column`, `model_equations$eq_group_string`,
#'   and `power_params`.
#' @param datasets_to_evaluate Character vector. Names of the datasets within
#'   `prepared_data_list` to use for evaluating SSE during ablation.
#' @param groups_to_evaluate Character vector. Groups ("CU", "CI") to evaluate
#'   within each dataset. Defaults to `c("CU", "CI")`.
#' @param features Character vector or NULL. A specific list of feature names
#'   (columns in `data_suv_bi`) to consider. If NULL (default), common features
#'   across all `datasets_to_evaluate` will be identified and used. The initial
#'   numerator/denominator regions must be a subset of these features.
#' @param verbose Logical. Print progress messages? Defaults to TRUE.
#'
#' @return A list containing:
#'   \item{final_numerator_regs}{Character vector of the refined numerator regions.}
#'   \item{final_denominator_regs}{Character vector of the refined denominator regions.}
#'   \item{initial_aggregated_sse}{Numeric. The aggregated SSE of the initial biomarker.}
#'   \item{final_aggregated_sse}{Numeric. The aggregated SSE of the refined biomarker.}
#'   \item{ablation_log}{A data frame tracking each removed region, the iteration,
#'     and the aggregated SSE before and after its removal.}
#'   \item{per_dataset_sse_initial}{Named list. SSE for each dataset/group with initial biomarker.}
#'   \item{per_dataset_sse_final}{Named list. SSE for each dataset/group with final biomarker.}
#'
#' @details
#' The function works as follows:
#' 1. Identifies common features across `datasets_to_evaluate` if `features` is NULL.
#'    Validates `initial_numerator_regs` and `initial_denominator_regs` against these.
#' 2. Calculates the baseline aggregated SSE for the initial biomarker.
#' 3. Enters an iterative loop:
#'    a. In each iteration, it temporarily removes one region at a time from the
#'       current numerator, recalculates the CVR, and then calculates the aggregated SSE.
#'    b. It does the same for regions in the current denominator.
#'    c. It identifies the single region (from num or den) whose removal results
#'       in the lowest new aggregated SSE.
#'    d. If this lowest new SSE is better (lower) than the current baseline SSE,
#'       the identified region is permanently removed, the baseline SSE is updated,
#'       and the removal is logged. The loop continues.
#'    e. If no single region removal improves the aggregated SSE, the loop stops.
#' The aggregation of SSE across datasets/groups is currently a simple mean.
#'
#' @export
#' @importFrom dplyr bind_rows filter select all_of left_join n_distinct
#' @importFrom rlang `%||%` .data sym
#' @importFrom stats setNames na.omit as.formula
#' @importFrom lme4 lmerControl
#' @importFrom methods is
region_ablation <- function(initial_numerator_regs,
                                         initial_denominator_regs,
                                         var_composition,
                                         prepared_data_list,
                                         config,
                                         datasets_to_evaluate,
                                         groups_to_evaluate = c("CU", "CI"),
                                         features = NULL,
                                         verbose = TRUE) {
  
  # --- Initial Setup & Validation ---
  stopifnot(
    is.character(initial_numerator_regs), length(initial_numerator_regs) >= 1, # Need at least one
    is.character(initial_denominator_regs), length(initial_denominator_regs) >= 1,
    is.numeric(var_composition), var_composition %in% 0:3,
    is.list(prepared_data_list), length(prepared_data_list) > 0,
    is.list(config),
    is.character(datasets_to_evaluate), length(datasets_to_evaluate) > 0,
    is.character(groups_to_evaluate), all(groups_to_evaluate %in% c("CU", "CI")),
    is.null(features) || is.character(features),
    is.logical(verbose)
  )
  missing_dsets <- setdiff(datasets_to_evaluate, names(prepared_data_list))
  if (length(missing_dsets) > 0) stop("Datasets not found: ", paste(missing_dsets, collapse=", "))
  id_col <- config$preprocessing$id_column %||% "RID"
  all_power_params <- config$power_params
  eq_group_str <- config$model_equations$eq_group_string
  if(is.null(eq_group_str)) stop("Missing 'eq_group_string' in config.")
  eq_group <- try(stats::as.formula(eq_group_str)); if(inherits(eq_group, "try-error")) stop("Invalid eq_group_str")
  lmer_control <- lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                                    check.conv.singular = "ignore",
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
  
  # --- Determine Features to Use ---
  if (is.null(features)) {
    if(verbose){message("Identifying common features across evaluation datasets...")}
    common_features_list <- list()
    for (dset_name in datasets_to_evaluate) {
      dset_data_list <- prepared_data_list[[dset_name]]
      if (!"data_suv_bi" %in% names(dset_data_list) || !is.data.frame(dset_data_list$data_suv_bi)) next
      current_dset_features <- setdiff(names(dset_data_list$data_suv_bi), id_col)
      if (length(current_dset_features) > 0) common_features_list[[dset_name]] <- current_dset_features
    }
    if (length(common_features_list) == 0) stop("No datasets with valid SUV features found.")
    features <- Reduce(intersect, common_features_list)
    if (length(features) == 0) stop("No common features found across datasets.")
    if(verbose){message(sprintf("Using %d common features for ablation.", length(features)))}
  }
  # Validate initial regions against the final feature set
  if (!all(initial_numerator_regs %in% features)) stop("Initial numerator regions not in provided/common features.")
  if (!all(initial_denominator_regs %in% features)) stop("Initial denominator regions not in provided/common features.")
  
  
  # --- Helper to evaluate SSE for a given set of regions ---
  # Returns a named list of SSEs: $Dataset_Group -> SSE_value
  .evaluate_sse_for_regions <- function(num_regs, den_regs) {
    if(length(num_regs) == 0 || length(den_regs) == 0) return(stats::setNames(Inf, "AggregatedSSE_Error")) # Penalize empty
    
    sse_results_list <- list()
    for (dset_name in datasets_to_evaluate) {
      current_eval_data_list <- prepared_data_list[[dset_name]]
      # Reconstruct chromosome for .calculate_cvr based on *all* available features in this dataset
      dset_available_features <- setdiff(names(current_eval_data_list$data_suv_bi), id_col)
      # Ensure the current num/den are subsets of what's available in this dataset
      if(!all(num_regs %in% dset_available_features) || !all(den_regs %in% dset_available_features)) {
        warning(sprintf("Skipping SSE eval for %s: current num/den regions not all available.", dset_name), call.=FALSE)
        next
      }
      
      chromosome <- rep(1.5, length(dset_available_features))
      chromosome[dset_available_features %in% num_regs] <- 0.5
      chromosome[dset_available_features %in% den_regs] <- 2.5
      
      cvr_val <- tryCatch(.calculate_cvr(chromosome, current_eval_data_list, dset_available_features, var_composition, verbose=FALSE), error=function(e) NULL)
      if(is.null(cvr_val) || all(is.na(cvr_val))) {
        warning(sprintf("CVR calculation failed for %s during SSE evaluation.", dset_name), call.=FALSE)
        next
      }
      
      eval_data <- current_eval_data_list$data
      if(nrow(eval_data) != length(cvr_val)) {
        warning(sprintf("Row mismatch CVR/clinical for %s.", dset_name), call.=FALSE); next
      }
      eval_data$value <- cvr_val
      essential_cols <- c("value", "time", id_col, "DX", "AB")
      eval_data <- stats::na.omit(eval_data[, essential_cols, drop = FALSE])
      
      for (grp in groups_to_evaluate) {
        sse_val <- .feval_group_sse_only(eval_data, grp, eq_group, all_power_params, lmer_control, id_col)
        sse_results_list[[paste(dset_name, grp, sep="_")]] <- sse_val
      }
    }
    return(sse_results_list)
  }
  
  # --- Initialize ---
  current_num_regs <- initial_numerator_regs
  current_den_regs <- initial_denominator_regs
  ablation_log_list <- list() # Store data frames
  iteration <- 0
  max_ablation_iterations <- length(initial_numerator_regs) + length(initial_denominator_regs) - 1 # Max possible removals
  
  if(verbose) message("--- Starting Iterative Ablation to Minimize Aggregated SSE ---")
  if(verbose) message(sprintf("Initial Num (%d): %s", length(current_num_regs), paste(current_num_regs, collapse=",")))
  if(verbose) message(sprintf("Initial Den (%d): %s", length(current_den_regs), paste(current_den_regs, collapse=",")))
  
  # --- Calculate Initial Baseline SSE ---
  per_dataset_sse_initial <- .evaluate_sse_for_regions(current_num_regs, current_den_regs)
  current_baseline_agg_sse <- mean(unlist(per_dataset_sse_initial), na.rm = TRUE)
  if(is.nan(current_baseline_agg_sse)) current_baseline_agg_sse <- Inf # Handle all NA case
  if(verbose) message(sprintf("Iteration 0: Baseline Aggregated SSE = %.2f", current_baseline_agg_sse))
  
  # --- Iteration Loop ---
  while (iteration < max_ablation_iterations) {
    iteration <- iteration + 1
    if(verbose) message(sprintf("\n--- Ablation Iteration %d (Current SSE: %.2f) ---", iteration, current_baseline_agg_sse))
    best_new_agg_sse_this_iter <- current_baseline_agg_sse # Start with current best
    region_to_remove_candidate <- NULL
    type_of_region_to_remove <- NULL # "Num" or "Den"
    sse_after_best_removal_list <- NULL
    
    # Regions eligible for removal in this iteration
    # (must leave at least one region in num and one in den)
    eligible_num_to_test <- if(length(current_num_regs) > 1) current_num_regs else character(0)
    eligible_den_to_test <- if(length(current_den_regs) > 1) current_den_regs else character(0)
    
    if(length(eligible_num_to_test) == 0 && length(eligible_den_to_test) == 0) {
      if(verbose) message("No more regions eligible for removal (num/den would become empty). Stopping.")
      break
    }
    
    # --- Test removing each numerator region ---
    for (region in eligible_num_to_test) {
      temp_num <- setdiff(current_num_regs, region)
      if(verbose) message(sprintf("  Testing removal of Num:%s...", region))
      ablated_sse_list <- .evaluate_sse_for_regions(temp_num, current_den_regs)
      ablated_agg_sse <- mean(unlist(ablated_sse_list), na.rm = TRUE)
      if(is.nan(ablated_agg_sse)) ablated_agg_sse <- Inf
      
      if(verbose) message(sprintf("    -> Aggregated SSE = %.2f", ablated_agg_sse))
      if (ablated_agg_sse < best_new_agg_sse_this_iter) {
        best_new_agg_sse_this_iter <- ablated_agg_sse
        region_to_remove_candidate <- region
        type_of_region_to_remove <- "Num"
        sse_after_best_removal_list <- ablated_sse_list
      }
    }
    
    # --- Test removing each denominator region ---
    for (region in eligible_den_to_test) {
      temp_den <- setdiff(current_den_regs, region)
      if(verbose) message(sprintf("  Testing removal of Den:%s...", region))
      ablated_sse_list <- .evaluate_sse_for_regions(current_num_regs, temp_den)
      ablated_agg_sse <- mean(unlist(ablated_sse_list), na.rm = TRUE)
      if(is.nan(ablated_agg_sse)) ablated_agg_sse <- Inf
      
      if(verbose) message(sprintf("    -> Aggregated SSE = %.2f", ablated_agg_sse))
      if (ablated_agg_sse < best_new_agg_sse_this_iter) {
        best_new_agg_sse_this_iter <- ablated_agg_sse
        region_to_remove_candidate <- region
        type_of_region_to_remove <- "Den"
        sse_after_best_removal_list <- ablated_sse_list
      }
    }
    
    # --- Decision: Remove region or stop ---
    improvement_threshold <- 1e-5 # Only remove if SSE decreases meaningfully
    # since we are removing regions, we are ok with a small improvement, as it simplified the biomarker
    # if we were aggregating regions, the improvement needed would be higher
    if (!is.null(region_to_remove_candidate) && (current_baseline_agg_sse - best_new_agg_sse_this_iter) > improvement_threshold) {
      if(verbose) message(sprintf("  => Removing %s:%s. New Aggregated SSE = %.2f (was %.2f)",
                                  type_of_region_to_remove, region_to_remove_candidate, best_new_agg_sse_this_iter, current_baseline_agg_sse))
      # Log before permanent removal
      ablation_log_list[[length(ablation_log_list) + 1]] <- data.frame(
        iteration = iteration,
        removed_region = region_to_remove_candidate,
        removed_from = type_of_region_to_remove,
        sse_before_removal = current_baseline_agg_sse,
        sse_after_removal = best_new_agg_sse_this_iter,
        stringsAsFactors = FALSE
      )
      # Permanently remove
      if(type_of_region_to_remove == "Num") current_num_regs <- setdiff(current_num_regs, region_to_remove_candidate)
      else current_den_regs <- setdiff(current_den_regs, region_to_remove_candidate)
      # Update baseline for next iteration
      current_baseline_agg_sse <- best_new_agg_sse_this_iter
      # Update per-dataset SSEs for the new baseline
      per_dataset_sse_initial <- sse_after_best_removal_list # This becomes the new baseline
    } else {
      if(verbose) message("No single region removal improved aggregated SSE significantly. Stopping iteration.")
      break # Exit while loop
    }
  } # End while loop
  
  # --- Prepare Final Output ---
  ablation_log_df <- if(length(ablation_log_list) > 0) dplyr::bind_rows(ablation_log_list) else data.frame()
  per_dataset_sse_final <- per_dataset_sse_initial # The last successful baseline
  
  if(verbose) message("\n--- Iterative Ablation Finished ---")
  if(verbose) message(sprintf("Final Num (%d): %s", length(current_num_regs), paste(current_num_regs, collapse=",")))
  if(verbose) message(sprintf("Final Den (%d): %s", length(current_den_regs), paste(current_den_regs, collapse=",")))
  
  return(list(
    final_numerator_regs = current_num_regs,
    final_denominator_regs = current_den_regs,
    initial_aggregated_sse = mean(unlist(.evaluate_sse_for_regions(initial_numerator_regs, initial_denominator_regs)), na.rm=TRUE), # Re-calc for safety
    final_aggregated_sse = current_baseline_agg_sse,
    ablation_log = ablation_log_df,
    per_dataset_sse_initial = .evaluate_sse_for_regions(initial_numerator_regs, initial_denominator_regs),
    per_dataset_sse_final = per_dataset_sse_final
  ))
}