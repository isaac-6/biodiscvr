# Various useful functions


#' Number of CVR combinations given a number of regions to choose from
#' It considers A/B equivalent to B/A. 
#' @param x number of regions to choose from
#' @noRd
.combos <- function(x) {
  (3^x-3*2^x+3)/6
}


# `%||%` <- function(lhs, rhs) {
#   if (!is.null(lhs)) {
#     lhs
#   } else {
#     rhs
#   }
# }


#' Append a data frame row to a CSV file
#'
#' Handles header writing correctly only if the file doesn't exist.
#' Uses readr::write_csv for robust writing.
#'
#' @param results_df A single-row data frame containing the results to append.
#' @param csv_path Path to the CSV file.
#'
#' @return Invisibly returns TRUE on success, FALSE on failure.
#' @noRd
#' @importFrom readr write_csv
.append_to_csv <- function(results_df, csv_path) {
  # Basic validation of inputs
  stopifnot(
    is.data.frame(results_df),
    nrow(results_df) == 1,
    is.character(csv_path),
    length(csv_path) == 1,
    nzchar(csv_path) # Ensure path is not empty
  )
  
  # Check if file exists *before* trying to write
  # This is safer than relying solely on append=TRUE for header logic with readr > 2.0
  file_already_exists <- file.exists(csv_path)
  
  tryCatch({
    readr::write_csv(
      x = results_df,
      file = csv_path,
      append = file_already_exists, # Only append if file exists
      col_names = !file_already_exists # Only write headers if file is new
    )
    invisible(TRUE) # Signal success
  }, error = function(e) {
    # Provide a warning if writing fails
    warning(sprintf("Failed to write/append results to CSV '%s': %s",
                    csv_path, conditionMessage(e)), call. = FALSE)
    invisible(FALSE) # Signal failure
  })
}



# --- Functions to aid the definition of literature biomarkers ---

# --- Helper function `tr2suvr` ---
# It needs to return a numeric vector of SUVR values of the same length as nrow(data_suv).
#' @noRd
.tr2suvr <- function(dataset_cohort_data, # folder with data.csv and data_suv_bi.csv
                     var_composition,
                     fixed_numerator_regs = NULL,
                     fixed_denominator_regs = NULL,
                     verbose = FALSE) {
  
  features <- names(dataset_cohort_data$data_suv_bi)
  chromosome <- rep(1.5,length(features))
  chromosome[features %in% fixed_numerator_regs] <- 0.1
  chromosome[features %in% fixed_denominator_regs] <- 2.9
  
  bio <- .calculate_cvr(chromosome,
                        dataset_cohort_data, # folder with data.csv and data_suv_bi.csv
                        features,
                        var_composition,
                        fixed_numerator_regs = NULL,
                        fixed_denominator_regs = NULL,
                        verbose = FALSE)
  return(bio) 
}

#' @importFrom mclust Mclust mclustBIC
#' @importFrom stats dnorm uniroot sd median quantile IQR na.omit aggregate as.formula
#' @importFrom dplyr filter select all_of left_join distinct n_distinct between

# --- Helper function for GMM intersection ---
# (Ensure this helper is defined before .calculate_DDS)
#' @noRd
.get_intersectionGMM <- function(gmm_fit) {
  # Check if fit is valid and has required components
  if (is.null(gmm_fit) || !is.list(gmm_fit$parameters) ||
      !all(c("mean", "variance") %in% names(gmm_fit$parameters)) ||
      length(gmm_fit$parameters$mean) != 2 ||
      length(gmm_fit$parameters$variance$sigma) != 2) {
    warning("Invalid gmm_fit object passed to .get_intersectionGMM.", call.=FALSE)
    return(NA_real_)
  }
  
  sd1 <- sqrt(gmm_fit$parameters$variance$sigma[1])
  sd2 <- sqrt(gmm_fit$parameters$variance$sigma[2])
  # Handle list vs vector for mean safely
  mu1 <- gmm_fit$parameters$mean[[1]][1] # Take first element if it's somehow a list/vector itself
  mu2 <- gmm_fit$parameters$mean[[2]][1]
  
  # Check for valid means/sds
  if(anyNA(c(sd1, sd2, mu1, mu2)) || sd1 <= 0 || sd2 <= 0) {
    warning("Invalid means or standard deviations in GMM fit.", call.=FALSE)
    return(NA_real_)
  }
  
  diff_gaussians <- function(x) {
    dnorm(x, mean = mu1, sd = sd1) - dnorm(x, mean = mu2, sd = sd2)
  }
  
  # Use tryCatch for uniroot as it can fail
  intersection_root <- NA_real_
  tryCatch({
    # Check if signs at means are different (suggests intersection between them)
    d1 <- diff_gaussians(mu1)
    d2 <- diff_gaussians(mu2)
    
    if (sign(d1) != sign(d2)) {
      # Find root between the means
      root_result <- uniroot(diff_gaussians, interval = sort(c(mu1, mu2)), extendInt = "no", check.conv=TRUE)
      intersection_root <- root_result$root
    } else {
      # Fallback: Find root closest to mu1 on the side towards mu2
      # Search interval needs careful definition - extend cautiously
      search_interval <- sort(c(mu1, mu2 + 3 * sd2)) # Extend 3 SDs from mu2
      if (sign(diff_gaussians(search_interval[1])) != sign(diff_gaussians(search_interval[2]))) {
        root_result <- uniroot(diff_gaussians, interval = search_interval, extendInt = "upX", check.conv=TRUE)
        intersection_root <- root_result$root
      } else {
        # Further fallback or warning if no root found in extended interval
        warning("Could not find GMM intersection reliably between/near means.", call.=FALSE)
        # Maybe return midpoint or NA? Let's return NA.
        intersection_root <- NA_real_
      }
    }
  }, error = function(e){
    warning("Error finding GMM intersection: ", conditionMessage(e), call.=FALSE)
    intersection_root <- NA_real_
  })
  
  return(intersection_root)
}


#' Calculate Biomarker Value based on Braak-like Staging and GMM Thresholds
#'
#' Internal function to compute a specific SUVR based on data-driven stages.
#'
#' @param dataset_data_list List containing `$data` and `$data_suv_bi`.
#' @param stage_definitions List of character vectors defining regions for stages (e.g., `$st1`, `$st2`).
#' @param reference_region Character vector of reference region column name(s).
#' @param id_col Character string for the unique subject identifier.
#' @param verbose Logical. Print progress messages?
#'
#' @return Numeric vector of calculated log-transformed biomarker values, or NULL on failure.
#' @noRd
.calculate_DDS <- function(dataset_data_list,
                                           stage_definitions,
                                           reference_region,
                                           id_col,
                                           verbose = FALSE) {
  
  if(verbose) message("    Calculating DDS...")
  data_clinical <- dataset_data_list$data
  data_suv <- dataset_data_list$data_suv_bi
  
  # --- Validate Inputs ---
  stopifnot(
    is.list(stage_definitions), length(stage_definitions) > 0,
    is.character(reference_region), length(reference_region) > 0,
    is.character(id_col), length(id_col) == 1
  )
  # Check reference region exists
  if(!all(reference_region %in% names(data_suv))) {
    warning("Reference region(s) not found in data_suv_bi: ", paste(setdiff(reference_region, names(data_suv)), collapse=", "), call.=FALSE); return(NULL)
  }
  # Check ID exists
  if(!id_col %in% names(data_clinical)) {
    warning("ID column '", id_col, "' not found in data", call.=FALSE); return(NULL)
  }
  
  # --- Calculate SUVR and Thresholds per Stage ---
  stage_suvr_list <- list()
  tau_posit_stages <- numeric(length(stage_definitions))
  names(tau_posit_stages) <- names(stage_definitions)
  
  for (stage_name in names(stage_definitions)) {
    target_regions <- stage_definitions[[stage_name]]
    if(verbose) message(sprintf("      Processing stage %s (%d regions)...", stage_name, length(target_regions)))
    
    # Check target regions exist
    missing_target <- setdiff(target_regions, names(data_suv))
    if(length(missing_target) > 0) {
      warning(sprintf("Stage '%s': Target regions not found: %s. Cannot calculate SUVR/threshold.",
                      stage_name, paste(missing_target, collapse=", ")), call.=FALSE)
      stage_suvr_list[[stage_name]] <- rep(NA_real_, nrow(data_suv)) # Store NA vector
      tau_posit_stages[stage_name] <- NA_real_
      next # Skip to next stage
    }
    
    # Calculate SUVR for this stage using the provided function
    # Assuming tr2suvr handles composition internally based on region lists
    # And returns a single vector of SUVRs
    stage_suvr <- tryCatch({
      .tr2suvr(dataset_cohort_data = dataset_data_list, # folder with data.csv and data_suv_bi.csv
               var_composition = 0,
               fixed_numerator_regs = target_regions,
               fixed_denominator_regs = reference_region,
               verbose = FALSE)
    }, error = function(e) {
      warning(sprintf("Stage '%s': tr2suvr failed. Error: %s", stage_name, conditionMessage(e)), call.=FALSE)
      return(NULL)
    })
    
    if(is.null(stage_suvr) || length(stage_suvr) != nrow(data_suv)) {
      warning(sprintf("Stage '%s': SUVR calculation failed or returned incorrect length.", stage_name), call.=FALSE)
      stage_suvr_list[[stage_name]] <- rep(NA_real_, nrow(data_suv))
      tau_posit_stages[stage_name] <- NA_real_
      next
    }
    stage_suvr_list[[stage_name]] <- stage_suvr # Store the calculated SUVR vector
    
    # Fit GMM on SUVR (handle potential errors/warnings)
    # Use only non-NA, finite, positive values for GMM fitting
    valid_exp_suvr <- stage_suvr
    valid_exp_suvr <- valid_exp_suvr[!is.na(valid_exp_suvr) & is.finite(valid_exp_suvr) & valid_exp_suvr > 0]
    
    if(length(valid_exp_suvr) < 10) { # Need sufficient points for GMM
      warning(sprintf("Stage '%s': Not enough valid data points (%d) to fit GMM.", stage_name, length(valid_exp_suvr)), call.=FALSE)
      tau_posit_stages[stage_name] <- NA_real_
      next
    }
    
    # Suppress messages/prints from Mclust if desired
    gmm_fit <- suppressMessages(try(mclust::Mclust(valid_exp_suvr, G = 2, verbose = FALSE), silent = TRUE))
    
    if (inherits(gmm_fit, "try-error") || is.null(gmm_fit$parameters)) {
      warning(sprintf("Stage '%s': GMM fitting failed. Error: %s", stage_name, ifelse(inherits(gmm_fit, "try-error"), conditionMessage(attr(gmm_fit, "condition")), "Fit invalid")), call.=FALSE)
      tau_posit_stages[stage_name] <- NA_real_
    } else {
      # Calculate intersection threshold
      tau_posit_stages[stage_name] <- .get_intersectionGMM(gmm_fit)
      if(is.na(tau_posit_stages[stage_name])) warning(sprintf("Stage '%s': Failed to find GMM intersection threshold.", stage_name), call.=FALSE)
      else if(verbose) message(sprintf("      Stage %s threshold = %.4f", stage_name, tau_posit_stages[stage_name]))
    }
  } # End stage loop
  
  # Check if any thresholds were successfully calculated
  if(all(is.na(tau_posit_stages))) {
    warning("Failed to calculate GMM thresholds for all stages.", call.=FALSE); return(NULL)
  }
  
  # --- Determine Highest Stage Reached ---
  df_stage_suvr <- as.data.frame(stage_suvr_list) # Combine stage SUVRs
  # Ensure column order matches stage order for indexing later
  df_stage_suvr <- df_stage_suvr[, names(stage_definitions), drop=FALSE]

  # Create positivity matrix (handle NA thresholds by treating comparison as FALSE)
  valid_thresholds <- tau_posit_stages[!is.na(tau_posit_stages)]
  valid_stage_names <- names(valid_thresholds)
  if(length(valid_stage_names) == 0) { warning("No valid thresholds available for positivity matrix.", call.=FALSE); return(NULL)}
  
  n_rows <- nrow(df_stage_suvr)
  positivity_matrix <- matrix(FALSE, nrow = n_rows, ncol = length(valid_stage_names),
                              dimnames = list(NULL, valid_stage_names))
  for(stage_name in valid_stage_names) {
    positivity_matrix[, stage_name] <- df_stage_suvr[[stage_name]] > valid_thresholds[stage_name]
  }
  # Handle NAs in positivity matrix if SUVRs were NA - treat as FALSE
  positivity_matrix[is.na(positivity_matrix)] <- FALSE
  
  # Get the index (1-based) of the rightmost positive stage
  rightmost_positivity_idx <- apply(positivity_matrix, 1, function(row) {
    pos_indices <- which(row)
    ifelse(length(pos_indices) > 0, max(pos_indices), 1) # Return 1 if none are positive
  })
  
  # --- Propagate Baseline Stage ---
  # Need to merge rightmost_positivity_idx with RID and time first
  stage_data <- data.frame(
    placeholder_id = data_clinical[[id_col]],
    time = data_clinical$time, # Assuming 'time' exists
    stage_idx = rightmost_positivity_idx
  )
  names(stage_data)[1] <- id_col # Use correct ID name
  
  # Find baseline stage for each ID (minimum time)
  stage_data <- stage_data[order(stage_data[[id_col]], stage_data$time), ]
  baseline_stage <- stage_data[!duplicated(stage_data[[id_col]]), c(id_col, "stage_idx")]
  names(baseline_stage)[2] <- "baseline_stage_idx"
  
  # Merge baseline stage back and propagate
  stage_data <- dplyr::left_join(stage_data, baseline_stage, by = id_col)
  # Original code implies using the baseline stage always. Let's follow that.
  # If baseline_stage_idx is NA (e.g., only NA times?), use original index? Default to 0.
  final_stage_idx <- stage_data$baseline_stage_idx
  final_stage_idx[is.na(final_stage_idx)] <- 0 # Default to 0 if baseline couldn't be determined
  
  # --- Select Value from Corresponding Stage SUVR ---
  final_value_vector <- rep(NA_real_, n_rows)
  valid_final_indices <- final_stage_idx > 0 & final_stage_idx <= ncol(df_stage_suvr)
  
  # Use matrix indexing for efficiency
  row_indices <- seq_len(n_rows)[valid_final_indices]
  col_indices <- final_stage_idx[valid_final_indices]
  # Create indexing matrix
  idx_matrix <- cbind(row_indices, col_indices)
  # Extract values
  final_value_vector[valid_final_indices] <- as.matrix(df_stage_suvr)[idx_matrix]
  
  # Handle cases where final_stage_idx was 0 (no positivity) -> assign NA or a baseline value?
  # Let's assign NA if no stage was positive at baseline.
  final_value_vector[final_stage_idx == 0] <- NA_real_
  
  if(verbose) message("    Finished calculating Braak-staging biomarker values.")
  return(final_value_vector) # Return the log-transformed value
}


