#' Calculate Composite Value Ratio (Internal Helper)
#'
#' Calculates CVR = (Target Composite / Reference Composite) for each row of data,
#' based on feature selection/weighting defined by a chromosome.
#' Used internally, typically as part of a GA fitness evaluation.
#'
#' @param chromosome Numeric vector. Encodes feature selection and weighting.
#'   Values < 1 select target features (weight = 1 - value).
#'   Values between 1 and 2 are ignored (?).
#'   Values > 2 select reference features (weight = value - 2).
#' @param dataset_cohort_data List. Should contain at least `$data_suv_bi`,
#'   a data frame with SUV values for different regions (columns).
#'   Optionally contains `$data_uv` and `$data_vol` if `var_composition = 0`.
#' @param features Character vector. Names of all possible features (columns in
#'   `data_suv_bi`) corresponding to the indices of `chromosome`.
#' @param fixed_numerator_regs fixed numerator regions.
#' @param fixed_denominator_regs fixed denominator regions.
#' @param var_composition Integer. Method for combining features:
#'   \itemize{
#'     \item `0`: Volume/UV weighted average (Requires `$data_uv`, `$data_vol`).
#'     \item `1`: Simple arithmetic mean (`rowMeans`).
#'     \item `2`: Weighted arithmetic mean (weights from chromosome).
#'     \item `3`: Geometric mean (`exp(mean(log(values)))`).
#'   }
#' @param verbose Logical. Print verbose messages (mainly for debugging)?
#'
#' @return A numeric vector of the calculated `(Target/Reference)` values,
#'   with the same length as `nrow(dataset_cohort_data$data_suv_bi)`.
#'   Returns a vector of `NA` values if critical errors occur (e.g., missing
#'   columns, invalid chromosome leading to no target/reference regions, division
#'   by zero). Returns `-Inf` only if specifically calculated (e.g., log(0)).
#'
#' @keywords internal
.calculate_cvr <- function(chromosome,
                           dataset_cohort_data,
                           features,
                           var_composition,
                           fixed_numerator_regs = NULL,
                           fixed_denominator_regs = NULL,
                           verbose = FALSE) {
  
  # random test
  # chromosome <- mapply(function(min, max) runif(1, min, max), min_bounds, max_bounds)
  
  # browser()
  
  # --- Input Validation ---
  stopifnot(
    is.numeric(chromosome),
    is.list(dataset_cohort_data),
    !is.null(dataset_cohort_data$data_suv_bi),
    is.character(features),
    length(chromosome) == length(features),
    var_composition %in% 0:3
  )
  
  d_suv <- dataset_cohort_data$data_suv_bi
  n_rows <- nrow(d_suv)
  # Create default NA vector for early return on error
  na_result <- rep(NA_real_, n_rows)
  
  # --- Feature Selection ---
  # Helper for verbose messaging
  message_verbose <- function(...) if(verbose) message(...)
  
  # 1. Determine default selection based on chromosome
  target_indices <- which(chromosome < 1)
  reference_indices <- which(chromosome > 2)
  target_features_chr <- features[target_indices]
  reference_features_chr <- features[reference_indices]
  
  # 2. Set initial values (will be overridden if fixed args are provided)
  target_features <- target_features_chr
  reference_features <- reference_features_chr
  
  # 3. Override numerator if fixed_numerator_regs is provided
  if (!is.null(fixed_numerator_regs)) {
    message_verbose(".calculate_cvr: Overriding target features with fixed_numerator_regs.")
    target_features <- fixed_numerator_regs
  }
  
  # 4. Override denominator if fixed_denominator_regs is provided
  if (!is.null(fixed_denominator_regs)) {
    message_verbose(".calculate_cvr: Overriding reference features with fixed_denominator_regs.")
    reference_features <- fixed_denominator_regs
  }
  
  # Remove any potential NAs from feature names themselves
  target_features <- target_features[!is.na(target_features)]
  reference_features <- reference_features[!is.na(reference_features)]
  
  n_target <- length(target_features)
  n_reference <- length(reference_features)
  
  # --- Check if valid target/reference features selected ---
  if (n_target == 0 || n_reference == 0) {
    if(verbose) message(".calculate_cvr: Chromosome resulted in no target (", n_target, ") or no reference (", n_reference, ") features.")
    return(na_result) # Invalid chromosome for CVR calculation
  }
  
  # --- Check if selected features exist in the data ---
  available_cols <- names(d_suv)
  missing_target_cols <- setdiff(target_features, available_cols)
  missing_reference_cols <- setdiff(reference_features, available_cols)
  
  if (length(missing_target_cols) > 0) {
    warning(".calculate_cvr: Required target feature(s) not found in data: ", paste(missing_target_cols, collapse=", "), call. = FALSE)
    return(na_result)
  }
  if (length(missing_reference_cols) > 0) {
    warning(".calculate_cvr: Required reference feature(s) not found in data: ", paste(missing_reference_cols, collapse=", "), call. = FALSE)
    return(na_result)
  }
  
  # --- Define Helper Function for Composite Calculation ---
  calculate_composite <- function(selected_features, indices, type) {
    len <- length(selected_features)
    if (len == 0) return(NA_real_) # Should be caught above, but safeguard
    
    # --- Single Feature Case ---
    if (len == 1) {
      comp_val <- d_suv[[selected_features]] # Use [[ ]] for single column vector
      return(comp_val)
    }
    
    # --- Multiple Feature Cases ---
    # Subset data safely now we know columns exist
    sub_data <- d_suv[, selected_features, drop = FALSE]
    
    if (var_composition == 0) {
      # --- Volume/UV Weighted Average ---
      # Requires d_uv and d_vol which might not be standard
      d_uv <- dataset_cohort_data$data_uv
      d_vol <- dataset_cohort_data$data_vol
      if (is.null(d_uv) || is.null(d_vol)) {
        warning("var_composition=0 requires 'data_uv' and 'data_vol' in dataset list, which are missing. Cannot calculate composite.", call. = FALSE)
        return(rep(NA_real_, n_rows))
      }
      # Check columns exist in d_uv/d_vol too
      if (!all(selected_features %in% names(d_uv)) || !all(selected_features %in% names(d_vol))) {
        warning("var_composition=0: Selected features missing from 'data_uv' or 'data_vol'.", call.=FALSE)
        return(rep(NA_real_, n_rows))
      }
      # Calculation (handle potential zero sum in denominator)
      sum_uv <- rowSums(d_uv[, selected_features, drop = FALSE], na.rm = TRUE) # Decide NA handling
      sum_vol <- rowSums(d_vol[, selected_features, drop = FALSE], na.rm = TRUE)
      comp_val <- ifelse(sum_vol == 0, NA_real_, sum_uv / sum_vol)
      
    } else if (var_composition == 1) {
      # --- Simple Arithmetic Mean ---
      comp_val <- rowMeans(sub_data, na.rm = TRUE) # Decide NA handling
      
    } else if (var_composition == 2) {
      # --- Weighted Arithmetic Mean ---
      if (type == "target") {
        weights <- 1 - chromosome[indices]
      } else { # reference
        weights <- chromosome[indices] - 2
      }
      sum_weights <- sum(weights)
      if (abs(sum_weights) < .Machine$double.eps) { # Check for near-zero sum
        warning("Sum of weights for weighted average is zero. Cannot calculate composite.", call.=FALSE)
        return(rep(NA_real_, n_rows))
      }
      # Use matrix multiplication for efficiency
      # Handle NAs in sub_data appropriately before multiplication? Assume non-NA for now or impute earlier.
      weighted_vals <- as.matrix(sub_data) %*% weights
      comp_val <- weighted_vals / sum_weights
      
    } else if (var_composition == 3) {
      # --- Geometric Mean ---
      # Log transform, take mean, exponentiate. Handles zeros/negatives better.
      log_sub_data <- log(as.matrix(sub_data))
      # Check for -Inf resulting from log(0) or NaN from log(negative)
      if (any(!is.finite(log_sub_data))) {
        warning("Non-positive values encountered in data for geometric mean calculation. Result may contain NA/-Inf.", call.=FALSE)
        # Replace non-finite with NA before rowMeans to avoid propagating Inf/NaN?
        log_sub_data[!is.finite(log_sub_data)] <- NA
      }
      mean_log_vals <- rowMeans(log_sub_data, na.rm = TRUE) # Decide NA handling
      comp_val <- exp(mean_log_vals)
    }
    return(comp_val)
  } # End of calculate_composite helper
  
  # --- Calculate Target and Reference Composites ---
  target_composite <- calculate_composite(target_features, target_indices, "target")
  reference_composite <- calculate_composite(reference_features, reference_indices, "reference")
  
  # --- Check if composite calculations failed ---
  if (anyNA(target_composite) || anyNA(reference_composite)) {
    if (verbose) message(".calculate_cvr: NA values produced during composite calculation.")
    # Returning NA result vector propagates the failure
    return(na_result)
  }
  
  # --- Calculate Ratio ---
  # Avoid division by zero or negative reference values
  ratio <- ifelse(reference_composite <= 0, NA_real_, target_composite / reference_composite)
  
  # # Take log, handle non-positive ratios resulting from negative target
  # log_ratio <- ifelse(ratio <= 0, NA_real_, log(ratio))
  # 
  # # Check final result for issues
  # if (anyNA(log_ratio)){
  #   warning(".calculate_cvr: NA values generated during final ratio/log calculation (division by zero, log of non-positive?).", call.=FALSE)
  # } else if (any(!is.finite(log_ratio))) {
  #   warning(".calculate_cvr: Non-finite values (-Inf) generated during final log calculation (likely target=0).", call.=FALSE)
  #   # Keep -Inf as it might be meaningful? Or convert to NA? Let's keep for now.
  #   # log_ratio[!is.finite(log_ratio)] <- NA_real_
  # }
  
  return(ratio)
}



# --- Internal Helper: Fitness Function for a Single Dataset ---
# This is the function the GA maximizes.
# It takes a 'candidate' solution (e.g., numeric vector of coefficients),
# the prepared data for *one* dataset, model formulas, and power parameters.
# It must return a SINGLE numeric fitness score. Higher is better.
# chromosome is the vector to be optimised (maximised)
# features are the column names of the selected variables (regions) to be included in the search
# var_composition determines how regions are merged ( 0 = volume-weighed, 1 = arithmetic mean, 2 = weighed mean, 3 = geometric mean)
# cohort data is the list with the cohort files (data, data_suv_bi)
# group is either "CU" or "CI", for cognitively unimpaired or impaired.
# eq_all is the linear mixed-effects model equation when including AB+ and AB-, so we can calculate the separation
# eq_group is the linear mixed-effects model equation when only considering a single group (AB+) and either CU or CI
.calculate_fitness <- function(chromosome, 
                               features, 
                               var_composition, 
                               fixed_numerator_regs = NULL,
                               fixed_denominator_regs = NULL,
                               cohort_data, 
                               group, 
                               eq_all, 
                               eq_group, 
                               all_power_params, 
                               lmer_control) {
  
  # get CVR
  ratio <- .calculate_cvr(chromosome,
                          cohort_data,
                          features,
                          var_composition,
                          fixed_numerator_regs,
                          fixed_denominator_regs,
                          verbose = FALSE)
  
  # Take log, handle non-positive ratios and NA
  if (any(is.na(ratio)) || any(ratio <= 0)) {
    return(-Inf)  # Return -Inf, the worst fitness
  } else {
    log_ratio <- log(ratio)  # Take the log element-wise
  }
  
  
  
  # --- 1. Apply candidate solution to SUV data ---
  # This is highly dependent on what 'candidate_solution' represents.
  # Example: If it's weights for columns in data_suv_bi
  suv_data <- cohort_data$data_suv_bi
  main_data <- cohort_data$data
  
  main_data$value <- log_ratio

  
  # --- 2. Get evaluation metrics ---
  # We don't need CIs for every fitness evaluation
  # Note: data expects 'value', 'time', 'ID', 'DX', 'AB' columns
  # Ensure 'time' is appropriately defined/scaled in main_data
  metrics <- try(.feval_group(data = main_data,
                              group,
                             eq_all,
                             eq_group,
                             all_power_params, 
                             lmer_control),
                 silent = TRUE)
  
  if (inherits(metrics, "try-error") || anyNA(metrics)) {
    warning(".feval_group failed or returned NA during fitness calculation. Returning poor fitness.", call. = FALSE)
    return(-Inf)
  }
  
  # --- 3. Combine metrics into a single fitness score ---
  # truncate separation
  
  if (metrics["SepAB"][[1]] > 2.6) {
    tsep = 2.6
  } else {
    tsep = metrics["SepAB"][[1]]
  }
  fitness_score <- tsep / metrics["SSE"][[1]]
  # Handle potential NA values from .feval_group if they weren't caught
  if(is.na(fitness_score)) fitness_score <- -Inf
  
  return(fitness_score)
}






#' Combine Multiple Fitness Scores
#'
#' Calculates an overall fitness score by combining a vector of individual
#' fitness scores. The combination involves the product of the scores and
#' their cosine similarity to a reference vector.
#'
#' @param fitnesses A non-empty numeric vector of individual fitness scores.
#' @param reference_vector A numeric vector of the same length as `fitnesses`
#'   representing the reference or "ideal" profile. If NULL (default), a
#'   vector of all ones is used.
#' @param na.rm Logical. Should missing values (NA) be removed before
#'   calculation? Defaults to `FALSE`. If `FALSE` and NAs are present, the
#'   result will be `NA`.
#'
#' @return A single numeric value representing the combined fitness score,
#'   calculated as `prod(fitnesses) * (cosine_similarity)^2`. Returns `NA`
#'   if `na.rm` is `FALSE` and `NA` values are present in inputs, or if
#'   cosine similarity cannot be computed (e.g., due to zero vectors when
#'   `na.rm=TRUE` removes all elements).
#'
#' @details
#' The function first calculates the base fitness as the product of all
#' non-NA elements in `fitnesses` (`prod(fitnesses, na.rm = na.rm)`).
#' It then calculates the cosine similarity between the `fitnesses` vector
#' and the `reference_vector`. Cosine similarity measures the cosine of the
#' angle between two vectors and ranges from -1 (exactly opposite) to 1
#' (exactly the same direction), with 0 indicating orthogonality.
#' The formula used is `cos(theta) = dot_product(A, B) / (norm(A) * norm(B))`.
#' The function handles potential division by zero if either vector has a
#' zero norm after NA removal (in which case cosine similarity is treated as 0).
#' The final combined fitness is `base_fitness * (cosine_similarity)^2`.
#' Squaring the cosine similarity makes the contribution always positive and
#' emphasizes alignment (values closer to 1 or -1 contribute more positively
#' than values closer to 0).
#'
#' Note: If `fitnesses` contains zeros, the `base_fitness` (and thus the final
#' result) will be zero. Consider if this is the desired behavior.
#' Also, if `fitnesses` can contain negative values, the sign of the
#' `base_fitness` will depend on the number of negative values.
#'
#' @export
#' @examples
#' fitness_scores <- c(0.8, 0.9, 0.7)
#' reference <- c(1, 1, 1) # Aim for high scores
#' .multi_fitness(fitness_scores, reference)
#' .multi_fitness(fitness_scores) # Uses default reference of all 1s
#'
#' fitness_scores_neg <- c(0.8, -0.9, 0.7)
#' .multi_fitness(fitness_scores_neg) # Note the impact of negative values on prod()
#'
#' fitness_scores_aligned <- c(2, 4, 6)
#' reference_aligned <- c(1, 2, 3) # Perfectly aligned direction
#' .multi_fitness(fitness_scores_aligned, reference_aligned) # Cosine squared will be 1
#'
#' fitness_scores_orthog <- c(1, 0)
#' reference_orthog <- c(0, 1) # Orthogonal vectors
#' .multi_fitness(fitness_scores_orthog, reference_orthog) # Cosine will be 0
#'
#' .multi_fitness(c(0, 0), c(1, 1)) # Base is 0
#' .multi_fitness(c(1, 1), c(0, 0)) # Cosine is 0 (due to zero norm)
#'
#' .multi_fitness(c(1, NA, 2)) # Returns NA
#' .multi_fitness(c(1, NA, 2), na.rm = TRUE) # Calculates based on c(1, 2)
#'
.multi_fitness <- function(fitnesses, reference_vector = NULL, na.rm = FALSE) {
  
  # --- Input Validation ---
  if (!is.numeric(fitnesses) || length(fitnesses) == 0) {
    stop("'fitnesses' must be a non-empty numeric vector.")
  }
  if (!is.null(reference_vector)) {
    if (!is.numeric(reference_vector)) {
      stop("'reference_vector', if provided, must be numeric.")
    }
    if (length(fitnesses) != length(reference_vector)) {
      stop("'fitnesses' and 'reference_vector' must have the same length.")
    }
  }
  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("'na.rm' must be TRUE or FALSE.")
  }
  
  # --- Handle NA ---
  if (na.rm) {
    if (!is.null(reference_vector)) {
      complete_cases <- stats::complete.cases(fitnesses, reference_vector)
      fit_clean <- fitnesses[complete_cases]
      ref_clean <- reference_vector[complete_cases]
    } else {
      fit_clean <- stats::na.omit(fitnesses)
      # ref_clean will be generated later based on length of fit_clean
    }
    # Check if any data remains after NA removal
    if (length(fit_clean) == 0) {
      warning("No complete cases remaining after NA removal.")
      return(NA_real_)
    }
  } else {
    # Check for NAs if na.rm is FALSE
    if (anyNA(fitnesses) || (!is.null(reference_vector) && anyNA(reference_vector))) {
      return(NA_real_)
    }
    fit_clean <- fitnesses
    ref_clean <- reference_vector # Will be set below if NULL
  }
  
  # --- Set Default Reference Vector ---
  if (is.null(reference_vector)) {
    # If original was NULL, generate based on length *after* NA removal
    ref_clean <- rep(1.0, length(fit_clean))
  }
  # If reference_vector was provided, ref_clean was already set during NA handling
  
  # --- Calculate Base Magnitude (Product of Absolute Values) ---
  # Check for zeros before calculating prod(abs())
  if (any(abs(fit_clean) < .Machine$double.eps)) {
    base_magnitude <- 0.0
  } else {
    base_magnitude <- prod(abs(fit_clean))
  }
  
  # --- Calculate Cosine Similarity ---
  # Use cleaned vectors 'fit_clean' and 'ref_clean'
  dot_product <- sum(fit_clean * ref_clean) # sum handles na.rm internally
  
  norm_fit <- sqrt(sum(fit_clean^2))
  norm_ref <- sqrt(sum(ref_clean^2))
  
  # Handle zero norms to avoid NaN/Inf
  if (norm_fit < .Machine$double.eps || norm_ref < .Machine$double.eps) {
    # If either vector has zero magnitude, cosine is undefined or 0
    cosine_sim <- 0.0
    if (norm_fit < .Machine$double.eps && norm_ref < .Machine$double.eps) {
      # Handle the 0/0 case - arguably similarity is 1 if both are zero vectors?
      # Or more commonly treated as 0 or NA. Let's use 0.
      cosine_sim <- 0.0 # Or perhaps 1.0 or NA_real_ depending on desired meaning
    }
  } else {
    cosine_sim <- dot_product / (norm_fit * norm_ref)
    # Clamp to [-1, 1] due to potential floating point inaccuracies
    cosine_sim <- max(-1.0, min(1.0, cosine_sim))
  }
  
  
  # --- Calculate Combined Fitness Magnitude ---
  combined_magnitude <- base_magnitude * (cosine_sim^2)
  
  # --- Determine Final Sign ---
  # Check if *any* non-zero element in the cleaned fitness vector was negative
  # Handle the case where the result is exactly zero separately
  final_sign <- 1.0
  if (combined_magnitude > .Machine$double.eps) { # Only assign sign if non-zero
    if (any(fit_clean < -(.Machine$double.eps))) { # Check for any negative values
      final_sign <- -1.0
    }
  }
  
  # --- Apply Sign ---
  final_combined_fitness <- final_sign * combined_magnitude
  
  return(final_combined_fitness)
}




#' Calculate Reference SSE Direction Vector
#'
#' Calculates the estimated sample size (SSE) for a specific reference
#' biomarker across multiple cohorts. The resulting vector is typically used as
#' the `direction_reference` argument in `calculate_multi_cohort_fitness`.
#'
#' @param reference_dna Numeric vector. The 'dna' representing the reference
#'   biomarker (e.g., a standard composite or a previously known good solution).
#' @param all_data List. The complete data structure returned by
#'   `check_and_prepare_data()`.
#' @param datasets_to_use Character vector. Names of the datasets within
#'   `all_data` to include (must match those used in the main fitness function).
#' @param get_suvr_function Function. The function to calculate the biomarker
#'   value from 'dna' and a dataset's data frame. Signature: `function(dna, data_df)`.
#' @param eq_group Formula. The formula for the within-group SSE models.
#' @param power_params List. Parameters for `longpower::lmmpower` containing:
#'   `pct_change`, `t`, `power`.
#' @param sse_group Integer. Defines the group for SSE calculation: 0 for CI
#'   (DX=1, A+), 1 for CU (DX=0, A+).
#' @param verbose Logical. If TRUE, print progress messages. Defaults to TRUE.
#'
#' @return A named numeric vector where names are the dataset names and values
#'   are the calculated reference SSEs. Returns `NA` for datasets where the
#'   calculation fails (e.g., insufficient data, model convergence issues).
#'   Returns `NULL` if any essential input is invalid.
#'
#' @details
#' This function iterates through the specified datasets, calculates the SUVR
#' for the *single* `reference_dna` provided, fits the `lmer` model based on
#' `eq_group` for the specified `sse_group`, and runs `lmmpower` to estimate
#' the required sample size (SSE). This provides a benchmark SSE vector representing
#' performance proportionality for the reference biomarker.
#'
#' @export
#' @importFrom lme4 lmer lmerControl
#' @importFrom longpower lmmpower
#' @importFrom stats update na.omit
calculate_reference_direction <- function(reference_dna,
                                          all_data,
                                          datasets_to_use,
                                          get_suvr_function,
                                          eq_group,
                                          power_params,
                                          sse_group = 0,
                                          verbose = TRUE) {
  
  # --- Basic Input Validation ---
  if (!is.numeric(reference_dna) || !is.list(all_data) ||
      !is.character(datasets_to_use) || !all(datasets_to_use %in% names(all_data)) ||
      !is.function(get_suvr_function) || !inherits(eq_group, "formula") ||
      !is.list(power_params) || !all(c("pct_change", "t", "power") %in% names(power_params)) ||
      !sse_group %in% c(0, 1)) {
    stop("Invalid input arguments provided to calculate_reference_direction.")
    return(NULL) # Should not be reached due to stop, but for safety
  }
  if (length(datasets_to_use) == 0) {
    stop("`datasets_to_use` cannot be empty.")
    return(NULL)
  }
  
  
  # Initialize result vector
  reference_sse_vector <- numeric(length(datasets_to_use))
  names(reference_sse_vector) <- datasets_to_use
  reference_sse_vector[] <- NA_real_ # Initialize with NA
  
  # Common model control
  lmer_control <- lme4::lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE,
                                    check.conv.singular = "ignore",
                                    optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
  
  if (verbose) cat("Calculating Reference SSE Direction for datasets:", paste(datasets_to_use, collapse=", "), "\n")
  
  # --- Loop through specified datasets ---
  for (dataset_name in datasets_to_use) {
    if (verbose) cat("  Processing dataset:", dataset_name, "...\n")
    
    # --- 1. Get Data & Calculate Reference SUVR ---
    if (!"data" %in% names(all_data[[dataset_name]]) || is.null(all_data[[dataset_name]]$data)) {
      warning("Dataset '", dataset_name, "' missing 'data' component. Skipping.", call. = FALSE)
      next # Skip to next dataset
    }
    current_data <- all_data[[dataset_name]]$data
    
    suvr_result <- try(get_suvr_function(reference_dna, current_data), silent = TRUE)
    if (inherits(suvr_result, "try-error") || length(suvr_result) != nrow(current_data)) {
      warning("get_suvr_function failed or returned incorrect length for reference_dna in dataset '", dataset_name, "'. Skipping.", call. = FALSE)
      next
    }
    current_data$value <- suvr_result # Add the reference biomarker value
    
    # --- 2. Calculate SSE ---
    sse_data <- NULL
    if (sse_group == 0) { # CI A+
      sse_data <- current_data[which(current_data$DX == 1 & current_data$AB == TRUE & !is.na(current_data$AB)), ]
    } else if (sse_group == 1) { # CU A+
      sse_data <- current_data[which(current_data$DX == 0 & current_data$AB == TRUE & !is.na(current_data$AB)), ]
    }
    
    if (!is.null(sse_data) && nrow(sse_data) > 5) { # Need sufficient data
      sse_model <- try(lme4::lmer(eq_group, data = sse_data, REML = TRUE, control = lmer_control), silent = TRUE)
      
      if (!inherits(sse_model, "try-error")) {
        pow_sse <- try(longpower::lmmpower(sse_model, # Correct package name
                                           pct.change = power_params$pct_change,
                                           t = power_params$t,
                                           power = power_params$power), silent = TRUE)
        
        if (!inherits(pow_sse, "try-error") && !is.null(pow_sse$n)) {
          ref_sse <- pow_sse$n[1]
          if(is.finite(ref_sse) && ref_sse > 0){
            reference_sse_vector[dataset_name] <- ref_sse # Store valid SSE
            if(verbose) cat("    Reference SSE:", round(ref_sse), "\n")
          } else {
            warning("lmmpower returned non-positive or non-finite SSE for reference_dna in '", dataset_name, "'.", call. = FALSE)
          }
        } else {
          warning("lmmpower failed for reference_dna SSE calculation in '", dataset_name, "'.", call. = FALSE)
        }
      } else {
        warning("lmer model failed for reference_dna SSE calculation in '", dataset_name, "'.", call. = FALSE)
      }
    } else {
      warning("Insufficient data (n=", nrow(sse_data), ") for reference_dna SSE calculation group ", sse_group, " in '", dataset_name, "'.", call.=FALSE)
    }
    
  } # End loop
  
  # Final check for any NAs
  if (any(is.na(reference_sse_vector))) {
    warning("Could not calculate reference SSE for all specified datasets: ",
            paste(names(reference_sse_vector)[is.na(reference_sse_vector)], collapse=", "),
            ". Returning vector with NAs.", call.=FALSE)
  }
  
  if(verbose) cat("Finished calculating reference direction.\n")
  return(reference_sse_vector)
}





