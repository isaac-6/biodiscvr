# --- Internal Helper Functions ---
# These are not exported for the user but used by fevalMethod_...
# No roxygen required unless you want detailed internal docs (@keywords internal)

# Repeatability as error% (Standard Deviation of Residuals)
# Input: A fitted model object (e.g., from lmer)
# Output: Standard deviation of residuals
.fRep <- function(fit) {
  return(stats::sd(stats::residuals(fit)))
}

# Group separation (t-statistic of the time:DX interaction)
# Input: A fitted model object (e.g., from lmer) with a 'time:DX' interaction
# Output: t-value for the 'time:DX' term, or NA if not found/error
.fSep <- function(fit) {
  res <- NA_real_
  tryCatch({
    model_summary <- summary(fit)
    coef_table <- stats::coef(model_summary)
    if ("time:DX" %in% rownames(coef_table)) {
      res <- coef_table["time:DX", "t value"]
    } else {
      warning("'.fSep': Term 'time:DX' not found in model summary.", call. = FALSE)
    }
  }, error = function(e) {
    warning("'.fSep': Error extracting 'time:DX' t-value: ", e$message, call. = FALSE)
  })
  return(res)
}

# Group separation, amyloid (t-statistic of the time:AB interaction)
# Input: A fitted model object (e.g., from lmer) with a 'time:AB' interaction
# Output: t-value for the 'time:AB' term, or NA if not found/error
.fSepAB <- function(fit) {
  res <- NA_real_
  tryCatch({
    model_summary <- summary(fit)
    coef_table <- stats::coef(model_summary)
    if ("time:AB" %in% rownames(coef_table)) {
      res <- coef_table["time:AB", "t value"]
    } else {
      warning("'.fSepAB': Term 'time:AB' not found in model summary.", call. = FALSE)
    }
  }, error = function(e) {
    warning("'.fSepAB': Error extracting 'time:AB' t-value: ", e$message, call. = FALSE)
  })
  return(res)
}

#' Internal helper to calculate metrics for a specific group (CU or CI)
#'
#' Calculates metrics using group-specific parameters retrieved from config.
#'
#' @param data Data frame containing prepared modeling data...
#' @param group Character string, either "CU" or "CI".
#' @param eq_all Formula object for the separation model.
#' @param eq_group Formula object for the within-group model.
#' @param all_power_params List. The *entire* power_params structure loaded
#'   from the config file (expected to have keys like 'CI', 'CU').
#' @param lmer_control Control object for lmer fitting.
#'
#' @return A named numeric vector with 'Rep', 'SepAB', 'SSE'.
#'
#' @noRd
.feval_group <- function(data, group, eq_all, eq_group, all_power_params, lmer_control) {
  
  # --- Input Validation ---
  stopifnot(
    is.data.frame(data),
    all(c("DX", "AB", "value", "time", "RID") %in% names(data)), # Ensure RID is present if used in formulas
    is.character(group), length(group) == 1, group %in% c("CU", "CI"), # Extend if more groups added
    inherits(eq_all, "formula"),
    inherits(eq_group, "formula"),
    is.list(all_power_params), # Check it's a list
    methods::is(lmer_control, "lmerControl")
  )
  
  # --- Get Group-Specific Power Parameters ---
  if (!group %in% names(all_power_params)) {
    stop(sprintf("Power parameters for group '%s' not found in the 'power_params' section of the configuration.", group))
  }
  group_power_params <- all_power_params[[group]]
  
  # --- Validate the retrieved group parameters ---
  required_keys <- c("pct_change", "trial_time_points", "power")
  if (!is.list(group_power_params) || !all(required_keys %in% names(group_power_params))) {
    stop(sprintf("Power parameters for group '%s' in config must be a list containing: %s",
                 group, paste(required_keys, collapse=", ")))
  }
  # Optionally add type checks for pct_change, power (numeric) and trial_time_points (numeric vector)
  
  # --- Set Group Specific DX value ---
  group_DX <- NA_integer_
  if (group == "CU") {
    group_DX <- 0L
  } else if (group == "CI") {
    group_DX <- 1L
  } # Add more groups if needed
  
  # --- Separation (Amyloid effect within the specified DX group) ---
  
  eSepAB <- NA_real_
  data_filtered_sep <- data[data$DX == group_DX & !is.na(data$AB), ]
  
  if (sum(data_filtered_sep$AB == 0) > 30) {
    model_sep <- try(lme4::lmer(eq_all, data = data_filtered_sep, REML = TRUE, control = lmer_control), silent = TRUE)
    if (!inherits(model_sep, "try-error")) {
      eSepAB <- .fSepAB(model_sep) # Assumes eq_all includes time:AB
      if (is.na(eSepAB)) {
        warning(sprintf("Group '%s': Could not extract 'time:AB' t-value from separation model.", group), call. = FALSE)
      }
    } else {
      warning(sprintf("Group '%s': Failed to fit separation model (eq_all): %s", group, conditionMessage(attr(model_sep, "condition"))), call. = FALSE)
    }
  } else {
    eSepAB <- 1
    # warning(sprintf("Group '%s': No non-NA AB data available for separation model. Setting default Sep = 1", group), call. = FALSE)
  }
  
  # --- Repeatability and SSE (within the specified DX group, Amyloid Positive only) ---
  eRep <- NA_real_
  eSSE <- NA_real_
  data_filtered_sse <- data[data$DX == group_DX & data$AB == TRUE & !is.na(data$AB), ]
  
  if (nrow(data_filtered_sse) > 0) {
    model_sse <- try(lme4::lmer(eq_group, data = data_filtered_sse, REML = TRUE, control = lmer_control), silent = TRUE)
    if (!inherits(model_sse, "try-error")) {
      eRep <- .fRep(model_sse) * 100
      
      # Power calculation using extracted parameters
      pow_res <- try(longpower::lmmpower(model_sse,
                                         pct.change = group_power_params$pct_change,
                                         t = unlist(group_power_params$trial_time_points), # Use the vector directly
                                         power = group_power_params$power), silent = TRUE)
      
      if (!inherits(pow_res, "try-error") && !is.null(pow_res$n)) {
        eSSE <- pow_res$n[1]
      } else {
        pow_err_msg <- if(inherits(pow_res, "try-error")) conditionMessage(attr(pow_res, "condition")) else "returned NULL or no 'n' element"
        warning(sprintf("Group '%s': longpower::lmmpower failed or returned unexpected result: %s", group, pow_err_msg), call. = FALSE)
      }
    } else {
      warning(sprintf("Group '%s': Failed to fit SSE model (eq_group): %s", group, conditionMessage(attr(model_sse, "condition"))), call. = FALSE)
    }
  } else {
    warning(sprintf("Group '%s': No AB positive data available for SSE model.", group), call. = FALSE)
  }
  
  # Return named vector
  return(stats::setNames(c(eRep, eSepAB, eSSE), c("Rep", "SepAB", "SSE")))
}







# Place this in R/model_evaluation.R along with .feval_group and others

#' @importFrom lme4 lmer lmerControl bootMer
#' @importFrom longpower lmmpower
#' @importFrom stats sd residuals coef quantile setNames na.omit
#' @importFrom methods is
#'
#' @param data Data frame containing prepared modeling data...
#' @param group Character string, either "CU" or "CI".
#' @param eq_all Formula object for the separation model.
#' @param eq_group Formula object for the within-group model.
#' @param all_power_params List. The *entire* power_params structure...
#' @param lmer_control Control object for lmer fitting.
#' @param nsim Integer. Number of bootstrap simulations for CI calculation.
#'
#' @return A named numeric vector with point estimates and CI bounds:
#'   'Rep', 'Rep.lower', 'Rep.upper', 'SepAB', 'SepAB.lower', 'SepAB.upper',
#'   'SSE', 'SSE.lower', 'SSE.upper'. Values are NA if models/bootstrapping fail.
#'
#' @noRd
.feval_group_95CI <- function(data, group, eq_all, eq_group, all_power_params, lmer_control, nsim) {
  
  # --- Input Validation ---
  stopifnot(
    is.data.frame(data),
    all(c("DX", "AB", "value", "time", "RID") %in% names(data)),
    is.character(group), length(group) == 1, group %in% c("CU", "CI"),
    inherits(eq_all, "formula"),
    inherits(eq_group, "formula"),
    is.list(all_power_params),
    methods::is(lmer_control, "lmerControl"),
    is.numeric(nsim), length(nsim) == 1, nsim > 100 # Ensure reasonable nsim
  )
  nsim <- as.integer(nsim) # Ensure it's integer
  
  # --- Get Group-Specific Power Parameters ---
  # (Identical to .feval_group)
  if (!group %in% names(all_power_params)) stop(sprintf("Power parameters for group '%s' not found.", group))
  group_power_params <- all_power_params[[group]]
  required_keys <- c("pct_change", "trial_time_points", "power")
  if (!is.list(group_power_params) || !all(required_keys %in% names(group_power_params))) stop(sprintf("Power parameters for group '%s' invalid.", group))
  
  # --- Set Group Specific DX value ---
  # (Identical to .feval_group)
  group_DX <- NA_integer_
  if (group == "CU") group_DX <- 0L else if (group == "CI") group_DX <- 1L
  
  # --- Initialize ALL results with NAs ---
  eRep = eRep.lower = eRep.upper = NA_real_
  eSepAB = eSepAB.lower = eSepAB.upper = NA_real_
  eSSE = eSSE.lower = eSSE.upper = NA_real_
  
  # --- Helper for safe quantile calculation ---
  safe_quantile <- function(x, probs) {
    res <- rep(NA_real_, length(probs))
    valid_x <- stats::na.omit(x) # Remove NAs from bootstrap results
    if(length(valid_x) > 99) { # Need at least 100 non-NA points for quantile
      res_try <- try(stats::quantile(valid_x, probs = probs, na.rm = FALSE), silent=TRUE) # na.rm=F as we already omitted
      if(!inherits(res_try, "try-error")) res <- res_try
    } else {
      warning("Not enough valid bootstrap replicates to calculate quantiles.", call. = FALSE)
    }
    return(res)
  }
  
  # --- Separation (Amyloid effect within the specified DX group) ---
  data_filtered_sep <- data[data$DX == group_DX & !is.na(data$AB), ]
  model_sep <- NULL # Initialize model variable
  
  if (nrow(data_filtered_sep) > 0) {
    model_sep <- try(lme4::lmer(eq_all, data = data_filtered_sep, REML = TRUE, control = lmer_control), silent = TRUE)
    
    if (!(sum(data_filtered_sep$AB == 0) > 30)) {
      eSepAB <- NA
      eSepAB.lower <- NA
      eSepAB.upper <- NA
    } else {
      if (!inherits(model_sep, "try-error")) {
        eSepAB <- .fSepAB(model_sep) # Calculate point estimate
        if (is.na(eSepAB)) {
          warning(sprintf("Group '%s': Could not extract point estimate for 'time:AB' t-value.", group), call. = FALSE)
        }
        
        # --- Bootstrap for Separation CI ---
        # Consider parallel options if needed: parallel = "multicore" or "snow"
        b_sep <- try(lme4::bootMer(model_sep, FUN = .fSepAB, nsim = nsim, type = "parametric", use.u = FALSE, re.form=NA), silent = TRUE)
        # use.u=FALSE, re.form=NA are often recommended for parametric bootstrap stability
        
        if (!inherits(b_sep, "try-error")) {
          ci_sep <- safe_quantile(b_sep$t, probs = c(0.025, 0.975))
          eSepAB.lower <- ci_sep[1]
          eSepAB.upper <- ci_sep[2]
        } else {
          warning(sprintf("Group '%s': bootMer failed for separation (SepAB). CI calculation skipped. Error: %s",
                          group, conditionMessage(attr(b_sep, "condition"))), call.=FALSE)
        }
        
      } else {
        warning(sprintf("Group '%s': Failed to fit separation model (eq_all), cannot calculate SepAB or CIs. Error: %s",
                        group, conditionMessage(attr(model_sep, "condition"))), call. = FALSE)
      }
    }
    
    
  } else {
    warning(sprintf("Group '%s': No non-NA AB data for separation model.", group), call. = FALSE)
  }
  
  
  # --- Repeatability and SSE (within the specified DX group, Amyloid Positive only) ---
  data_filtered_sse <- data[data$DX == group_DX & data$AB == TRUE & !is.na(data$AB), ]
  model_sse <- NULL # Initialize model variable
  
  if (nrow(data_filtered_sse) > 0) {
    model_sse <- try(lme4::lmer(eq_group, data = data_filtered_sse, REML = TRUE, control = lmer_control), silent = TRUE)
    
    if (!inherits(model_sse, "try-error")) {
      eRep <- .fRep(model_sse) * 100 # Calculate point estimate
      
      # --- Power calculation for SSE (including potential CIs from longpower) ---
      pow_res <- try(longpower::lmmpower(model_sse,
                                         pct.change = group_power_params$pct_change,
                                         t = unlist(group_power_params$trial_time_points),
                                         power = group_power_params$power), # Request CI from lmmpower
                     silent = TRUE)
      
      if (!inherits(pow_res, "try-error")) {
        if (!is.null(pow_res$n)) eSSE <- pow_res$n[1]
        # Extract CI if available from lmmpower results
        if (!is.null(pow_res$n.CI)) {
          eSSE.lower <- pow_res$n.CI[1] # Assuming [lower, lower, upper, upper]
          eSSE.upper <- pow_res$n.CI[3]
        } else {
          warning(sprintf("Group '%s': longpower::lmmpower did not return confidence intervals for SSE (n.CI is NULL).", group), call.=FALSE)
        }
      } else {
        pow_err_msg <- conditionMessage(attr(pow_res, "condition"))
        warning(sprintf("Group '%s': longpower::lmmpower failed: %s", group, pow_err_msg), call. = FALSE)
      }
      
      # --- Bootstrap for Repeatability CI ---
      b_rep <- try(lme4::bootMer(model_sse, FUN = .fRep, nsim = nsim, type = "parametric", use.u = FALSE, re.form=NA), silent = TRUE)
      
      if (!inherits(b_rep, "try-error")) {
        ci_rep <- safe_quantile(b_rep$t, probs = c(0.025, 0.975)) * 100 # Multiply by 100
        eRep.lower <- ci_rep[1]
        eRep.upper <- ci_rep[2]
      } else {
        warning(sprintf("Group '%s': bootMer failed for repeatability (Rep). CI calculation skipped. Error: %s",
                        group, conditionMessage(attr(b_rep, "condition"))), call.=FALSE)
      }
      
    } else {
      warning(sprintf("Group '%s': Failed to fit SSE model (eq_group), cannot calculate Rep, SSE or CIs. Error: %s",
                      group, conditionMessage(attr(model_sse, "condition"))), call. = FALSE)
    }
  } else {
    warning(sprintf("Group '%s': No AB positive data for SSE model.", group), call. = FALSE)
  }
  
  # --- Return named vector with all results ---
  return(stats::setNames(
    c(eRep, eRep.lower, eRep.upper,
      eSepAB, eSepAB.lower, eSepAB.upper,
      eSSE, eSSE.lower, eSSE.upper),
    c("Rep", "Rep.lower", "Rep.upper",
      "SepAB", "SepAB.lower", "SepAB.upper",
      "SSE", "SSE.lower", "SSE.upper")
  ))
}
