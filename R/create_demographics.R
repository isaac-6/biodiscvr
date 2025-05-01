#' Create Demographics Summary Table
#'
#' Generates a summary table of key demographic and clinical variables for the
#' total cohort, Cognitively Unimpaired (CU, DX=0), and Cognitively Impaired
#' (CI, DX=1) groups within a processed dataset.
#'
#' @param processed_dataset_list List. A single element from the list returned
#'   by `preprocess_data()`, containing the filtered 'data' data frame and
#'   potentially other processed data frames for ONE dataset (e.g., `preprocessed_data$ADNI`).
#'   The 'data' data frame is expected to contain columns like RID, DX, AGE, sex,
#'   PTEDUCAT, AB, SITE.
#' @param dataset_name Character string. The name of the dataset being summarized
#'   (used if writing output).
#' @param output_csv_path Character string or NULL. If provided, the full path
#'   to save the demographics table as a CSV file. If the file exists, it will
#'   be overwritten. e.g., output/demo.ADNI.csv
#' @param id_col Character string. Name of the unique individual identifier column.
#'   Defaults to "RID".
#'
#' @return A data frame containing the formatted demographic summary table.
#'   Rows represent characteristics (N, sex, Age, Education, AB status, Sites),
#'   and columns represent groups (Total, CU, CI). Returns NULL if the input
#'   'data' dataframe is missing or essential columns are absent.
#'
#' @details
#'   The function calculates:
#'   - N: Number of unique individuals.
#'   - sex: N (%) for each level (assuming codes 1 and 2).
#'   - Age (baseline): Median \\[Q1, Q2\\] (Min-Max). Uses the 'AGE' column.
#'   - Education (years): Median \\[Q1, Q2\\] (Min-Max). Uses the 'PTEDUCAT' column.
#'   - Amyloid Positive (%): Percentage of individuals with AB=1. Uses the 'AB' column.
#'   - N Sites: Count of unique sites. Uses the 'SITE' column.
#'
#'   It handles missing columns gracefully by reporting NA and issuing warnings.
#'
#' @export
#' @importFrom stats median quantile IQR sd na.omit na.pass aggregate
#' @importFrom readr write_csv
create_demographics_table <- function(processed_dataset_list,
                                      dataset_name,
                                      output_csv_path = NULL,
                                      id_col = "RID") {
  
  # --- Input Validation ---
  if (!is.list(processed_dataset_list) || !"data" %in% names(processed_dataset_list) ||
      !is.data.frame(processed_dataset_list$data)) {
    warning(sprintf("Dataset '%s': Input 'processed_dataset_list' is invalid or missing the 'data' data frame. Cannot create demographics table.", dataset_name), call. = FALSE)
    return(NULL)
  }
  data_df <- processed_dataset_list$data
  
  if (nrow(data_df) == 0) {
    warning(sprintf("Dataset '%s': Input 'data' data frame has 0 rows. Cannot create demographics table.", dataset_name), call. = FALSE)
    return(NULL)
  }
  
  # --- Check required columns ---
  required_cols <- c(id_col, "DX", "AGE", "sex", "PTEDUCAT", "AB", "SITE")
  available_cols <- names(data_df)
  missing_cols <- setdiff(required_cols, available_cols)
  if (length(missing_cols) > 0) {
    warning(sprintf("Dataset '%s': Missing required columns for demographics table: %s. Results for these variables will be NA.",
                    dataset_name, paste(missing_cols, collapse=", ")), call.=FALSE)
  }
  
  # --- Internal Helper Functions for Formatting ---
  
  # Format Numeric: Median [Q1, Q3] (Min-Max)
  .format_numeric <- function(x) {
    x_clean <- stats::na.omit(x)
    if (length(x_clean) == 0) return("NA")
    med <- stats::median(x_clean)
    q1q3 <- stats::quantile(x_clean, probs = c(0.25, 0.75))
    minmax <- range(x_clean)
    sprintf("%.1f [%.1f, %.1f] (%.1f-%.1f)", med, q1q3[1], q1q3[2], minmax[1], minmax[2])
  }
  
  # Format Categorical: N (%) - expects a frequency table
  .format_categorical <- function(freq_table, total_n) {
    if(total_n == 0) return("NA")
    paste(sprintf("%d (%.1f%%)", freq_table, (freq_table / total_n) * 100), collapse = " | ")
  }
  
  # Format Percentage Positive (for AB status)
  .format_percent_positive <- function(x) {
    x_clean <- stats::na.omit(x)
    n_total <- length(x_clean)
    if (n_total == 0) return("NA")
    n_pos <- sum(x_clean == 1) # Assuming 1 means positive
    sprintf("%.1f%% (%d/%d)", (n_pos / n_total) * 100, n_pos, n_total)
  }
  
  
  # --- Calculate Summaries for Each Group ---
  groups <- list(
    Total = data_df,
    CU = if("DX" %in% available_cols) data_df[data_df$DX == 0 & !is.na(data_df$DX), ] else data.frame(), # DX=0
    CI = if("DX" %in% available_cols) data_df[data_df$DX == 1 & !is.na(data_df$DX), ] else data.frame()  # DX=1
  )
  
  summary_list <- list()
  
  for (group_label in names(groups)) {
    group_data <- groups[[group_label]]
    results <- list()
    
    # N (Unique Individuals)
    n_unique <- if(id_col %in% names(group_data)) length(unique(group_data[[id_col]])) else NA
    results$"N (Unique IDs)" <- n_unique
    
    # Check if group has data before proceeding
    if (nrow(group_data) == 0) {
      results$"sex (1 | 2)" <- "NA"
      results$"Age (years)" <- "NA"
      results$"Education (years)" <- "NA"
      results$"Amyloid Positive (%)" <- "NA"
      results$"N Sites" <- "NA"
    } else {
      # sex (assuming codes 1 and 2)
      if ("sex" %in% available_cols) {
        sex_table <- table(factor(group_data$sex, levels=c(1, 2))) # Ensure both levels present
        results$"sex (1 | 2)" <- .format_categorical(sex_table, nrow(group_data))
      } else { results$"sex (1 | 2)" <- "NA" }
      
      # Age at Baseline (Use unique per ID if multiple rows per ID exist)
      if ("AGE" %in% available_cols && id_col %in% available_cols) {
        age_per_id <- stats::aggregate(as.formula(paste("AGE ~", id_col)), data = group_data, FUN = function(x) x[1], na.action=na.pass) # Take first age if multiple
        results$"Age (years)" <- .format_numeric(age_per_id$AGE)
      } else { results$"Age (years)" <- "NA" }
      
      # Education (Use unique per ID)
      if ("PTEDUCAT" %in% available_cols && id_col %in% available_cols) {
        educ_per_id <- stats::aggregate(as.formula(paste("PTEDUCAT ~", id_col)), data = group_data, FUN = function(x) x[1], na.action=na.pass)
        results$"Education (years)" <- .format_numeric(educ_per_id$PTEDUCAT)
      } else { results$"Education (years)" <- "NA" }
      
      # Amyloid Positivity (%)
      if ("AB" %in% available_cols && id_col %in% available_cols) {
        ab_per_id <- stats::aggregate(as.formula(paste("AB ~", id_col)), data = group_data, FUN = function(x) x[1], na.action=na.pass)
        results$"Amyloid Positive (%)" <- .format_percent_positive(ab_per_id$AB)
      } else { results$"Amyloid Positive (%)" <- "NA" }
      
      # Number of Sites
      if ("SITE" %in% available_cols) {
        results$"N Sites" <- length(unique(stats::na.omit(group_data$SITE)))
      } else { results$"N Sites" <- "NA" }
    }
    
    summary_list[[group_label]] <- as.data.frame(results)
  }
  
  # Combine into final table
  demographics_table <- do.call(cbind, summary_list)
  # Add characteristics as row names (or first column)
  demographics_table <- cbind(Characteristic = rownames(demographics_table), demographics_table)
  rownames(demographics_table) <- NULL
  
  
  # --- Save to CSV (Optional) ---
  if (!is.null(output_csv_path)) {
    tryCatch({
      output_dir <- dirname(output_csv_path)
      if (!dir.exists(output_dir)) {
        message("Creating output directory for demographics CSV: ", output_dir)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      if(dir.exists(output_dir)){
        readr::write_csv(demographics_table, output_csv_path)
        message("Demographics table saved to: ", output_csv_path)
      } else {
        warning("Failed to create output directory '", output_dir, "'. Cannot write demographics CSV.", call.=FALSE)
      }
    }, error = function(e) {
      warning(sprintf("Failed to write demographics table to '%s': %s", output_csv_path, e$message), call. = FALSE)
    })
  }
  
  # Return the table
  return(demographics_table)
}