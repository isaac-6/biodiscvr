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