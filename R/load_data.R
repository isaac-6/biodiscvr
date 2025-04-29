#' Load Datasets from a Structured Directory
#'
#' Scans a root directory for subdirectories, each representing a dataset.
#' Within each dataset subdirectory, it attempts to load specific CSV files
#' ('data.csv' and 'data_suv_bi.csv').
#'
#' @param root_path Character string. The path to the main directory containing
#'   the dataset subdirectories.
#' @param verbose Logical. If TRUE, print messages about which folders are being
#'   processed and files loaded. Defaults to TRUE.
#'
#' @return A named list where:
#'   - Each name corresponds to a subdirectory name (dataset identifier).
#'   - Each element is itself a named list containing the loaded data frames.
#'     Currently, the possible names within this inner list are 'data' (from
#'     'data.csv') and 'data_suv_bi' (from 'data_suv_bi.csv'). If a file
#'     is not found or fails to load for a given dataset, the corresponding
#'     element will be missing from the inner list for that dataset.
#'   Returns an empty list (`list()`) if `root_path` does not exist, is not a
#'   directory, or contains no subdirectories with the target files.
#'
#' @examples
#' # --- Setup Example Directory Structure ---
#' # (Usually you wouldn't do this inside examples, but it's needed here
#' #  for a self-contained, runnable example)
#' \dontrun{
#'   # Create a temporary root directory
#'   temp_root <- file.path(tempdir(), "my_datasets")
#'   if (dir.exists(temp_root)) unlink(temp_root, recursive = TRUE) # Clean up first
#'   dir.create(temp_root)
#'
#'   # Create Dataset A
#'   dir_a <- file.path(temp_root, "dataset_A")
#'   dir.create(dir_a)
#'   write.csv(data.frame(id = 1:3, value = rnorm(3)),
#'             file.path(dir_a, "data.csv"), row.names = FALSE)
#'   write.csv(data.frame(id = 1:3, score = runif(3)),
#'             file.path(dir_a, "data_suv_bi.csv"), row.names = FALSE)
#'
#'   # Create Dataset B (missing data_suv_bi.csv)
#'   dir_b <- file.path(temp_root, "dataset_B")
#'   dir.create(dir_b)
#'   write.csv(data.frame(id = 10:12, value = rnorm(3)),
#'             file.path(dir_b, "data.csv"), row.names = FALSE)
#'
#'   # --- Load the Datasets ---
#'   all_data <- load_datasets(temp_root)
#'
#'   # --- Access the Data ---
#'   # List loaded datasets
#'   print(names(all_data))
#'
#'   # Access data.csv from dataset_A
#'   print(head(all_data$dataset_A$data))
#'
#'   # Access data_suv_bi.csv from dataset_A
#'   print(head(all_data$dataset_A$data_suv_bi))
#'
#'   # Access data.csv from dataset_B
#'   print(head(all_data$dataset_B$data))
#'
#'   # Attempt to access missing file (will be NULL)
#'   print(all_data$dataset_B$data_suv_bi)
#'
#'   # --- Clean Up ---
#'   unlink(temp_root, recursive = TRUE)
#' }
#'
#' @export
#' @importFrom readr read_csv
#' @importFrom tools file_path_sans_ext
load_datasets <- function(root_path, verbose = TRUE) {
  
  # --- Input Validation ---
  if (!is.character(root_path) || length(root_path) != 1) {
    stop("'root_path' must be a single character string.")
  }
  if (!dir.exists(root_path)) {
    warning("Directory not found: ", root_path)
    return(list())
  }
  
  # --- Identify Potential Dataset Folders ---
  # List only immediate subdirectories, not recursive, return only names
  dataset_folders <- list.dirs(root_path, full.names = FALSE, recursive = FALSE)
  
  if (length(dataset_folders) == 0) {
    warning("No subdirectories found in: ", root_path)
    return(list())
  }
  
  if (verbose) {
    message("Found potential dataset folders: ", paste(dataset_folders, collapse = ", "))
  }
  
  # --- Define Files to Load (could be an argument later) ---
  files_to_load <- c("data.csv", "data_suv_bi.csv")
  # Use file names without extension as keys in the inner list
  file_keys <- tools::file_path_sans_ext(files_to_load)
  
  # --- Initialize Results List ---
  loaded_data <- list()
  
  # --- Loop Through Each Dataset Folder ---
  for (folder_name in dataset_folders) {
    dataset_path <- file.path(root_path, folder_name)
    if (verbose) message("Processing: ", folder_name)
    
    dataset_content <- list() # To store data frames for this specific dataset
    
    # --- Attempt to Load Each Target File ---
    for (i in seq_along(files_to_load)) {
      file_name <- files_to_load[i]
      file_key <- file_keys[i]
      file_path <- file.path(dataset_path, file_name)
      
      if (file.exists(file_path)) {
        if (verbose) message("  Attempting to load: ", file_name)
        tryCatch({
          # Using readr::read_csv for speed and consistency (tibbles)
          # Suppress column type messages unless debugging is needed
          df <- readr::read_csv(file_path, show_col_types = FALSE, progress = FALSE)
          dataset_content[[file_key]] <- df
          if (verbose) message("    -> Success.")
        }, error = function(e) {
          warning("Failed to read '", file_name, "' in '", folder_name, "'. Error: ", e$message, call. = FALSE)
        })
      } else {
        if (verbose) message("  File not found: ", file_name)
        # Optionally, issue a warning here if a file *must* exist
        # warning("Required file '", file_name, "' not found in '", folder_name, "'.", call. = FALSE)
      }
    } # End loop through files_to_load
    
    # --- Add to Main List if any data was loaded for this folder ---
    if (length(dataset_content) > 0) {
      loaded_data[[folder_name]] <- dataset_content
    } else {
      if (verbose) message("  -> No target files found or loaded for this folder.")
    }
    
  } # End loop through dataset_folders
  
  # --- Final Checks ---
  if (length(loaded_data) == 0) {
    warning("No data successfully loaded from any subdirectory in: ", root_path)
  } else {
    if (verbose) message("Finished loading. Found data for ", length(loaded_data), " dataset(s).")
  }
  
  return(loaded_data)
}

