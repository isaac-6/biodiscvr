library(testthat)
library(biodiscvr)

# Setup paths [cite: 12, 16]
synth_data_root_dir <- system.file("synthdata", package = "biodiscvr", mustWork = TRUE)
path_to_pkg_files <- system.file("files", package = "biodiscvr", mustWork = TRUE)
pkg_synth_config <- "config_synth.yaml"
pkg_synth_dict <- "dict_suv_synth.csv"

test_that("End-to-end synthetic workflow functions correctly", {
  
  # 1. Test Loading [cite: 18]
  loaded_data <- load_datasets(root_path = synth_data_root_dir)
  expect_type(loaded_data, "list")
  expect_true(length(loaded_data) > 0)
  
  # 2. Test Preprocessing [cite: 18, 19]
  preprocessed_output <- preprocess_data(
    loaded_data = loaded_data,
    files_path = path_to_pkg_files,
    config_filename = pkg_synth_config,
    dict_suv_filename = pkg_synth_dict,
    verbose = FALSE
  )
  # Fix: Use expect_type for lists
  expect_type(preprocessed_output$data, "list")
  expect_type(preprocessed_output$config, "list")
  
  # 3. Test Demographics [cite: 20]
  dset_name <- names(preprocessed_output$data)[1]
  demo_tab <- create_demographics_table(
    processed_dataset_list = preprocessed_output$data[[dset_name]],
    dataset_name = dset_name,
    id_col = "RID"
  )
  expect_s3_class(demo_tab, "data.frame")
  
  # 4. Test Single-Cohort Discovery [cite: 23, 24, 25]
  # Fix: Provide a specific CSV name and use tempdir() correctly
  test_csv <- "test_results.csv"
  
  single_res <- biodiscvr_single(
    dataset_data = preprocessed_output$data[[dset_name]],
    dataset_name = dset_name,
    group = "CI",
    config = preprocessed_output$config,
    var_composition = 1,
    output_csv_name = test_csv,
    output_dir = tempdir(), 
    ga_seed = 42
  )
  
  expect_type(single_res, "list")
  expect_s3_class(single_res$result_row, "data.frame")
  
  
  # 5. Check that a final fitness value exists
  final_fitness <- single_res$result_row$fitness_value
  expect_true(!is.na(final_fitness))
})