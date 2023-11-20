library(rprojroot)
working_directory <- is_r_package$find_file()

data_dir = file.path(working_directory, "tests", "datasets", "external-1")
info_dir = file.path(working_directory, "information", "AMAP")

call_read_data <- function(compartment = "biota", data_format = "external", purpose = "AMAP",
                           contaminants = "EXTERNAL_FO_PW_DATA.csv") {
  read_data(
    compartment = compartment,
    purpose = purpose,
    contaminants = contaminants,
    stations = "EXTERNAL_AMAP_STATIONS.csv",
    data_dir = data_dir,
    data_format = data_format,
    info_dir = info_dir
  )
}
  
test_that('create timeseries with missing species', {

  biota_data <- call_read_data(contaminants = "EXTERNAL_FO_PW_DATA_missing_species.csv")

  biota_data <- tidy_data(biota_data)

  # Check the default for the missing_species control
  expect_equal(biota_data$info$missing_species, 'error')

  expect_error(create_timeseries(biota_data, get_basis = get_basis_most_common),
    "The following values are not in the species reference table.\nPlease add them to the reference table or edit your data to continue.\nAnlobicephala melas"
  )

  # If, on the other hand, the control value is TRUE, should pass and remove data
  #biota_data$info$missing_species <- 'warning'

  #biota_timeSeries <- create_timeseries(biota_data, get_basis = get_basis_most_common)
})
