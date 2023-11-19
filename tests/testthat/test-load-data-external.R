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
  
test_that('reading external data works', {
  biota_data <- call_read_data()

  # We can do more checks here, but this is a decent start
  expect_type(biota_data, 'list')
  expect_type(biota_data$info, 'list')
  expect_type(biota_data$data, 'list')
  expect_equal(biota_data$info$compartment, 'biota')
  expect_equal(biota_data$info$purpose, 'AMAP')
  expect_equal(biota_data$info$data_format, 'external')
  expect_equal(biota_data$info$max_year, 2020)

  expect_type(biota_data$info$determinand, 'list')
  expect_type(biota_data$info$species, 'list')
  expect_type(biota_data$info$thresholds, 'list')
  expect_type(biota_data$info$matrix, 'list')

  # Validate a few options by breaking and testing for errors
  expect_error(call_read_data(compartment = "wibble"), 'should be one of "biota", "sediment", "water"')
  expect_error(call_read_data(data_format = "wibble"), 'should be one of "ICES", "external"')
  expect_error(call_read_data(purpose = "wibble"), 'should be one of "OSPAR", "HELCOM", "AMAP", "custom"')

  # Check we report errors for missing/invalid uncertainty units; see #352
  expect_error(call_read_data(contaminants = "EXTERNAL_FO_PW_DATA_invalid_uncertainty.csv"), 
    'Missing or invalid uncertainty units for specified uncertainty values')
})
