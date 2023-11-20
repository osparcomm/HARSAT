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
  
test_that('create timeseries', {

  biota_data <- call_read_data()

  biota_data <- tidy_data(biota_data)

  biota_timeseries <- create_timeseries(biota_data, get_basis = get_basis_most_common)

  expect_type(biota_timeseries, 'list')

  expect_type(biota_timeseries$info, 'list')
  expect_type(biota_timeseries$data, 'list')
  expect_equal(biota_timeseries$info$compartment, 'biota')
  expect_equal(biota_timeseries$info$purpose, 'AMAP')
  expect_equal(biota_timeseries$info$data_format, 'external')
  expect_equal(biota_timeseries$info$max_year, 2020)

  expect_type(biota_timeseries$info$determinand, 'list')
  expect_type(biota_timeseries$info$species, 'list')
  expect_type(biota_timeseries$info$thresholds, 'list')
  expect_type(biota_timeseries$info$matrix, 'list')

  expect_type(biota_timeseries$stations, 'list')

  expect_type(biota_timeseries$timeSeries, 'list')
})
