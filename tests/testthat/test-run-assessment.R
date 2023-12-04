library(rprojroot)
working_directory <- is_testthat$find_file()

data_dir = file.path(working_directory, "..", "datasets", "external-1", "data")
info_dir = file.path(working_directory, "..", "datasets", "external-1", "information")

call_read_data <- function(compartment = "biota", data_format = "external", purpose = "AMAP",
                           contaminants = "EXTERNAL_FO_PW_DATA_TINY.csv") {
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

test_that('run_assessment works', {
    biota_data <- call_read_data()
    biota_data <- tidy_data(biota_data)
    biota_timeseries <- create_timeseries(biota_data, get_basis = get_basis_most_common)
    biota_assessment <- run_assessment(biota_timeseries)

    expect_equal(biota_assessment$info$compartment, 'biota')
    expect_type(biota_assessment$timeSeries, 'list')
    expect_type(biota_assessment$assessment, 'list')
    expect_type(biota_assessment$assessment$`A902 PFOA Globicephala melas LI JV`, 'list')
    expect_equal(biota_assessment$assessment$`A902 PFOA Globicephala melas LI JV`$method, 'linear')
    expect_equal(biota_assessment$assessment$`A902 PFOA Globicephala melas LI JV`$convergence, 0)
    expect_equal(biota_assessment$assessment$`A902 PFOA Globicephala melas BB JV`$method, 'none')
    expect_equal(biota_assessment$assessment$`A902 PFOA Globicephala melas BB JV`$convergence, 0)
})
