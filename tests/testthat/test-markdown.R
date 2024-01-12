library(rprojroot)
working_directory <- is_testthat$find_file()

data_dir = file.path(working_directory, "..", "datasets", "external-1", "data")
info_dir = file.path(working_directory, "..", "datasets", "external-1", "information")

library(rmarkdown)

package_dir = system.file(package = "harsat")
template_dir = file.path(package_dir, "markdown")

output_dir <- file.path(tempdir(), 'report-tests')
dir.create(output_dir, showWarnings = FALSE)

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

test_that('markdown reporting works', {

  biota_data <- call_read_data()
  biota_data <- tidy_data(biota_data)
  biota_timeseries <- create_timeseries(biota_data, get_basis = get_basis_most_common)
  biota_assessment <- run_assessment(biota_timeseries)

  report_file <- file.path(template_dir, "report_assessment.Rmd")
  output_file <- tempfile("plot", tmpdir = output_dir, fileext = c(".htm"))

  params = list(
    assessment_object = biota_assessment,
    series = 'A902 PFOA Globicephala melas LI JV'
  )

  rmarkdown::render(
    report_file,
    output_file = output_file,
    params = params,
    envir = new.env()
  )

  file.body <- readr::read_file(output_file)
  
  expect_true(base::grepl('Perfluorooctanoic', file.body, fixed = TRUE))
  expect_true(base::grepl('long-finned pilot whale', file.body, fixed = TRUE))
  expect_true(base::grepl('liver', file.body, fixed = TRUE))
  expect_true(base::grepl('Føroyskar hvalvágir', file.body, fixed = TRUE))

  expect_true(base::grepl('Timeseries metadata', file.body, fixed = TRUE))
  expect_true(base::grepl('Assessment plot', file.body, fixed = TRUE))
  expect_true(base::grepl('Trend with data', file.body, fixed = TRUE))
  expect_true(base::grepl('Statistical analysis', file.body, fixed = TRUE))

  ## In addition, we should also be able to run `report_assessment`

  report_assessment(
    biota_assessment,
    subset = NULL,
    output_dir = output_dir
  )

  report.files <- list.files(output_dir, pattern = '\\.html$')
  expect_setequal(report.files, c(
    "A902 Faroe Islands Føroyskar hvalvágir PFOA Globicephala melas BB JV.html",
    "A902 Faroe Islands Føroyskar hvalvágir PFOA Globicephala melas LI JV.html",
    "A902 Faroe Islands Føroyskar hvalvágir PFOS Globicephala melas BB JV.html",
    "A902 Faroe Islands Føroyskar hvalvágir PFOS Globicephala melas LI JV.html"
  ))
})
