## A very simple script, that should be run on an installed
## package, more or less as follows:
##
## R --quiet --vanilla < test_serial.r
##
## This serves to check the namespacing usage for the parallel 
## system, but we check serial too to be safe.

library(harsat)
library(here)
working.directory <- here()


water_data <- read_data(
  compartment = "water",
  purpose = "OSPAR",
  contaminants = "water.txt",
  stations = "stations.txt",
  data_dir = file.path(working.directory, "data", "example_OSPAR"),
  info_dir = file.path(working.directory, "information", "OSPAR_2022"),
  extraction = "2023/08/23"
)

water_data <- tidy_data(water_data)

water_timeseries <- create_timeseries(
  water_data,
  determinands.control = list(
    CHR = list(det = "CHRTR", action = "replace"),
    BBKF = list(det = c("BBF", "BKF", "BBJF", "BBJKF"), action = "bespoke"),
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum"),
    CB138 = list(det = c("CB138+163"), action = "replace"),
    HCEPX = list(det = c("HCEPC", "HCEPT"), action = "sum"), 
    HCH = list(det = c("HCHA", "HCHB", "HCHG"), action = "sum")
  )
)

water_assessment_serial <- run_assessment(
  water_timeseries, 
  AC = "EQS", 
  parallel = FALSE
)

water_assessment_parallel <- run_assessment(
  water_timeseries, 
  AC = "EQS", 
  parallel = TRUE
)
