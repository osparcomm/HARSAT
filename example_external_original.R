# Intro ----

# Original example to demonstrate the use of external data
# It is based on an assessment of mercury concentrations with some data 
# provided in an ICES extraction and other data provided from an external file


# Setup ----

# Sources functions (folder R) and reference tables (folder information)
# Only those functions needed for this example are sourced here
# The functions and reference tables folders are assumed to be in the current
# R project folder

rm(list = objects())

function_path <- file.path("R")

source(file.path(function_path, "import_functions.R"))
source(file.path(function_path, "import_check_functions.R"))
source(file.path(function_path, "import_external_data.R"))
source(file.path(function_path, "assessment_functions.R"))
source(file.path(function_path, "ctsm_lmm.R"))
source(file.path(function_path, "reporting_functions.R"))
source(file.path(function_path, "support_functions.R"))
source(file.path(function_path, "information_functions.R"))




# Read data from ICES extraction ----

# mercury data with supporting variables, station dictionary, and 
# quality assurance (methods) file

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "OSPAR",                               
  contaminants = "mercury_data.csv", 
  stations = "station_dictionary.csv", 
  QA = "quality_assurance.csv",
  data_path = file.path("data", "example_external_data"),
  info_files = list(determinand = "determinand_external_data.csv"),
  info_path = "information",
  extraction = "2022/01/11",
  max_year = 2020L)  

# saveRDS(biota_data, file.path("RData", "biota data.rds"))


# Adjust ICES extraction and import external AMAP data ----

# corrects known errors in data
# makes adjustments to data that can't be done easily in the core code

# import external AMAP data 
# uses functions in 'R/import external data.r' 
# the functions need to be generalised to deal with more determinands and to 
# deal with the case where there is no ICES extraction; i.e. only external data

# the markdown also resolves some inconsistencies where the 'same data' are in
# both the ICES extraction and the external data file

# the report of the adjustments can be sent anywhere

rmarkdown::render(
  "example_external_data_adjustments.Rmd", 
  output_file = "mercury_adjustments.html",
  output_dir = file.path("output", "example_external_original") 
)

# saveRDS(biota_data, file.path("RData", "biota data adjusted.rds"))


# Prepare data for next stage ----

# gets correct variable and streamlines some of the data files

biota_data <- ctsm_tidy_data(biota_data)



# Construct timeseries ----

biota_timeSeries <- ctsm_create_timeSeries(
  biota_data,
  determinands.control = list(
    "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
  ),
  get_basis = get_basis_biota_OSPAR
)

# saveRDS(biota_timeSeries, file.path("RData", "biota timeSeries.rds"))


# Assessment ----

## main runs ----

biota_assessment <- ctsm_assessment(
  biota_timeSeries, 
  AC = c("BAC", "EQS.OSPAR", "HQS"),
  get_assessment_criteria = get.AC.OSPAR,
  parallel = TRUE 
)


## check convergence ----

(wk_check <- check_convergence_lmm(biota_assessment$assessment))

# this time series has missing standard errors   
# "3371 HG Mytilus edulis SB Not_applicable"                       

# refit with different numerical differencing arguments
biota_assessment <- ctsm_update_assessment(
  biota_assessment, series == wk_check, hess.d = 0.01, hess.r = 8
)

# saveRDS(biota_assessment, file.path("RData", "biota assessment.rds"))


# Summary files ----

biota_web <- biota_assessment
biota_web$info$AC <- c("BAC", "EQS.OSPAR")

summary_table(
  biota_web,
  determinandGroups = list(levels = "Metals", labels = "Metals"),
  classColour = list(
    below = c("BAC" = "blue", "EQS.OSPAR" = "green"),
    above = c("BAC" = "orange", "EQS.OSPAR" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS.OSPAR"),
  output_dir = file.path("output", "example_external_original")
)



