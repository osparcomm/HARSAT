# Intro ----

# This is a very simple example that assesses a small subset of the biota data
# used in the OSPAR 2022 CEMP assessment


# Setup ----

# Sources functions (folder functions) and reference tables (folder information)
# Only those functions needed for this example are sourced here
# The functions and reference tables folders are assumed to be in the current
# R project folder


rm(list = objects())

function_path <- file.path("functions")

source(file.path(function_path, "import_functions.R"))
source(file.path(function_path, "import_check_functions.R"))
source(file.path(function_path, "assessment_functions.R"))
source(file.path(function_path, "ctsm_lmm.R"))
source(file.path(function_path, "reporting_functions.R"))
source(file.path(function_path, "support_functions.R"))

# source reference tables and associated information functions

info_AC_type <- "OSPAR"

source(file.path(function_path, "information_functions.R"))

info.determinand <- ctsm_read_determinand("determinand_simple_OSPAR.csv")

info.species <- ctsm_read_species("species_OSPAR_2022.csv")

info.assessment.criteria <- ctsm_read_assessment_criteria(
  list(biota = "assessment criteria biota.csv")
)


# Read data from ICES extraction ----

# There are three input data sets: 
# - the contaminant data
# - the station dictionary
# - the quality assurance data (more accurately called a chemical methods file)

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "OSPAR",                               
  contaminants = "test_data.csv", 
  stations = "station_dictionary.csv", 
  QA = "quality_assurance.csv",
  data_path = file.path("data", "example_simple_OSPAR"), 
  extraction = "2022/01/11",
  max_year = 2020L  
)  


# Prepare data for next stage ----

# gets correct variable and streamlines some of the data files

biota_data <- ctsm_tidy_data(biota_data)


# Construct timeseries ----

# identifies groups of data that form a coherent timeseries
# also does a lot of data cleaning and processing (creates oddities folder)

biota_timeSeries <- ctsm_create_timeSeries(
  biota_data,
  determinands = c("CD", "CB153", "HBCD","HBCDA", "HBCDG", "PYR1OH"), 
  determinands.control = list(
    HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum"),
    "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
  ), 
  get_basis = get_basis_biota_OSPAR
)

# identical (apart from call) to: 
# 
# ctsm_create_timeSeries(
#   biota_data,
#   determinands = ctsm_get_determinands("biota"),
#   determinands.control = list(
#     HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum"),
#     "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
#   )
# )
# 
# ctsm_create_timeSeries(
#   biota_data,
#   determinands.control = list(
#     HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum"),
#     "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
#   )
# )


# Assessment ----

# do the statistical analysis

biota_assessment <- ctsm_assessment(
  biota_timeSeries, 
  AC = c("BAC", "EAC", "EQS.OSPAR", "HQS")
)


# check convergence - no errors this time

ctsm_check_convergence(biota_assessment$assessment)


# Summary files ----

webGroups <- list(
  levels = c("Metals", "Metabolites", "Organobromines", "Chlorobiphenyls"),  
  labels = c(
    "Metals", "PAH metabolites", "Organobromines",  "Polychlorinated biphenyls"
  )
)

classColour <- list(
  below = c(
    "BAC" = "blue", 
    "EAC" = "green", 
    "EQS.OSPAR" = "green",
    "HQS" = "green"
  ),
  above = c(
    "BAC" = "orange", 
    "EAC" = "red", 
    "EQS.OSPAR" = "red",
    "HQS" = "red"
  ), 
  none = "black"
)

ctsm_summary_table(
  biota_assessment, 
  determinandGroups = webGroups,
  classColour = classColour,
  collapse_AC = list(EAC = c("EAC", "EQS.OSPAR")),
  output_dir = file.path("output", "example_simple_OSPAR"), 
)


