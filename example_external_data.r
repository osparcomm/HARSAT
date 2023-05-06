# Intro ----

# An example to demonstrate the use of external data
# Mercury is the only determinand assessed


# Setup ----

# Sources functions (folder functions) and reference tables (folder information)
# Only those functions needed for this example are sourced here
# The functions and reference tables folders are assumed to be in the current
# R project folder

rm(list = objects())

function_path <- file.path("functions")

source(file.path(function_path, "import_functions.R"))
source(file.path(function_path, "import_check_functions.R"))
source(file.path(function_path, "import_external_data.R"))
source(file.path(function_path, "assessment_functions.R"))
source(file.path(function_path, "ctsm_lmm.R"))
source(file.path(function_path, "reporting_functions.R"))
source(file.path(function_path, "support_functions.R"))
source(file.path(function_path, "graphics_functions.R"))


# source reference tables and associated information functions

info_AC_type <- "EXTERNAL"
if(tolower(info_AC_type) != tolower("OSPAR") && tolower(info_AC_type) != tolower("HELCOM") && tolower(info_AC_type) != tolower("EXTERNAL")){
  stop("info_AC_type can only be OSPAR, HELCOM, or EXTERNAL")
}

source(file.path(function_path, "information_functions.R"))

info.assessment.criteria <- ctsm_read_assessment_criteria(
  list(biota = "assessment criteria biota.csv")
)


# Read data ----

# mercury data with supporting variables and station dictionary

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "AMAP",
  contaminants = "AMAP_external_data_new_data_only_CAN_MarineMammals.csv", 
  stations = "AMAP_external_new_stations_only.csv", 
  data_path = file.path("data", "example_external_data"), 
  data_format = "external",
  info_path = "information",
  control = list(region_id = "AMAP_region")
)  


# Prepare data for next stage ----

# get correct variables and streamline the data files

biota_data <- ctsm_tidy_data(biota_data)



# Construct timeseries ----

# uses OSPAR basis choice for mercury

biota_timeSeries <- ctsm_create_timeSeries(
  biota_data,
  get_basis = get_basis_biota_OSPAR
)


# Assessment ----

## main runs ----

biota_assessment <- ctsm_assessment(biota_timeSeries)


# use the code below if it takes a long time to run

# biota_assessment <- ctsm_assessment(biota_timeSeries, parallel = TRUE)


## check convergence ----

(wk_check <- ctsm_check_convergence(biota_assessment$assessment))


# Summary files ----

ctsm_summary_table(
  biota_assessment,
  output_dir = file.path("output", "example_external_data")
)


# Graphics output ----

# plots assessment with either data (file_type = "data") or annual index 
# (file_type = "index") or both (default)

# output can be png or pdf

# can subset assessment based on variables in either timeSeries or stations 
# components of object: commonly by determinand, matrix, species, station_code 
# or station_name; can also use the series identifier in row.names(timeSeries)
# if subset is NULL (default), all timeseries are plotted (can take some time)

ctsm_plot_assessment(
  biota_assessment,
  subset = species %in% "Phoca hispida",
  output_dir = file.path("output", "graphics"), 
  file_type = "data",
  file_format = "png"
)

ctsm_plot_assessment(
  biota_assessment,
  subset = matrix %in% "LI",
  output_dir = file.path("output", "graphics"), 
  file_type = "index",
  file_format = "pdf"
)

ctsm_plot_assessment(
  biota_assessment,
  subset = station_code %in% "A1",
  output_dir = file.path("output", "graphics"), 
  file_type = "data",
  file_format = "pdf"
)

ctsm_plot_assessment(
  biota_assessment, 
  subset = series == "A1 HG Phoca hispida LI adult",
  output_dir = file.path("output", "graphics"), 
  file_format = "pdf"
)

