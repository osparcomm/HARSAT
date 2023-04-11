# Intro ----

# An example to demonstrate the use of external data
# It is based on an assessment of mercury concentrations with some data 
# provided in an ICES extraction and other data provided from an external file


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

info_species_file_id <- "species_2020.csv"

info_AC_type <- "EXTERNAL"
if(tolower(info_AC_type) != tolower("OSPAR") && tolower(info_AC_type) != tolower("HELCOM") && tolower(info_AC_type) != tolower("EXTERNAL")){
  stop("info_AC_type can only be OSPAR, HELCOM, or EXTERNAL")
}

info_AC_infile <- list(
  biota = "assessment criteria biota.csv",
  sediment = "assessment criteria sediment.csv",
  water = "assessment criteria water.csv"
)

info_determinand_infile <- "determinand_external_data.csv"

source(file.path(function_path, "information_functions.R"))

# check if species values are within pre-defined range
# constants
min_value = 0L # min value for species values range check
max_value = 100L # max value for species values range check
species <- read.csv(file.path("information",info_species_file_id),header=TRUE,colClasses="character")
values_range_check_species(species, min_value, max_value)


# Read data from ICES extraction ----

# mercury data with supporting variables, station dictionary, and 
# quality assurance (methods) file

biota_data <- ctsm_read_data(
  compartment = "biota", 
  # purpose = "OSPAR",                               
  contaminants = "AMAP_external_data_new_data_only_CAN_MarineMammals.csv", 
  stations = "AMAP_external_new_stations_only.csv", 
  path = file.path("data", "example_external_data"), 
  extraction = "2022/01/11",
  max_year = 2020L, 
  data_format = "external", 
  control = list(region_id = "AMAP_region")
)  

# saveRDS(biota_data, file.path("RData", "biota data.rds"))


# Adjust ICES extraction and import external AMAP data ----

# corrects known errors in data
# makes adjustments to data that can't be done easily in the core code

# import external AMAP data 
# uses functions in 'functions/import external data.r' 
# the functions need to be generalised to deal with more determinands and to 
# deal with the case where there is no ICES extraction; i.e. only external data

# the markdown also resolves some inconsistencies where the 'same data' are in
# both the ICES extraction and the external data file

# the report of the adjustments can be sent anywhere

if(info_AC_type != "EXTERNAL") {
  rmarkdown::render(
    "example_external_data_adjustments.Rmd", 
    output_file = "mercury_adjustments.html",
    output_dir = file.path("output", "example_external_data") 
  )
}

# saveRDS(biota_data, file.path("RData", "biota data adjusted.rds"))


# Prepare data for next stage ----

# gets correct variables and streamlines some of the data files

biota_data <- ctsm_tidy_data(biota_data)



# Construct timeseries ----

if (info_AC_type == "EXTERNAL") {
  biota_timeSeries <- ctsm_create_timeSeries(
    biota_data,
    get_basis = get_basis_biota_OSPAR
  )
}

# only need the liipid weight control structure if lipidwt% has been reported in 
# multiple ways - typically an ICES extraction issue

if(info_AC_type != "EXTERNAL") {
  biota_timeSeries <- ctsm_create_timeSeries(
    biota_data,
    determinands.control = list(
      "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
    ),
    get_basis = get_basis_biota_OSPAR
  )
}  


# Assessment ----

## main runs ----

biota_assessment <- ctsm.assessment.setup(
  biota_timeSeries, 
  recent.trend = 20
)

biota_assessment$assessment <- ctsm.assessment(biota_assessment)


# use the code below if it takes a long time to run

# biota_assessment$assessment <- ctsm.assessment(biota_assessment, parallel = TRUE)



## check convergence ----

(wk_check <- ctsm_check_convergence(biota_assessment$assessment))


# this time series has missing standard errors   
# "3371 HG Mytilus edulis SB Not_applicable"                       

# refit with different numerical differencing arguments
# biota_assessment$assessment[wk_check] <- 
#   ctsm.assessment(biota_assessment, seriesID = wk_check, hess.d = 0.01, hess.r = 8)

# saveRDS(biota_assessment, file.path("RData", "biota assessment.rds"))


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
# or station_name
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


