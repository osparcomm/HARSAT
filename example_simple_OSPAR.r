# Intro ----

# This is a very simple example that assesses a small subset of the biota data
# used in the OSPAR 2022 CEMP assessment


# Setup ----

# Sources functions (folder R) and reference tables (folder information)
# The functions and reference tables folders are assumed to be in the current
# R project folder


rm(list = objects())

devtools::load_all()


# Read data from ICES extraction ----

# There are three input data sets: 
# - the contaminant data
# - the station dictionary
# - the quality assurance data (more accurately called a chemical methods file);
#     this will disappear before release

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "OSPAR",                               
  contaminants = "test_data.csv", 
  stations = "station_dictionary.csv", 
  QA = "quality_assurance.csv",
  data_dir = file.path("data", "example_simple_OSPAR"),
  info_files = list(
    determinand = "determinand_simple_OSPAR.csv", 
    thresholds = "thresholds_biota_simple_OSPAR.csv"
  ),
  info_dir = "information", 
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

biota_assessment <- run_assessment(
  biota_timeSeries, 
  AC = c("BAC", "EAC", "EQS", "HQS")
)


# can supply own function for calculating AC - in the example below it
# will generate exactly the same results

# my_get_AC <- get_AC$biota
# 
# biota_assessment <- ctsm_assessment(
#   biota_timeSeries, 
#   AC = c("BAC", "EAC", "EQS", "HQS"), 
#   get_AC_fn = my_get_AC
# )



# check convergence - no errors this time

check_assessment(biota_assessment)


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
    "EQS" = "green",
    "HQS" = "green"
  ),
  above = c(
    "BAC" = "orange", 
    "EAC" = "red", 
    "EQS" = "red",
    "HQS" = "red"
  ), 
  none = "black"
)

write_summary_table(
  biota_assessment, 
  determinandGroups = webGroups,
  classColour = classColour,
  collapse_AC = list(EAC = c("EAC", "EQS")),
  output_dir = file.path("output", "example_simple_OSPAR"), 
)


