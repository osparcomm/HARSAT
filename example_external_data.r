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

# source reference tables and associated information functions

info_species_file_id <- "species_2020.csv"
info_uncertainty_file_id <- "uncertainty_2020.csv"

info_AC_type <- "OSPAR"

info_AC_infile <- list(
  biota = "assessment criteria biota.csv",
  sediment = "assessment criteria sediment.csv",
  water = "assessment criteria water.csv"
)

info_determinand_infile <- "determinand_external_data.csv"

source(file.path(function_path, "information_functions.R"))


# Read data from ICES extraction ----

# mercury data with supporting variables, station dictionary, and 
# quality assurance (methods) file

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "OSPAR",                               
  contaminants = "mercury_data.txt", 
  stations = "station_dictionary.txt", 
  QA = "quality_assurance.txt",
  path = file.path("data", "example_external_data"), 
  extraction = "2022/01/11",
  max_year = 2020L)  

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

rmarkdown::render(
  "example_external_data_adjustments.Rmd", 
  output_file = "mercury_adjustments.html",
  output_dir = file.path("output", "example_external_data") 
)

# saveRDS(biota_data, file.path("RData", "biota data adjusted.rds"))


# Construct timeseries ----

biota_timeSeries <- ctsm_create_timeSeries(
  biota_data,
  determinands.control = list(
    "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
  )
)

# saveRDS(biota_timeSeries, file.path("RData", "biota timeSeries.rds"))


# Assessment ----

## main runs ----

biota_assessment <- ctsm.assessment.setup(
  biota_timeSeries, 
  AC = c("BAC", "EQS.OSPAR", "HQS"), 
  recent.trend = 20
)


library("parallel")
library("pbapply")

wk.cores <- detectCores()
wk.cluster <- makeCluster(wk.cores - 1)

clusterExport(
  wk.cluster, 
  c("biota_assessment", "negTwiceLogLik", "convert.basis", 
    unique(
      c(objects(pattern = "ctsm*"), 
        objects(pattern = "get*"), 
        objects(pattern = "info*")
      )
    )
  )
)

clusterEvalQ(wk.cluster, {
  library("lme4")
  library("tidyverse")
})  


biota_assessment$assessment <- ctsm.assessment(
  biota_assessment, 
  clusterID = wk.cluster
)

stopCluster(wk.cluster)


## check convergence ----

(wk_check <- ctsm_check_convergence(biota_assessment$assessment))

# this time series has missing standard errors   
# "Ireland_Lee K Estuary HG Mytilus edulis SB Not_applicable"                       

# refit with different numerical differencing arguments
biota_assessment$assessment[wk_check] <- 
  ctsm.assessment(biota_assessment, seriesID = wk_check, hess.d = 0.01, hess.r = 8)


## tidy up ----

# make timeSeries factor levels more interpretable

biota_assessment$timeSeries <- biota_assessment$timeSeries %>% 
  rownames_to_column(".rownames") %>% 
  mutate(
    .matrix = ctsm_get_info("matrix", matrix, "name") %>% 
      str_to_sentence() %>% 
      recode(
        "Erythrocytes (red blood cells in vertebrates)" = "Red blood cells",
        "Egg homogenate of yolk and albumin" = "Egg yolk & albumin"
      ), 
    level6name = case_when(
      level6element %in% "matrix" ~ .matrix,
      level6element %in% "sex" ~ recode(level6name, "F" = "Female", "M" = "Male"),
      TRUE ~ level6name
    ),
    level7name = case_when(
      level7element %in% "matrix" ~ .matrix,
      level7element %in% "AMAP_group" ~ gsub("_", " ", level7name),
      TRUE ~ level7name
    ),
    .matrix = NULL,
    level6element = recode(
      level6element, 
      matrix = "Tissue", 
      sex = "Sex", 
      "METOA" = "Chemical analysis"
    ),
    level7element = recode(
      level7element, 
      matrix = "Tissue", 
      AMAP_group = "Mammal group"
    )
  ) %>% 
  column_to_rownames(".rownames")


# saveRDS(biota_assessment, file.path("RData", "biota assessment.rds"))


# Summary files ----

webGroups = list(
  levels = "Metals", 
  labels = "Metals"
)


# only use environmental AC (not the health quality standard)

biota_web <- biota_assessment
biota_web$info$AC <- c("BAC", "EQS.OSPAR")

biota_web <- ctsm_web_initialise(
  biota_web,
  classColour = list(
    below = c(
      "BAC" = "blue", 
      "EQS.OSPAR" = "green"
    ),
    above = c(
      "BAC" = "orange", 
      "EQS.OSPAR" = "red"
    ), 
    none = "black"
  ),
  determinandGroups = webGroups)

# saveRDS(biota_web, file.path("RData", "biota web.rds"))


# write table 

ctsm_summary_table(
  assessments = list(Biota = biota_web),
  determinandGroups = webGroups,
  path = file.path("output", "example_external_data")
)

