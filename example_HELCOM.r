# Intro ----

# Essentially this replicates the HELCOM 2022 assessment
# The biota data file has been trimmed to manage its size, and the 'initial' 
# data, unique to HELCOM, has not been implemented


# Setup ----

## functions and reference tables ----

# Sources functions (folder functions) and reference tables (folder information)
# Only those functions needed for this example are sourced here
# The functions and reference tables folders are assumed to be in the current
# R project folder

# Test comment

rm(list = objects())

function_path <- file.path("functions")

source(file.path(function_path, "import_functions.R"))
source(file.path(function_path, "import_check_functions.R"))
source(file.path(function_path, "assessment_functions.R"))
source(file.path(function_path, "ctsm_lmm.R"))
source(file.path(function_path, "proportional_odds_functions.R"))
source(file.path(function_path, "imposex_functions.R"))
source(file.path(function_path, "imposex_clm.R"))
source(file.path(function_path, "reporting_functions.R"))
source(file.path(function_path, "support_functions.R"))


# source reference tables and associated information functions

info_AC_type <- "HELCOM"

source(file.path(function_path, "information_functions.R"))

info.determinand <- ctsm_read_determinand("determinand_HELCOM.csv")

info.species <- ctsm_read_species("species_HELCOM_2023.csv")

info.assessment.criteria <- ctsm_read_assessment_criteria(
  list(
    biota = "assessment criteria biota HELCOM.csv",
    sediment = "assessment criteria sediment HELCOM.csv",
    water = "assessment criteria water.csv"
  )
)



# Read data and make adjustments ----

## biota ----

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "HELCOM",                               
  contaminants = "biota_data.csv", 
  stations = "station_dictionary.csv", 
  QA = "quality_assurance.csv",
  data_path = file.path("data", "example_HELCOM"), 
  extraction = "2022/10/06",
  max_year = 2021L
)  

# saveRDS(biota_data, file.path("RData", "biota data.rds"))


## sediment ----

sediment_data <- ctsm_read_data(
  compartment = "sediment",
  purpose = "HELCOM",
  contaminants = file.path("example_HELCOM_new_format", "sediment_data.csv"),
  stations = file.path("example_HELCOM", "station_dictionary.csv"),
  data_path = "data",
  extraction = "2022/10/06",
  max_year = 2021L,
  data_format = "ICES_new"
)

# saveRDS(sediment_data, file = file.path("RData", "sediment data.rds"))


## water ----

water_data <- ctsm_read_data(
  compartment = "water", 
  purpose = "HELCOM",                               
  contaminants = file.path("example_HELCOM_new_format", "water_data.csv"), 
  stations = file.path("example_HELCOM", "station_dictionary.csv"), 
  data_path = "data", 
  extraction = "2022/10/06",
  max_year = 2021L, 
  data_format = "ICES_new"
)  

# saveRDS(water_data, file.path("RData", "water data.rds"))


## adjustments ----

# correct known errors in the data

rmarkdown::render(
  "example_HELCOM_data_adjustments.Rmd", 
  output_file = "HELCOM_adjustments.html",
  output_dir = file.path("output", "example_HELCOM") 
)

# saveRDS(sediment_data, file = file.path("RData", "sediment data adjusted.rds"))
# saveRDS(water_data, file = file.path("RData", "water data adjusted.rds"))


# Prepare data for next stage ----

# gets correct variable and streamlines some of the data files

biota_data <- ctsm_tidy_data(biota_data)
sediment_data <- ctsm_tidy_data(sediment_data)
water_data <- ctsm_tidy_data(water_data)



# Construct timeseries ----

## biota ----

# ad-hoc change to merge MU and MU&EP data for organics for Finnish perch
# only need to do this for years up to and including 2013 when there are 
# no MU&EP measurements, so no risk of mixing up MU and MU&EP data.

biota_data$data <- left_join(
  biota_data$data, 
  biota_data$stations[c("station_code", "country")],
  by = "station_code"
)
  
biota_data$data <- mutate(
  biota_data$data,
  .id = country == "Finland" & 
    species == "Perca fluviatilis" &  
    year <= 2013 & 
    !(determinand %in% c("CD", "HG", "PB", "PFOS")),
  matrix = if_else(
    .id & !(determinand %in% c("DRYWT%", "FATWT%")), 
    "MU&EP", 
    matrix
  ) 
)

wk <- biota_data$data %>% 
  filter(.id & determinand %in% c("DRYWT%", "FATWT%")) %>% 
  mutate(
    replicate = max(biota_data$data$replicate) + 1:n(),
    matrix = "MU&EP",
    .id = NULL
  )

biota_data$data <- mutate(biota_data$data, .id = NULL)

biota_data$data <- bind_rows(biota_data$data, wk)


# ad-hoc change to merge methods of analysis for Poland for PYR10H

biota_data$data <- mutate(
  biota_data$data, 
  method_analysis = if_else(
    alabo %in% "IMWP" & 
      determinand %in% "PYR1OH" &
      year %in% 2020:2021,
    "HPLC-FD", 
    method_analysis
  )
)  


# ad_hoc change to info_TEQ to make it appropriate for human health QS

info_TEQ["CDFO"] <- 0.0003

biota_data$data$country <- NULL


biota_timeSeries <- ctsm_create_timeSeries(
  biota_data,
  determinands.control = list(
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum"),
    SBDE6 = list(
      det = c("BDE28", "BDE47", "BDE99", "BD100", "BD153", "BD154"), 
      action = "sum"
    ),
    HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum"),
    CB138 = list(det = "CB138+163", action = "replace"),
    SCB6 = list(
      det = c("CB28", "CB52", "CB101", "CB138", "CB153", "CB180"), 
      action = "sum"
    ),
    TEQDFP = list(det = names(info_TEQ), action = "bespoke"),
    VDS = list(det = "VDSI", action = "bespoke"), 
    INTS = list(det = "INTSI", action = "bespoke"),
    "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
  ),
  normalise = ctsm_normalise_biota_HELCOM,
  normalise.control = list(
    lipid = list(method = "simple", value = 5), 
    other = list(method = "none") 
  )
)


# resolve Finnish perch changes

biota_timeSeries$data <- mutate(
  biota_timeSeries$data, 
  matrix = if_else(
    year <= 2013 & matrix == "MU&EP", 
    "MU", 
    matrix
  )
)

# resolve Polish metoa changes

biota_timeSeries$data <- left_join(
  biota_timeSeries$data, 
  biota_data$stations[c("station_code", "country")],
  by = "station_code"
)


biota_timeSeries$data <- mutate(
  biota_timeSeries$data, 
  method_analysis = if_else(
    country == "Poland" &  
      determinand %in% "PYR1OH" &
      year %in% 2020,
    "HPLC-ESI-MS-MS", 
    method_analysis
  ), 
  method_analysis = if_else(
    country == "Poland" &  
      determinand %in% "PYR1OH" &
      year %in% 2021,
    "GC-MS-MS", 
    method_analysis
  ), 
)  

biota_timeSeries$data$country <- NULL

# saveRDS(biota_timeSeries, file.path("RData", "biota timeSeries.rds"))


## sediment ----

sediment_timeSeries <- ctsm_create_timeSeries(
  sediment_data,
  determinands.control = list(
    SBDE6 = list(
      det = c("BDE28", "BDE47", "BDE99", "BD100", "BD153", "BD154"), 
      action = "sum"
    ),
    HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum")
  ),
  normalise = ctsm_normalise_sediment_HELCOM,
  normalise.control = list(
    metals = list(method = "pivot", normaliser = "AL"), 
    copper = list(method = "hybrid", normaliser = "CORG", value = 5),
    organics = list(method = "simple", normaliser = "CORG", value = 5) 
  )
)

# saveRDS(sediment_timeSeries, file = file.path("RData", "sediment timeSeries.rds"))


## water ----

water_timeSeries <- ctsm_create_timeSeries(
  water_data,
  determinands.control = list(
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum")
  )
)

# saveRDS(water_timeSeries, file.path("RData", "water timeSeries.rds"))


# Assessment ----

## sediment ----

### main runs ----

sediment_assessment <- ctsm_assessment(
  sediment_timeSeries, AC = "EQS", parallel = TRUE)


### check convergence ----

ctsm_check_convergence(sediment_assessment$assessment)


# saveRDS(sediment_assessment, file.path("RData", "sediment assessment.rds"))



## biota ----

### main runs ----

# preliminary analysis required for imposex assessment 
# takes a long time to run!!!!!
# I need to turn this into a function

source("example_HELCOM_imposex_preparation.R")


# can sometimes be useful to split up the assessment because of size limitations
# not really needed here, but done to illustrate

wk_determinands <- ctsm_get_determinands("biota")
wk_group <- info.determinand[wk_determinands, "biota_group"]

biota_assessment <- ctsm_assessment(
  biota_timeSeries, 
  AC = c("BAC", "EAC", "EQS", "MPC"), 
  subset = determinand %in% wk_determinands[wk_group == "Metals"], 
  parallel = TRUE
)

wk_organics <- c(
  "PAH_parent", "PBDEs", "Organobromines", "Organofluorines", 
  "Chlorobiphenyls", "Dioxins"
)  

biota_assessment <- ctsm_update_assessment(
  biota_assessment, 
  subset = determinand %in% wk_determinands[wk_group %in% wk_organics], 
  parallel = TRUE
)

biota_assessment <- ctsm_update_assessment(
  biota_assessment, 
  subset = determinand %in% wk_determinands[wk_group %in% "Metabolites"]
)

biota_assessment <- ctsm_update_assessment(
  biota_assessment, 
  subset = determinand %in% wk_determinands[wk_group %in% "Imposex"]
)


### check convergence ----

# no checks for Imposex

wk_id <- biota_assessment$timeSeries$determinand %in% 
  wk_determinands[wk_group != "Imposex"] 
wk_id <- row.names(biota_assessment$timeSeries)[wk_id]
(wk_check <- ctsm_check_convergence(biota_assessment$assessment[wk_id]))


# two time series need to be refitted

# "2109 PB Perca fluviatilis MU" - fixed bounds
biota_assessment <- ctsm_update_assessment( 
  biota_assessment, series == wk_check[1], fixed_bound = 20
)

# "2299 PYR1OH Limanda limanda BI HPLC-FD" - standard errors
biota_assessment <- ctsm_update_assessment( 
  biota_assessment, series == wk_check[2], hess.d = 0.0001, hess.r = 8
)


# check it has worked

ctsm_check_convergence(biota_assessment$assessment[wk_id])


# saveRDS(biota_assessment, file.path("RData", "biota assessment.rds"))


## water ----

### main runs ----

water_assessment <- ctsm_assessment(
  water_timeSeries, AC = "EQS", parallel = TRUE
)


### check convergence ----

(wk_check <- ctsm_check_convergence(water_assessment$assessment))

# refit a couple of time series

# "5190 CD Yes" - fixed effects on bounds
# "5192 CD Yes" - fixed effects on bounds

water_assessment <- ctsm_update_assessment(
  water_assessment, series %in% wk_check, fixed_bound = 20
)


ctsm_check_convergence(water_assessment$assessment)


# saveRDS(water_assessment, file.path("RData", "water assessment.rds"))


# Summary files ----

## web objects ----

webGroups <- list(
  levels = c(
    "Metals", "Organotins", 
    "PAH_parent", "Metabolites", 
    "PBDEs", "Organobromines", 
    "Organofluorines", 
    "Chlorobiphenyls", "Dioxins", 
    "Imposex" 
  ),  
  labels = c(
    "Metals", "Organotins", 
    "PAH parent compounds", "PAH metabolites", 
    "Polybrominated diphenyl ethers", "Organobromines (other)", 
    "Organofluorines", 
    "Polychlorinated biphenyls", "Dioxins", 
    "Imposex"
  )
)

ctsm_summary_table(
  biota_assessment,
  determinandGroups = webGroups,
  classColour = list(
    below = c("BAC" = "green", "EAC" = "green", "EQS" = "green", "MPC" = "green"),
    above = c("BAC" = "red", "EAC" = "red", "EQS" = "red", "MPC" = "red"),
    none = "black"
  ),
  collapse_AC = list(EAC = c("EAC", "EQS", "MPC")),
  output_dir = file.path("output", "example_HELCOM")
)

ctsm_summary_table(
  sediment_assessment,
  determinandGroups = webGroups,
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_HELCOM")
)

ctsm_summary_table(
  water_assessment,
  determinandGroups = webGroups,
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_HELCOM")
)

