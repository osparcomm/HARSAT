# Intro ----

# Essentially this replicates the HELCOM 2022 assessment
# The biota data file has been trimmed to manage its size, and the 'initial' 
# data, unique to HELCOM, has not been implemented


# Setup ----

## functions and reference tables ----

# Sources functions (folder R) and reference tables (folder information)
# Only those functions needed for this example are sourced here
# The functions and reference tables folders are assumed to be in the current
# R project folder

# Test comment

rm(list = objects())

devtools::load_all()


# Read data and make adjustments ----

## biota ----

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "HELCOM",                               
  contaminants = "biota_data.csv", 
  stations = "station_dictionary.csv", 
  QA = "quality_assurance.csv",
  data_dir = file.path("data", "example_HELCOM"),
  info_dir = "information",
  extraction = "2022/10/06",
  max_year = 2021L
)  



## sediment ----

sediment_data <- ctsm_read_data(
  compartment = "sediment",
  purpose = "HELCOM",
  contaminants = file.path("example_HELCOM_new_format", "sediment_data.csv"),
  stations = file.path("example_HELCOM", "station_dictionary.csv"),
  data_dir = "data",
  info_dir = "information",
  extraction = "2022/10/06",
  max_year = 2021L,
  data_format = "ICES_new"
)



## water ----

water_data <- ctsm_read_data(
  compartment = "water", 
  purpose = "HELCOM",                               
  contaminants = file.path("example_HELCOM_new_format", "water_data.csv"), 
  stations = file.path("example_HELCOM", "station_dictionary.csv"), 
  data_dir = "data", 
  info_dir = "information", 
  extraction = "2022/10/06",
  max_year = 2021L, 
  data_format = "ICES_new"
)  



## adjustments ----

# correct known errors in the data

info_TEQ <- c(
  "CB77" = 0.0001, "CB81" = 0.0003, "CB105" = 0.00003, "CB118" = 0.00003, 
  "CB126" = 0.1, "CB156" = 0.00003, "CB157" = 0.00003, "CB167" = 0.00003, 
  "CB169" = 0.03, "CDD1N" = 1, "CDD4X" = 0.1, "CDD6P" = 0.01, "CDD6X" = 0.1, 
  "CDD9X" = 0.1, "CDDO" = 0.0003, "CDF2N" = 0.3, "CDF2T" = 0.1, "CDF4X" = 0.1, 
  "CDF6P" = 0.01, "CDF6X" = 0.1, "CDF9P" = 0.01,
  "CDF9X" = 0.1, "CDFO" = 0.00003, "CDFP2" = 0.03, "CDFX1" = 0.1, "TCDD" = 1
)


rmarkdown::render(
  "example_HELCOM_data_adjustments.Rmd", 
  output_file = "HELCOM_adjustments.html",
  output_dir = file.path("output", "example_HELCOM") 
)



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



## water ----

water_timeSeries <- ctsm_create_timeSeries(
  water_data,
  determinands.control = list(
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum")
  )
)



# Assessment ----

## sediment ----

### main runs ----

sediment_assessment <- run_assessment(
  sediment_timeSeries, 
  AC = "EQS",
  parallel = TRUE
)


### check convergence ----

check_convergence(sediment_assessment)



## biota ----

### main runs ----

# preliminary analysis required for imposex assessment 
# takes a long time to run!!!!!
# I need to turn this into a function

rmarkdown::render(
  file.path("man", "fragments", "example_HELCOM_imposex_preparation.Rmd"), 
  output_file = "HELCOM_imposex_preparation.html",
  output_dir = file.path("output", "example_HELCOM") 
)


# source("example_HELCOM_imposex_preparation.R")


# can sometimes be useful to split up the assessment because of size limitations
# not really needed here, but done to illustrate

wk_determinands <- ctsm_get_determinands(biota_timeSeries$info)
wk_group <- ctsm_get_info(
  biota_timeSeries$info$determinand, wk_determinands, "biota_group"
)

biota_assessment <- run_assessment(
  biota_timeSeries, 
  AC = c("BAC", "EAC", "EQS", "MPC"), 
  subset = determinand %in% wk_determinands[wk_group == "Metals"],
  parallel = TRUE
)

wk_organics <- c(
  "PAH_parent", "PBDEs", "Organobromines", "Organofluorines", 
  "Chlorobiphenyls", "Dioxins"
)  

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = determinand %in% wk_determinands[wk_group %in% wk_organics], 
  parallel = TRUE
)

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = determinand %in% wk_determinands[wk_group %in% "Metabolites"]
)

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = determinand %in% wk_determinands[wk_group %in% "Imposex"]
)


### check convergence ----

check_assessment(biota_assessment)


# two time series need to be refitted

# "2109 PB Perca fluviatilis MU" - fixed bounds
biota_assessment <- update_assessment( 
  biota_assessment, 
  series == "2109 PB Perca fluviatilis MU", 
  fixed_bound = 20
)

# "2299 PYR1OH Limanda limanda BI HPLC-FD" - standard errors
biota_assessment <- update_assessment( 
  biota_assessment, 
  series == "2299 PYR1OH Limanda limanda BI HPLC-FD", 
  hess.d = 0.0001, hess.r = 8
)


# check it has worked

check_assessment(biota_assessment)



## water ----

### main runs ----

water_assessment <- run_assessment(
  water_timeSeries, 
  AC = "EQS", 
  parallel = TRUE
)


### check convergence ----

check_convergence(water_assessment)


# refit a couple of time series where the fixed effects are on their bounds

# "5190 CD Yes" - fixed effects on bounds
# "5192 CD Yes" - fixed effects on bounds

wk_id <- check_convergence(water_assessment, save_result = TRUE)

water_assessment <- update_assessment(
  water_assessment, 
  series %in% wk_id$not_converged, 
  fixed_bound = 20
)


# check that refitted timeseries have converged

check_convergence(water_assessment)




# Summary files ----

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

write_summary_table(
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

write_summary_table(
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

write_summary_table(
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

