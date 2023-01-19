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

info_species_file_id <- "species_2020.csv"
info_uncertainty_file_id <- "uncertainty_2020.csv"

info_AC_type <- "HELCOM"

info_AC_infile <- list(
  biota = "assessment criteria biota HELCOM.csv",
  sediment = "assessment criteria sediment HELCOM.csv",
  water = "assessment criteria water.csv"
)

info_determinand_infile <- "determinand_HELCOM.csv"

source(file.path(function_path, "information_functions.R"))

# check if species values are within pre-defined range
# constants
min_value = 0L # min value for species values range check
max_value = 100L # max value for species values range check
species <- read.csv(file.path("information",info_species_file_id),header=TRUE,colClasses="character")
values_range_check_species(species, min_value, max_value)

## determinands ----

# determinands <- list(
#   Biota = list(
#     Metals = c("CD", "PB", "HG"),
#     PAH_parent = c("FLU", "BAP"), 
#     Metabolites = "PYR1OH",
#     PBDEs = "SBDE6",
#     Organobromines = "HBCD", 
#     Organofluorines = "PFOS",
#     Chlorobiphenyls = "SCB6",
#     Dioxins = "TEQDFP",
#     Imposex = c("IMPS", "INTS", "VDS")
#   ),  
#   Sediment = list(
#     Metals = c("CD", "CU", "PB"),
#     Organotins = "TBSN+",
#     PAH_parent = c("ANT", "FLU"),
#     PBDEs = "SBDE6",
#     Organobromines = "HBCD" 
#   ),
#   Water = list(
#     Metals = c("CD", "PB"), 
#     Organotins = "TBSN+",
#     Organofluorines = "PFOS"
#   )
# )



# Read data and make adjustments ----

## biota ----

biota_data <- ctsm_read_data(
  compartment = "biota", 
  purpose = "HELCOM",                               
  contaminants = "biota_data.txt", 
  stations = "station_dictionary.txt", 
  QA = "quality_assurance.txt",
  path = file.path("data", "example_HELCOM"), 
  extraction = "2022/10/06",
  max_year = 2021L
)  

# saveRDS(biota_data, file.path("RData", "biota data.rds"))


## sediment ----

sediment_data <- ctsm_read_data(
  compartment = "sediment",
  purpose = "HELCOM",
  contaminants = file.path("data", "example_HELCOM_new_format", "sediment_data.txt"),
  stations = file.path("data", "example_HELCOM", "station_dictionary.txt"),
  extraction = "2022/10/06",
  max_year = 2021L,
  data_format = "new"
)

# sediment_data <- ctsm_read_data(
#   compartment = "sediment",
#   purpose = "HELCOM",
#   contaminants = file.path("data", "example_HELCOM", "sediment_data.txt"),
#   stations = file.path("data", "example_HELCOM", "station_dictionary.txt"),
#   QA = file.path("data", "example_HELCOM", "quality_assurance.txt"),
#   extraction = "2022/10/06",
#   max_year = 2021L,
#   data_format = "old"
# )

# saveRDS(sediment_data, file = file.path("RData", "sediment data.rds"))


## water ----

water_data <- ctsm_read_data(
  compartment = "water", 
  purpose = "HELCOM",                               
  contaminants = file.path("data", "example_HELCOM_new_format", "water_data.txt"), 
  stations = file.path("data", "example_HELCOM", "station_dictionary.txt"), 
  extraction = "2022/10/06",
  max_year = 2021L, 
  data_format = "new"
)  

# saveRDS(water_data, file.path("RData", "water data.rds"))


## adjustments ----

# correct known errors in the data

rmarkdown::render(
  "example_HElCOM_data_adjustments.Rmd", 
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

# biota_data$QA <- mutate(
#   biota_data$QA, 
#   metoa = if_else(
#     alabo %in% "IMWP" & 
#       determinand %in% "PYR1OH" &
#       year %in% 2020:2021,
#     "HPLC-FD", 
#     metoa
#   )
# )  

biota_data$data <- mutate(
  biota_data$data, 
  metoa = if_else(
    alabo %in% "IMWP" & 
      determinand %in% "PYR1OH" &
      year %in% 2020:2021,
    "HPLC-FD", 
    metoa
  )
)  


  
# ad_hoc change to info_TEQ to make it appropriate for human health QS

info_TEQ["CDFO"] <- 0.0003


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

biota_timeSeries$data <- mutate(
  biota_timeSeries$data, 
  metoa = if_else(
    grepl("Poland", station) &  
      determinand %in% "PYR1OH" &
      year %in% 2020,
    "HPLC-ESI-MS-MS", 
    metoa
  ), 
  metoa = if_else(
    grepl("Poland", station) &  
      determinand %in% "PYR1OH" &
      year %in% 2021,
    "GC-MS-MS", 
    metoa
  ), 
)  

# saveRDS(biota_timeSeries, file.path("RData", "biota timeSeries.rds"))


## sediment ----

# ad-hoc change to bring in LOIGN

wk_id <- c(
  ctsm_get_determinands("sediment"), 
  "BDE28", "BDE47", "BDE99", "BD100", "BD153", "BD154", 
  "HBCDA", "HBCDB", "HBCDG"
)
names(wk_id) <- NULL

wk_id <- setdiff(wk_id, c("CD", "PB"))

info.determinand[wk_id, "sediment_auxiliary"] <- 
  paste0(info.determinand[wk_id, "sediment_auxiliary"], ", LOIGN") 

info.uncertainty["LOIGN", c("sediment.sd_constant", "sediment.sd_variable")] <-
  c(0, 0.1)


# and uncertainty for SBDE6 (estimated from OSPAR data)

info.uncertainty["SBDE6", c("sediment.sd_constant", "sediment.sd_variable")] <-
  c(0.0433, 0.1727)


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


sediment_assessment <- ctsm.assessment.setup(
  sediment_timeSeries, 
  AC = "EQS", 
  recent.trend = 20
)

library("parallel")
library("pbapply")

wk.cores <- detectCores()
wk.cluster <- makeCluster(wk.cores - 1)

clusterExport(
  wk.cluster, 
  c("sediment_assessment", "negTwiceLogLik",
    unique(c(
      objects(pattern = "ctsm*"), 
      objects(pattern = "get*"), 
      objects(pattern = "info*")
    ))
  )
)

clusterEvalQ(wk.cluster, {
  library("lme4")
  library("tidyverse")
})  


sediment_assessment$assessment <- ctsm.assessment(
  sediment_assessment, 
  clusterID = wk.cluster
)

stopCluster(wk.cluster)


### check convergence ----

ctsm_check_convergence(sediment_assessment$assessment)


# saveRDS(sediment_assessment, file.path("RData", "sediment assessment.rds"))



## biota ----

### main runs ----

biota_assessment <- ctsm.assessment.setup(
  biota_timeSeries, 
  AC = c("BAC", "EAC", "EQS", "MPC"), 
  recent.trend = 20
)


# preliminary analysis required for imposex assessment 
# takes a long time to run!!!!!
# I need to turn this into a function

source("example_HELCOM_imposex_preparation.R")


# can sometimes be useful to split up the assessment because of size limitations
# not really needed here, but done to illustrate

library("parallel")
library("pbapply")

wk.cores <- detectCores()
wk.cluster <- makeCluster(wk.cores - 1)

clusterExport(
  wk.cluster, 
  c("biota_assessment", "negTwiceLogLik", "convert.basis", 
    "assess_imposex", "imposex.assess.index", "imposex_class", "imposex.family", 
    "cuts6.varmean", "biota.VDS.cl", "biota.VDS.estimates", "imposex_assess_clm", 
    "imposex.clm.fit", "imposex.clm.X", "imposex.clm.loglik.calc", 
    "imposex.VDS.p.calc", "imposex.clm.predict", "imposex.clm.cl",
    "imposex.clm.contrast", 
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


wk_determinands <- ctsm_get_determinands("biota")
wk_group <- info.determinand[wk_determinands, "biota_group"]

biota_Metals <- ctsm.assessment(
  biota_assessment, 
  determinandID = wk_determinands[wk_group == "Metals"], 
  clusterID = wk.cluster
)


wk_organics <- c(
  "PAH_parent", "PBDEs", "Organobromines", "Organofluorines", 
  "Chlorobiphenyls", "Dioxins"
)  

biota_Organics <- ctsm.assessment(
  biota_assessment, 
  determinandID = wk_determinands[wk_group %in% wk_organics], 
  clusterID = wk.cluster
)


biota_Metabolites <- ctsm.assessment(
  biota_assessment, 
  determinandID = wk_determinands[wk_group %in% "Metabolites"], 
  clusterID = wk.cluster
)


biota_Imposex <- ctsm.assessment(
  biota_assessment, 
  determinandID = wk_determinands[wk_group %in% "Imposex"],
  clusterID = wk.cluster
)

stopCluster(wk.cluster)




### check convergence ----

# no checks for Imposex

wk_group <- c("Metals", "Organics", "Metabolites")

(wk_check <- sapply(
  wk_group,
  simplify = FALSE, FUN = function(id) {
    assessment <- get(paste("biota", id, sep = "_"))
    ctsm_check_convergence(assessment)
  }))


# two time series need to be refitted

# Germany_OMSH PB Perca fluviatilis MU" - fixed bounds
(wk_id <- wk_check$Metals[1])
biota_Metals[wk_id] <- 
  ctsm.assessment(biota_assessment, seriesID = wk_id, fixed_bound = 20)

# "Germany_FOE-B01 PYR1OH Limanda limanda BI HPLC-FD" - standard errors
(wk_id <- wk_check$Metabolites[1])
biota_Metabolites[wk_id] <- 
  ctsm.assessment(biota_assessment, seriesID = wk_id, hess.d = 0.0001, hess.r = 8)


# check it has worked

(wk_check <- sapply(
  wk_group,
  simplify = FALSE, FUN = function(id) {
    assessment <- get(paste("biota", id, sep = "_"))
    ctsm_check_convergence(assessment)
  }))


# put each compenent back into assessment object

biota_assessment$assessment <- local({
  out <- biota_assessment$assessment
  for (group in c(wk_group, "Imposex")) {
    assessment <- paste0("biota_", group)  
    assessment <- get(assessment)
    out[names(assessment)] <- assessment
  }
  out
})


### tidy up ----

# make timeSeries factor levels more interpretable

biota_assessment$timeSeries <- biota_assessment$timeSeries %>% 
  rownames_to_column(".rownames") %>% 
  mutate(
    .matrix = ctsm_get_info("matrix", matrix, "name"), 
    .matrix = str_to_sentence(.matrix),
    level6name = case_when(
      level6element %in% "matrix" ~ .matrix,
      level6element %in% "sex" ~ recode(level6name, "F" = "Female", "M" = "Male"),
      TRUE ~ level6name
    ),
    .matrix = NULL,
    level6element = recode(
      level6element, 
      matrix = "Tissue", 
      sex = "Sex", 
      "METOA" = "Chemical analysis"
    )
  ) %>% 
  column_to_rownames(".rownames")


# saveRDS(biota_assessment, file.path("RData", "biota assessment.rds"))


## water ----

### main runs ----

water_assessment <- ctsm.assessment.setup(
  water_timeSeries, 
  AC = "EQS", 
  recent.trend = 20
)


library(parallel)
library(pbapply)

wk.cores <- detectCores()
wk.cluster <- makeCluster(wk.cores - 1)

clusterExport(
  wk.cluster, 
  c("water_assessment", "negTwiceLogLik",
    unique(c(
      objects(pattern = "ctsm*"), 
      objects(pattern = "get*"), 
      objects(pattern = "info*")
    ))
  )
)

clusterEvalQ(wk.cluster, {
  library(lme4)
  library(tidyverse)
})  

water_assessment$assessment <- ctsm.assessment(
  water_assessment, 
  clusterID = wk.cluster
)

stopCluster(wk.cluster)


### check convergence ----

(wk_check <- ctsm_check_convergence(water_assessment$assessment))

# refit a couple of time series

# "Poland_CZP CD Yes" - fixed effects on bounds
(wk_id <- wk_check[1])
water_assessment$assessment[wk_id] <- 
  ctsm.assessment(water_assessment, seriesID = wk_id, fixed_bound = 20)

# "Poland_HZP CD Yes" - fixed effects on bounds
(wk_id <- wk_check[2])
water_assessment$assessment[wk_id] <- 
  ctsm.assessment(water_assessment, seriesID = wk_id, fixed_bound = 20)

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

biota_web <- biota_assessment
biota_web$info$AC <- c("BAC", "EAC", "EQS", "MPC")

biota_web <- ctsm_web_initialise(
  biota_web,
  classColour = list(
    below = c(
      "BAC" = "green", 
      "EAC" = "green", 
      "EQS" = "green", 
      "MPC" = "green"
    ),
    above = c(
      "BAC" = "red", 
      "EAC" = "red", 
      "EQS" = "red", 
      "MPC" = "red"
    ), 
    none = "black"
  ),
  determinandGroups = webGroups)

# saveRDS(biota_web, file.path("RData", "biota web.rds"))


sediment_web <- ctsm_web_initialise(
  sediment_assessment,
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  determinandGroups = webGroups
)

# saveRDS(sediment_web, file.path("RData", "sediment web.rds"))


water_web <- ctsm_web_initialise(
  water_assessment, 
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  determinandGroups = webGroups
)

# saveRDS(water_web, file.path("RData", "water web.rds"))


## write tables ----

ctsm_summary_table(
  assessments = list(
    Biota = biota_web, 
    Sediment = sediment_web, 
    Water = water_web
  ),
  determinandGroups = webGroups,
  path = file.path("output", "example_HELCOM")
)

