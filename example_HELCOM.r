# Introduction ----

# This vignette shows how to do an assessment (mostly) following the approach 
# taken in HELCOM HOLAS3. 

# The data were extracted from the ICES data base using the XHAT facilities on
# the [ICES webservice](https://dome.ices.dk/api/swagger/index.html). The data
# were extracted on 28 August 2023 and were filtered using is_helcom_area = TRUE
# anad maxYear = 2020. The data were subsequently reduced in size to make them
# more manageable for this example.

# There are two exceptions to the HOLAS3 approach. First, imposex data are not
# assessed here - these have an additional level of complexity that will be
# explained in a subsequent vignette (not written yet). Second, the method for
# dealing with 'initial' data, unique to HELCOM, has not been implemented (this
# is not yet available in harsat).

# We'll begin with contaminants in water which are the simplest data to assess.
# We'll then move on to sediment and biota which have more features to consider.


# Water assessment ----

rm(list = objects())

devtools::load_all()


# First, we use read_data to read in the contaminant data, the station
# dictionary, and two important reference tables: the determinand and thresholds
# reference tables. There are several things to note about the function call:
# * purpose = "HELCOM" means that the assessment configuration is set to
# mirror that of the HOLAS3 assessment; you can change the assessment
# configuration in many ways using the control argument, but that is not
# explained here
# * data_dir identifies the directory storing the contaminant data (water.txt)
# and station dictionary (stations.txt); using file.path prevents any
# difficulties using forward or backward slashes when writing file paths
# * info_dir identifies the directory storing the reference tables. These must
# be named determinand.csv and thresholds_water.csv. The determinand table
# identifies which determinands are going to be assessed. If there are no
# thresholds, then there doesn't need to be a file called thresholds_water.csv.
# * you don't have to specify the extraction date, but it can help to keep track 
# of things

# As well as reading in the contaminant data and station dictionary, the
# function allocates each record in the contaminant data to a station in the
# station dictionary. The process for doing this is quite complicated (and it
# can take a few minutes to run) and we don't go into details year.

# Finally, a message is printed out saying that the argument max_year is set to
# 2020. (This is as it should be since we set maxYear to be 2020 in the data
# extraction.) But an important consequence is that a contaminant time series
# will only be assessed if it has some data in the period 2015 to 2020 (i.e. in
# the last six monitoring years).

water_data <- read_data(
  compartment = "water", 
  purpose = "HELCOM",                               
  contaminants = "water.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_HELCOM"),
  info_dir = "information", 
  extraction = "2023/08/23"
)  


# We next simplify the data so that they are in the correct format for running
# the assessment. This also involves deleting some data that do not meet the
# conditions for the assessment.

# Notice that a message appears talking about 'oddities'. Both tidy_data and
# create_timeseries (the next function call) do a lot of checking of the data
# and strange values are written to the oddities folder for you to have a look
# at (in the hope that, if there are errors, they will get corrected and
# resubmitted to the ICES database). It turns out there are no strange value at
# this stage, but there are in the step that follows.

water_data <- tidy_data(water_data)


# We now do some more data cleaning and then group the data into time series.
# Each time series consists of measurements of a single determinand at a single
# monitoring station. 

# The determinands.control is an important argument to consider when running an
# assessment. Here, PFOS is one of the determinands to be assessed. However,
# PFOS can also be submitted as N-PFOS and BR-PFOS, its linear and branched
# components. The argument below tells the code to sum records of N-PFOS and
# BR-PFOS from the same sample and relabel them as PFOS. There are more
# complicated examples of the use of determinands.control in the sediment and
# biota examples below

water_timeseries <- create_timeseries(
  water_data,
  determinands.control = list(
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum")
  )
)


# If you want to see a list of all the time series, then you can run code along
# the following lines:

get_timeseries(water_timeseries) |> head(10) 


# At last it is time to run the assessment. You need to specify which thresholds
# to use, otherwise the code will ignore with any of them. For water there is
# only the EQS, but you still need to specify it. Look at the thresholds
# reference table to see what is available if you are unsure. Note that AC
# stands for Assessment Criteria which is what thresholds are often called. The
# parallel argument tells the code to use parallel processing. This usually
# speeds things up considerably. The assessment took about 1.5 minutes to run on
# my laptop.

water_assessment <- run_assessment(
  water_timeseries, 
  AC = "EQS", 
  parallel = TRUE
)

# We now need to check whether there were any convergence issues. Lack of
# convergence often occurs because there are errors in the data (e.g. reported
# in incorrect units) or because there are large outliers, so the best thing to
# do is first check your data. However, convergence can also be difficult if
# there are a lot of less-than measurements in the time series. Describing how
# to tweak the control arguments to get convergence is beyond the scope of this
# vignette (need to write another vignette to discuss this). Fortunately, there
# were no convergence issues here.

check_assessment(water_assessment)


# It is time to look at the results! The vignette for external data shows how
# you can plot the data for each time series along with the model fitted to the
# data. The code below prints out a csv file giving summary information about
# the assessment of each time series. This includes:
# - meta-data such as the monitoring location and number of years of data for
# each time series
# - the fitted values in the last monitoring year with associated upper
# one-sided 95% confidence limits
# - the trend assessments (p-values and trend estimates)
# - the status assessments (if there any thresholds)
# - (optionally) a symbology summarising the trend (shape) and status (colour)
# of each time series

# This function is being actively developed and the function arguments are 
# likely to evolve, so we'll leave their explanation for the next release.

write_summary_table(
  water_assessment,
  determinandGroups = list(
    levels = c("Metals", "Organotins", "Organofluorines"), 
    labels = c("Metals", "Organotins", "Organofluorines")
  ),
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_HELCOM")
)


# Sediment assessment ----

# The sediment assessment is very similar, but has a few extra features related
# to normalisation (to account for differences in grain size) which are describe
# below

sediment_data <- read_data(
  compartment = "sediment", 
  purpose = "HELCOM",                               
  contaminants = "sediment.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_HELCOM"),
  info_dir = "information", 
  extraction = "2023/08/23"
)  

sediment_data <- tidy_data(sediment_data)


# The create_timeseries call for sediment differs from the call for water in two
# ways. First, determinands.control identifies two groups of determinands that
# need to be summed. Second, the arguments normalise and normalise.control
# specify how the normalisation for grain size should be carried out. There are
# default functions for normalisation that will work in many cases. However, the
# process for HELCOM is more complicated (because unlike other metals, copper is
# normalised to organic carbon) so a customised function
# normalise_sediment_HELCOM is provided. The argument normalise.control
# specifies that metals (apart from copper) will be normalised to 5% aluminium
# and copper and organics will be normalised to 5% organic carbon. The normalise
# functions need a bit of work, so expect them to change.

sediment_timeseries <- create_timeseries(
  sediment_data,
  determinands.control = list(
    SBDE6 = list(
      det = c("BDE28", "BDE47", "BDE99", "BD100", "BD153", "BD154"), 
      action = "sum"
    ),
    HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum")
  ),
  normalise = normalise_sediment_HELCOM,
  normalise.control = list(
    metals = list(method = "pivot", normaliser = "AL", value = 5), 
    copper = list(method = "hybrid", normaliser = "CORG", value = 5),
    organics = list(method = "simple", normaliser = "CORG", value = 5) 
  )
)


# Now run the assessment. Again there is only one threshold, the EQS.  This only
# takes about a minute to run on my laptop.

sediment_assessment <- run_assessment(
  sediment_timeseries, 
  AC = "EQS",
  parallel = TRUE
)

# Everything has converged tis time.

check_assessment(sediment_assessment)


# Finally, we can plot individual time series assessments (see vignett for
# external data) or print out the summary table

write_summary_table(
  sediment_assessment,
  determinandGroups = webGroups <- list(
    levels = c("Metals", "Organotins", "PAH_parent", "PBDEs", "Organobromines"),  
    labels = c(
      "Metals", "Organotins", "Polycyclic aromatic hydrocarbons",  
      "Organobromines", "Organobromines" 
    )
  ),
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_HELCOM")
)




# Biota assessment ----

# STOP HERE!!!


biota_data <- read_data(
  compartment = "biota", 
  purpose = "HELCOM",                               
  contaminants = "biota_data.csv", 
  stations = "station_dictionary.csv", 
  QA = "quality_assurance.csv",
  data_dir = file.path("data", "example_HELCOM"),
  data_format = "ICES_old",
  info_dir = "information",
  extraction = "2022/10/06",
  max_year = 2021L,
  control = list(
    region = list(id = c("HELCOM_subbasin", "HELCOM_L3", "HELCOM_L4"))
  )
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





# Prepare data for next stage ----

# gets correct variable and streamlines some of the data files

biota_data <- tidy_data(biota_data)



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


biota_timeseries <- create_timeseries(
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

biota_timeseries$data <- mutate(
  biota_timeseries$data, 
  matrix = if_else(
    year <= 2013 & matrix == "MU&EP", 
    "MU", 
    matrix
  )
)

# resolve Polish metoa changes

biota_timeseries$data <- left_join(
  biota_timeseries$data, 
  biota_data$stations[c("station_code", "country")],
  by = "station_code"
)


biota_timeseries$data <- mutate(
  biota_timeseries$data, 
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

biota_timeseries$data$country <- NULL





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

