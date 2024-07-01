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
  info_dir = "./information/HELCOM_2023", 
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
  symbology = list(
    colour = list(
      below = c("EQS" = "green"), 
      above = c("EQS" = "red"), 
      none = "black"
    )
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
  info_dir = "./information/HELCOM_2023", 
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

# Everything has converged this time.

check_assessment(sediment_assessment)


# Finally, we can plot individual time series assessments (see vignette for
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
  symbology = list(
    colour = list(
      below = c("EQS" = "green"), 
      above = c("EQS" = "red"), 
      none = "black"
    )
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_HELCOM")
)




# Biota assessment ----

# The main difference in the biota assessment is the inclusion of effects data.
# This example has some PAH metabolite time series, but all imposex data have
# been excluded to keep things relatively simple. Imposex assessments have an
# additional modelling stage and this will be described in another vignette (not
# yet available). 

# The first two stages are just as before

biota_data <- read_data(
  compartment = "biota", 
  purpose = "HELCOM",                               
  contaminants = "biota.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_HELCOM"),
  info_dir = "./information/HELCOM_2023", 
  extraction = "2023/08/23"
)

biota_data <- tidy_data(biota_data)


# The determinands.control argument does rather more here. There are four summed
# variables: PFOS, SBDE6, HBCD and SCB6. There is also one variable CB138+163
# which needs to be relabeled as (replaced by) CB138. For the purposes of the
# assessment, the contribution of CB163 is regarded as small, and CB138+163 is
# taken to be a good proxy for CB138. Note that the replacements must be done
# before the six PCBs are summed to give SCB6 in order for them to be included
# in the sum.

# There are also two 'bespoke' actions in determinands.control. These are
# customised functions that do more complicated (non-standard things). One of
# them computes the TEQDFP. The other deals with the three different ways in
# which lipid weight measurements can be submitted.

# Finally, normalise_biota_HELCOM is a customised function that determines which
# measurements are normalised to 5% lipid in a HELCOM assessment. Again, the 
# normalisation functions are under active development and might well change 
# before the first release.

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
    TEQDFP = list(
      det = names(info_TEF$DFP_HOLAS3),
      action = "sum",
      weights = info_TEF$DFP_HOLAS3
    ),
    "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
  ),
  normalise = normalise_biota_HELCOM,
  normalise.control = list(
    lipid = list(method = "simple", value = 5), 
    other = list(method = "none") 
  )
)


# The asssessment took about 3.5 minutes on my laptop

biota_assessment <- run_assessment(
  biota_timeseries, 
  AC = c("BAC", "EAC", "EQS", "MPC"),
  parallel = TRUE
)

check_assessment(biota_assessment)

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = series == "2299 PYR1OH Limanda limanda BI HPLC-FD",
  hess.d = 0.0001, hess.r = 8
)

check_assessment(biota_assessment)



# And that's it :)

write_summary_table(
  biota_assessment,
  determinandGroups = list(
    levels = c(
      "Metals", "PAH_parent", "Metabolites", "PBDEs", "Organobromines", 
      "Organofluorines", "Chlorobiphenyls", "Dioxins"
    ),  
    labels = c(
      "Metals", "PAH compounds and metabolites", "PAH compounds and metabolites",
      "Organobromines", "Organobromines", "Organofluorines", 
      "PCBs and dioxins", "PCBs and dioxins"
    )
  ),
  symbology = list(
    colour = list(
      below = c("BAC" = "green", "EAC" = "green", "EQS" = "green", "MPC" = "green"),
      above = c("BAC" = "red", "EAC" = "red", "EQS" = "red", "MPC" = "red"),
      none = "black"
    )
  ),
  collapse_AC = list(EAC = c("EAC", "EQS", "MPC")),
  output_dir = file.path("output", "example_HELCOM")
)

