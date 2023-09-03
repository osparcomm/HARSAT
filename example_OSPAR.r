# Intro ----

# Commentrary to be added

rm(list = objects())

devtools::load_all()


water_data <- read_data(
  compartment = "water", 
  purpose = "OSPAR",                               
  contaminants = "water.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_OSPAR"),
  info_dir = file.path("information", "OSPAR_2022"), 
  extraction = "2023/08/23"
)  

water_data <- tidy_data(water_data)

water_timeseries <- create_timeseries(
  water_data,
  determinands.control = list(
    CHR = list(det = "CHRTR", action = "replace"),
    BBKF = list(det = c("BBF", "BKF", "BBJF", "BBJKF"), action = "bespoke"),
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum"),
    CB138 = list(det = c("CB138+163"), action = "replace"),
    HCEPX = list(det = c("HCEPC", "HCEPT"), action = "sum"), 
    HCH = list(det = c("HCHA", "HCHB", "HCHG"), action = "sum")
  )
)


water_assessment <- run_assessment(
  water_timeseries, 
  AC = "EQS", 
  parallel = TRUE
)

check_assessment(water_assessment)


write_summary_table(
  water_assessment,
  determinandGroups = list(
    levels = c(
      "Metals", "Organotins", "PAH_parent",  "Organofluorines", 
      "Chlorobiphenyls", "Organochlorines", "Pesticides"
    ),  
    labels = c(
      "Metals", "Organotins", "PAH parent compounds", "Organofluorines", 
      "Polychlorinated biphenyls", "Organochlorines (other)", "Pesticides"
    )
  ),
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_OSPAR")
)
