# A very simple example (20s)

rm(list = objects())

devtools::load_all()

water_data <- read_data(
  compartment = "water", 
  purpose = "HELCOM",                               
  contaminants = "water.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_simple"),
  info_dir = "information", 
  extraction = "2023/08/23"
)  

water_data <- tidy_data(water_data)

water_timeseries <- create_timeseries(
  water_data,
  determinands.control = list(
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum")
  )
)

water_assessment <- run_assessment(
  water_timeseries, 
  AC = "EQS"
)

check_assessment(water_assessment)

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
  output_dir = file.path("output", "example_simple")
)

