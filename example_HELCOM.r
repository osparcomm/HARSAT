# Introduction ----

# Water assessment ----

rm(list = objects())

devtools::load_all()

water_data <- read_data(
  compartment = "water", 
  purpose = "HELCOM",                               
  contaminants = "water.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_HELCOM"),
  info_dir = "./information/HELCOM_2023", 
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
  AC = "EQS", 
  parallel = TRUE
)

check_assessment(water_assessment)

write_summary_table(
  water_assessment,
  symbology = "default",
  symbology_control = list(
    colour = list(EQS = list(below = "green", above = "red"))
  ),
  output_dir = file.path("output", "example_HELCOM")
)


# Sediment assessment ----

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
    organics = list(method = "simple", normaliser = "CORG", value = 5),
    normalisers = list(method = "none")  
  )
)


sediment_assessment <- run_assessment(
  sediment_timeseries, 
  AC = "EQS",
  parallel = TRUE
)

check_assessment(sediment_assessment)


write_summary_table(
  sediment_assessment,
  determinandGroups = webGroups <- list(
    levels = c("Metals", "Organotins", "PAH_parent", "PBDEs", "Organobromines", 
               "Organic_constituents"),  
    labels = c(
      "Metals", "Organotins", "Polycyclic aromatic hydrocarbons",  
      "Organobromines", "Organobromines" , "Organic constituents"
    )
  ),
  symbology = "default",
  symbology_control = list(
    colour = list(EQS = list(below = "green", above = "red"))
  ),
  output_dir = file.path("output", "example_HELCOM")
)


# Biota assessment ----

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

biota_timeseries <- create_timeseries(
  biota_data,
  #get_basis = get_basis_most_common,
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


biota_assessment <- run_assessment(
  biota_timeseries, 
  AC = c("BAC", "EAC", "EQS", "MPC")
)

check_assessment(biota_assessment)

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = series == "2299 PYR1OH Limanda limanda BI HPLC-FD",
  hess.d = 0.0001, hess.r = 8
)

check_assessment(biota_assessment)


write_summary_table(
  biota_assessment,
  threshold_groups = list(EQS = c("BAC", "EAC", "EQS", "MPC")),
  symbology = "default",
  symbology_control = list(
    colour = list(
      EQS = list(below = "green", above = "red")
    )
  ),
  determinandGroups = list(
    levels = c(
      "Metals", "PAH_parent", "Metabolites", "PBDEs", "Organobromines", 
      "Organofluorines", "Chlorobiphenyls", "Dioxins", "Biological"
    ),  
    labels = c(
      "Metals", "PAH compounds and metabolites", "PAH compounds and metabolites",
      "Organobromines", "Organobromines", "Organofluorines", 
      "PCBs and dioxins", "PCBs and dioxins", "Biological"
    )
  ),
  output_dir = file.path("output", "example_HELCOM")
)

