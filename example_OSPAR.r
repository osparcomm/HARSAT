# Intro ----

# Commentrary to be added

rm(list = objects())

devtools::load_all()


# Water ----

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
  symbology = list(
    colour = list(
      below = c("EQS" = "green"), 
      above = c("EQS" = "red"), 
      none = "black"
    )
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_OSPAR")
)


report_assessment(
  water_assessment, 
  subset = determinand %in% "ZN",
  output_dir = file.path("output", "reports")
)



# Sediment ----

sediment_data <- read_data(
  compartment = "sediment", 
  purpose = "OSPAR",                               
  contaminants = "sediment.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_OSPAR"),
  info_dir = file.path("information", "OSPAR_2022"), 
  extraction = "2023/08/23"
)  

sediment_data <- tidy_data(sediment_data)


sediment_timeseries <- create_timeseries(
  sediment_data,
  determinands.control = list(
    CHR = list(det = "CHRTR", action = "replace"),
    BBKF = list(det = c("BBF", "BKF", "BBJF", "BBJKF"), action = "bespoke"),
    NAPC1 = list(det = c("NAP1M", "NAP2M"), action = "sum"),
    BD154 = list(det = "PBB153+BD154", action = "replace"),
    HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum"),
    CB138 = list(det = c("CB138+163"), action = "replace"),
    CB156 = list(det = c("CB156+172"), action = "replace"),
    TEQDFP = list(
      det = names(info_TEF$DFP_CEMP), 
      action = "sum", 
      weights = info_TEF$DFP_CEMP
    ),
    HCEPX = list(det = c("HCEPC", "HCEPT"), action = "sum")
  ),
  normalise = normalise_sediment_OSPAR,
  normalise.control = list(
    metals = list(method = "pivot", normaliser = "AL"), 
    organics = list(method = "simple", normaliser = "CORG", value = 2.5),
    exclude = expression(ospar_subregion %in% c("Iberian Coast", "Gulf of Cadiz"))
  )
)


sediment_assessment <- run_assessment(
  sediment_timeseries, 
  AC = c("BAC", "EAC", "EQS", "ERL", "FEQG"),
  parallel = TRUE
)

# 02 00s

check_assessment(sediment_assessment)


write_summary_table(
  sediment_assessment,
  determinandGroups = list(
    levels = c(
      "Metals", "Organotins", "PAH_parent", "PAH_alkylated",  
      "PBDEs", "Organobromines", "Chlorobiphenyls", "Dioxins", 
      "Organochlorines"
    ),
    labels = c(
      "Metals", "Organotins", "PAH parent compounds", "PAH alkylated compounds", 
      "Polybrominated diphenyl ethers", "Organobromines (other)", 
      "Polychlorinated biphenyls", "Dioxins", "Organochlorines (other)"
    )
  ),
  symbology = list(
    colour = list(
      below = c(
        "BAC" = "blue", 
        "ERL" = "green", 
        "EAC" = "green", 
        "EQS" = "green", 
        "FEQG" = "green"
      ),
      above = c(
        "BAC" = "orange", 
        "ERL" = "red", 
        "EAC" = "red", 
        "EQS" = "red", 
        "FEQG" = "red"
      ),
      none = "black"
    )
  ),
  collapse_AC = list(BAC = "BAC", EAC = c("EAC", "ERL", "EQS", "FEQG")),
  output_dir = file.path("output", "example_OSPAR")
)


report_assessment(
  sediment_assessment, 
  subset = series == "10744 BD153 SEDTOT",
  output_dir = file.path("output", "reports")
)




# Biota ----

biota_data <- read_data(
  compartment = "biota", 
  purpose = "OSPAR",                               
  contaminants = "biota.txt", 
  stations = "stations.txt", 
  data_dir = file.path("data", "example_OSPAR"),
  info_dir = file.path("information", "OSPAR_2022"), 
  extraction = "2023/08/23"
)  

biota_data <- tidy_data(biota_data)


biota_timeseries <- create_timeseries(
  biota_data,
  determinands.control = list(
    CHR = list(det = "CHRTR", action = "replace"),
    BBKF = list(det = c("BBF", "BKF", "BBJF", "BBJKF"), action = "bespoke"),
    NAPC1 = list(det = c("NAP1M", "NAP2M"), action = "sum"),
    BD154 = list(det = "PBB153+BD154", action = "replace"),
    SBDE6 = list(
      det = c("BDE28", "BDE47", "BDE99", "BD100", "BD153", "BD154"), 
      action = "sum"
    ),
    HBCD = list(det = c("HBCDA", "HBCDB", "HBCDG"), action = "sum"),
    PFOS = list(det = c("N-PFOS", "BR-PFOS"), action = "sum"),
    PFHXS = list(det = c("N-PFHXS", "BR-PFHXS"), action = "sum"),
    CB138 = list(det = c("CB138+163"), action = "replace"),
    CB156 = list(det = c("CB156+172"), action = "replace"),
    SCB6 = list(
      det = c("CB28", "CB52", "CB101", "CB138", "CB153", "CB180"), 
      action = "sum"
    ),
    SCB7 = list(
      det = c("CB28", "CB52", "CB101", "CB118", "CB138", "CB153", "CB180"), 
      action = "sum"
    ),
    TEQDFP = list(
      det = names(info_TEF$DFP_CEMP), 
      action = "sum", 
      weights = info_TEF$DFP_CEMP
    ),
    HCEPX = list(det = c("HCEPC", "HCEPT"), action = "sum"),
    HCH = list(det = c("HCHA", "HCHB", "HCHG"), action = "sum"),
    "LIPIDWT%" = list(det = c("EXLIP%", "FATWT%"), action = "bespoke")
  ), 
  get_basis = get_basis_biota_OSPAR
)


# biota_assessment <- run_assessment(
#   biota_timeseries, 
#   AC = c("BAC", "EAC", "EQS", "HQS"),
#   parallel = TRUE
# )


wk_metals <- 
  c("AG", "AS", "CD", "CO", "CR", "CU", "HG", "NI", "PB", "SE", "SN", "ZN") 

biota_assessment <- run_assessment(
  biota_timeseries, 
  subset = determinand %in% wk_metals,
  AC = c("BAC", "NRC", "EAC", "FEQG", "LRC", "QSsp", "MPC", "QShh"),
  get_AC_fn = get_AC_biota_OSPAR,
  parallel = TRUE
)

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = !determinand %in% wk_metals, 
  parallel = TRUE
)


check_assessment(biota_assessment)

biota_assessment <- update_assessment(
  biota_assessment, 
  subset = series == "5031 BBKF Mytilus edulis SB",
  hess.d = 0.0001, hess.r = 8
)

check_assessment(biota_assessment)


# environmental summary

wk_groups <- list(
  levels = c(
    "Metals", "Organotins", 
    "PAH_parent", "PAH_alkylated", "Metabolites", 
    "PBDEs", "Organobromines", 
    "Organofluorines", 
    "Chlorobiphenyls", "Dioxins", "Organochlorines",
    "Effects"
  ),  
  labels = c(
    "Metals", "Organotins", 
    "PAH parent compounds", "PAH alkylated compounds", "PAH metabolites", 
    "Polybrominated diphenyl ethers", "Organobromines (other)", 
    "Organofluorines", 
    "Polychlorinated biphenyls", "Dioxins", "Organochlorines (other)",
    "Biological effects (other)"
  )
)

write_summary_table(
  biota_assessment,
  determinandGroups = wk_groups,
  symbology = list(
    colour = list(
      below = c(
        "BAC" = "blue",
        "NRC" = "blue",
        "EAC" = "green", 
        "FEQG" = "green",
        "LRC" = "green", 
        "QSsp" = "green"
      ),
      above = c(
        "BAC" = "orange", 
        "NRC" = "orange", 
        "EAC" = "red", 
        "FEQG" = "red",
        "LRC" = "red", 
        "QSsp" = "red"
      ),
      none = "black"
    )
  ),
  collapse_AC = list(
    BAC = c("BAC", "NRC"),
    EAC = c("EAC", "FEQG", "LRC", "QSsp"), 
    HQS = c("MPC", "QShh")
  ),
  output_file = "biota_summary_env.csv",
  output_dir = file.path("output", "example_OSPAR")
)


# health summary

write_summary_table(
  biota_assessment,
  determinandGroups = wk_groups,
  symbology = list(
    colour = list(
      below = c(
        "MPC" = "green", 
        "QShh" = "green"
      ),
      above = c(
        "MPC" = "red", 
        "QShh" = "red"
      ),
      none = "black"
    )
  ),
  collapse_AC = list(
    BAC = c("BAC", "NRC"),
    EAC = c("EAC", "FEQG", "LRC", "QSsp"), 
    HQS = c("MPC", "QShh")
  ),
  output_file = "biota_summary_health.csv",
  output_dir = file.path("output", "example_OSPAR")
)



report_assessment(
  biota_assessment, 
  subset = determinand %in% "ZN",
  output_dir = file.path("output", "reports")
)
