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
  classColour = list(
    below = c("EQS" = "green"), 
    above = c("EQS" = "red"), 
    none = "black"
  ),
  collapse_AC = list(EAC = "EQS"),
  output_dir = file.path("output", "example_OSPAR")
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


info_TEQ <- c(
  "CB77" = 0.0001, "CB81" = 0.0003, "CB105" = 0.00003, "CB118" = 0.00003, 
  "CB126" = 0.1, "CB156" = 0.00003, "CB157" = 0.00003, "CB167" = 0.00003, 
  "CB169" = 0.03, "CDD1N" = 1, "CDD4X" = 0.1, "CDD6P" = 0.01, "CDD6X" = 0.1, 
  "CDD9X" = 0.1, "CDDO" = 0.0003, "CDF2N" = 0.3, "CDF2T" = 0.1, "CDF4X" = 0.1, 
  "CDF6P" = 0.01, "CDF6X" = 0.1, "CDF9P" = 0.01,
  "CDF9X" = 0.1, "CDFO" = 0.00003, "CDFP2" = 0.03, "CDFX1" = 0.1, "TCDD" = 1
)


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
    TEQDFP = list(det = names(info_TEQ), action = "bespoke"),
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
  classColour = list(
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
  ),
  collapse_AC = list(BAC = "BAC", EAC = c("EAC", "ERL", "EQS", "FEQG")),
  output_dir = file.path("output", "example_OSPAR")
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

info_TEQ <- c(
  "CB77" = 0.0001, "CB81" = 0.0003, "CB105" = 0.00003, "CB118" = 0.00003, 
  "CB126" = 0.1, "CB156" = 0.00003, "CB157" = 0.00003, "CB167" = 0.00003, 
  "CB169" = 0.03, "CDD1N" = 1, "CDD4X" = 0.1, "CDD6P" = 0.01, "CDD6X" = 0.1, 
  "CDD9X" = 0.1, "CDDO" = 0.0003, "CDF2N" = 0.3, "CDF2T" = 0.1, "CDF4X" = 0.1, 
  "CDF6P" = 0.01, "CDF6X" = 0.1, "CDF9P" = 0.01,
  "CDF9X" = 0.1, "CDFO" = 0.00003, "CDFP2" = 0.03, "CDFX1" = 0.1, "TCDD" = 1
)

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
    TEQDFP = list(det = names(info_TEQ), action = "bespoke"),
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
  AC = c("BAC", "NRC", "EAC", "LRC", "QSsp", "MPC", "QShh"),
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



write_summary_table(
  biota_assessment,
  determinandGroups = list(
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
  ),
  classColour = list(
    below = c(
      "BAC" = "blue",
      "NRC" = "blue",
      "EAC" = "green", 
      "LRC" = "green", 
      "QSsp" = "green", 
      "MPC" = "green",
      "QShh" = "green"
    ),
    above = c(
      "BAC" = "orange", 
      "NRC" = "orange", 
      "EAC" = "red", 
      "LRC" = "red", 
      "QSsp" = "red", 
      "MPC" = "red",
      "QShh" = "green"
    ),
    none = "black"
  ),
  collapse_AC = list(
    BAC = c("BAC", "NRC"),
    EAC = c("EAC", "LRC", "QSsp"), 
    HQS = c("MPC", "QShh")
  ),
  output_dir = file.path("output", "example_OSPAR")
)

