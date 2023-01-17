# function to check if species numerical values are within pre-defined range 

values_range_check_species <- function(species_data, min_value, max_value) {
  library(tidyverse)
  
  # select subset of columns
  numerical_columns <- select(species_data, ends_with("WT."))
  
  # convert to numerical type
  numerical_columns <- transform(numerical_columns, 
                                 MU_DRYWT.   = as.numeric(MU_DRYWT.), 
                                 LI_DRYWT.   = as.numeric(LI_DRYWT.),
                                 SB_DRYWT.   = as.numeric(SB_DRYWT.), 
                                 TM_DRYWT.   = as.numeric(TM_DRYWT.), 
                                 EH_DRYWT.   = as.numeric(EH_DRYWT.), 
                                 HA_DRYWT.   = as.numeric(HA_DRYWT.), 
                                 FE_DRYWT.   = as.numeric(FE_DRYWT.), 
                                 BL_DRYWT.   = as.numeric(BL_DRYWT.), 
                                 MU_LIPIDWT. = as.numeric(MU_LIPIDWT.),
                                 LI_LIPIDWT. = as.numeric(LI_LIPIDWT.),
                                 SB_LIPIDWT. = as.numeric(SB_LIPIDWT.), 
                                 TM_LIPIDWT. = as.numeric(TM_LIPIDWT.), 
                                 EH_LIPIDWT. = as.numeric(EH_LIPIDWT.), 
                                 HA_LIPIDWT. = as.numeric(HA_LIPIDWT.), 
                                 FE_LIPIDWT. = as.numeric(FE_LIPIDWT.), 
                                 BL_LIPIDWT. = as.numeric(BL_LIPIDWT.)
  )
  
  # check the condition if values in pre-defined range
  condition <- (numerical_columns >= min_value) & (numerical_columns <= max_value)
  
  # convert NA to TRUE to only have boolean values
  condition[is.na(condition)] <- TRUE
  
  # check if all values are TRUE
  check <- all(condition == TRUE)
  
  if(!check)
    print("Warning: not all values in species reference table are within range [min_value, max_value] !!!")
  else
    print("OK, all values in species reference table are within range [min_value, max_value]")
}


# Information extractor function ----

ctsm_get_info <- function(
  info_type, 
  input, 
  output, 
  compartment = NULL, 
  na_action = c("fail", "input_ok", "output_ok", "ok"), 
  sep = ".") {
  
  # information_functions.R
  # gets data from standard reference tables 
  
  # na_action: 
  #   fail doesn't allow any missing values
  #   input.ok allows missing values in input, but all non-missing values must be 
  #     recognised, and must have output
  #   output.ok requires all input values to be present, but allows missing values in 
  #     output (for e.g. dryweight by species)
  #   ok allows missing values everywhere
  
  na_action <- match.arg(na_action)
  
  # construct input variables and check that all input elements are recognised in 
  # information files
  
  info_file <- get(paste("info", info_type, sep = "."))
  input <- as.character(input)
  
  # check whether input is a combination of values - sometimes used when e.g. there are 
  # two methods used in the extraction of a chemical 
  
  split_input <- any(grepl("~", na.omit(input)))
  if (split_input) {
    input2 <- strsplit(input, "~", fixed = TRUE)
  }
  
  
  # check for failure due to missing values
  
  wk <- unique(if(split_input) unlist(input2) else input)
  ok <- switch(
    na_action, 
    fail = wk %in% rownames(info_file),
    input_ok = wk %in% c(rownames(info_file), NA),
    output_ok = wk %in% rownames(info_file),
    TRUE
  )
  
  if (any(!ok)) {
    stop('Not found in ', info_type, ' information file: ', 
         paste(wk[!ok], collapse = ", "))
  }
  
  
  # construct output variables and check that all information is present 
  
  if (!is.null(compartment)) {
    output <- paste(compartment, output, sep = sep)
  }
  
  if (!(output %in% names(info_file))) { 
    stop('Incorrect specification of output variable in function ctsm_get_info')
  }
  
  
  # check that if the input has multiple values (i.e. has had to be split) each element
  # gives the same output - then simplify input to just one of the relevant values
  
  if (split_input) {
    ok <- sapply(input2, function(i) length(unique(info_file[i, output])) == 1)
    if (any(!ok)) {
      stop('Incompatible data found in ', info_type, ' information file: ', 
           paste(input[!ok], collapse = ", "))
    }
    input <- sapply(input2, "[", 1)
  }
  
  out <- info_file[input, output]
  
  ok <- switch(
    na_action,
    fail = !is.na(out),
    input_ok = is.na(input) | (!is.na(input) & !is.na(out)),
    TRUE
  )
  
  if (any(!ok)) { 
    stop ('Missing values for following in ', info_type, ' information file: ', 
          paste(unique(input[!ok]), collapse = ", "))
  }

  out
}  




# Species information ----

info.path <- sub("functions", "information", function_path)
info.file <- function(file.name) file.path(info.path, file.name)

info.species <- read.csv(
    info.file(info_species_file_id), 
    row.names = "submitted_species", 
    na.strings = "", 
    check.names = FALSE
)


# Uncertainty estimates ----

info.uncertainty <- read.csv(
  info.file(info_uncertainty_file_id), 
  row.names = "determinand", 
  na.strings = ""
)


# Determinand information and functions ----

info.determinand <- read.csv(
  info.file(info_determinand_infile),
  na.strings = "", 
  colClasses = c(
    determinand = "character",
    common_name = "character",
    pargroup = "character", 
    biota_group = "character",
    sediment_group = "character",
    water_group = "character",
    biota_assess = "logical",
    sediment_assess = "logical",
    water_assess = "logical",
    biota_unit = "character",
    sediment_unit = "character", 
    water_unit = "character",       
    biota_auxiliary = "character", 
    sediment_auxiliary = "character", 
    water_auxiliary = "character", 
    distribution = "character", 
    good_status = "character"
  )
)

info.determinand <- dplyr::mutate(
  info.determinand, 
  common_name = ifelse(
    is.na(.data$common_name), 
    .data$determinand, 
    .data$common_name
  )
)

info.determinand <- tibble::column_to_rownames(info.determinand, "determinand")



# extractor functions

ctsm_get_determinands <- function(compartment = c("biota", "sediment", "water")) {
  
  # information_functions.R
  # gets determinands to be assessed from determinand reference table
  
  compartment <- match.arg(compartment)
  
  assess_id <- paste0(compartment, "_assess")
  ok <- info.determinand[[assess_id]]
  
  row.names(info.determinand)[ok]
}  


ctsm_get_auxiliary <- function(
  determinands, 
  compartment = c("biota", "sediment", "water")) {
  
  # information_functions.R
  # gets required auxiliary variables for determinands
  
  # in case determinands is a factor
  determinands <- as.character(determinands)
  
  determinands <- unique(determinands)
  
  # auxiliary_id <- paste0(compartment , "_auxiliary")
  # auxiliary <- info.determinand[determinands, auxiliary_id]
  
  auxiliary <- ctsm_get_info(
    "determinand", 
    determinands, 
    "auxiliary", 
    compartment, 
    na_action = "output_ok", 
    sep = "_"
  )
  
  auxiliary <- strsplit(auxiliary, ", ")
  auxiliary <- unlist(auxiliary)
  
  unique(c(na.omit(auxiliary)))
}


# check all auxiliary variables in info.determinand are recognised as 
# determinands in their own right

lapply(c("biota", "sediment", "water"), function(compartment) {
  
  determinands <- row.names(info.determinand)
  auxiliary <- ctsm_get_auxiliary(determinands, compartment)

  ok <- auxiliary %in% row.names(info.determinand)
  if(!all(ok)) {
    stop(
      'Not found in determinand information file: ', 
      paste(auxiliary[!ok], collapse = ", ")
    )
  }
})



# Toxic EQuivalents for WHO_DFP (health)

info_TEQ <- c(
  "CB77" = 0.0001, "CB81" = 0.0003, "CB105" = 0.00003, "CB118" = 0.00003, "CB126" = 0.1, 
  "CB156" = 0.00003, "CB157" = 0.00003, "CB167" = 0.00003, "CB169" = 0.03, 
  "CDD1N" = 1, "CDD4X" = 0.1, "CDD6P" = 0.01, "CDD6X" = 0.1, "CDD9X" = 0.1, "CDDO" = 0.0003,
  "CDF2N" = 0.3, "CDF2T" = 0.1, "CDF4X" = 0.1, "CDF6P" = 0.01, "CDF6X" = 0.1, "CDF9P" = 0.01,
  "CDF9X" = 0.1, "CDFO" = 0.00003, "CDFP2" = 0.03, "CDFX1" = 0.1, "TCDD" = 1
)





# Assessment criteria ----

read.assessment.criteria <- function(infile)  {

  sediment <- read.csv(
    info.file(infile$sediment), 
    na.strings = ""
  )
  sediment <- within(sediment, country[is.na(country)] <- "")

  biota <- read.csv(
    info.file(infile$biota), 
    na.strings = ""
  )
  wk <- strsplit(biota$sub.family, ",", fixed = TRUE)
  n <- sapply(wk, length)
  biota <- biota[rep(1:nrow(biota), times = n),]
  biota$sub.family <- unlist(wk)

  water <- read.csv(
    info.file(infile$water), 
    na.strings = ""
  )
  
  list(sediment = sediment, biota = biota, water = water)
}

info.assessment.criteria <- read.assessment.criteria(info_AC_infile)
rm(read.assessment.criteria)

# gets Assessment Criteria
# determinand is a vector of length n
# info can be either a data frame or a list of appropriate supporting variables
# if a list, then each element can be either a scalar (replicated to length n), or a vector of length n
# AC is a character vector of assessment concentration types

get.AC <- function(compartment, determinand, info, AC) {

  # check elements of info are of correct length
 
  stopifnot(sapply(info, length) %in% c(1, length(determinand)))
  

  # remove duplicate determinand information - need to fix this better
  
  info$determinand <- NULL
  
    
  # turn info into a dataframe if necessary

  data <- cbind(determinand, as.data.frame(info))
  
                
  # split by determinand groupings
  
  group <- ctsm_get_info(
    "determinand", data$determinand, "group", compartment, sep = "_"
  )
  
  data <- split(data, group, drop = TRUE)

  
  # get assessment concentrations
   
  out <- lapply(names(data), function(i) {
    do.call(paste("get.AC", compartment, i, sep = "."), list(data = data[[i]], AC = AC))
  }) 
  
  unsplit(out, group, drop = TRUE)
}


get.AC.biota.contaminant <- function(data, AC, export_cf = FALSE) {    
  
  AC_data <- info.assessment.criteria$biota
  stopifnot(AC %in% names(AC_data))
  
  AC_data <- mutate_if(AC_data, is.factor, as.character)
  
  data <- data %>% 
    rownames_to_column("rownames") %>% 
    mutate_if(is.factor, as.character) %>%
    mutate(
      species_group = ctsm_get_info("species", species, "species_group"),
      sub.family = ctsm_get_info("species", species, "species_subgroup")
    ) %>% 
    mutate_at(c("species_group", "sub.family"), as.character)
  
  lipid_info <- info.species %>% 
    rownames_to_column("species") %>% 
    select(.data$species, contains("LIPIDWT%")) %>% 
    gather(key = "matrix", value = "lipid_wt", contains("LIPIDWT%"), na.rm = TRUE) %>% 
    separate(matrix, c("matrix", NA), sep = "_") 
  
  drywt_info <- info.species %>% 
    rownames_to_column("species") %>% 
    select(.data$species, contains("DRYWT%")) %>% 
    gather(key = "matrix", value = "dry_wt", contains("DRYWT%"), na.rm = TRUE) %>% 
    separate(matrix, c("matrix", NA), sep = "_") 
  
  data <- left_join(data, lipid_info, by = c("species", "matrix"))
  data <- left_join(data, drywt_info, by = c("species", "matrix"))
  
  data <- data[c("rownames", "determinand", "sub.family", "basis", "dry_wt", "lipid_wt")]
  
  data <- left_join(data, AC_data, by = c("determinand", "sub.family"))
  
  out <- sapply(AC, simplify = FALSE, FUN = function(i) {
    basis_AC <- paste("basis", i, sep = ".")
    convert.basis(data[[i]], from = data[[basis_AC]], to = data$basis, data$dry_wt, "", data$lipid_wt, "")
  })
  
  out <- data.frame(out)
  
  if (export_cf) {
    out <- bind_cols(out, data[c("dry_wt", "lipid_wt")])
  }
  
  rownames(out) <- data$rownames
  out
}                           


if (info_AC_type == "OSPAR") {

  get.AC.biota.PAH_parent <- get.AC.biota.contaminant
  get.AC.biota.PAH_alkylated <- get.AC.biota.contaminant
  get.AC.biota.PBDEs <- get.AC.biota.contaminant
  get.AC.biota.Organotins <- get.AC.biota.contaminant
  get.AC.biota.Organobromines <- get.AC.biota.contaminant
  
  get.AC.biota.Metals <- function(data, AC, lipid_high = 3) {
    
    out <- get.AC.biota.contaminant(data, AC, export_cf = TRUE)
    
    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "EQS.OSPAR", "HQS") %in% names(AC)
    )
    
    out <- bind_cols(out, data)
    
    out <- out %>%
      rownames_to_column() %>%
      mutate(
        species_group = ctsm_get_info("species", .data$species, "species_group"),
        sub.family = ctsm_get_info("species", .data$species, "species_subgroup")
      )
    
    
    # mercury
    # only BAC is for mussels
    # no MPC (HQS) for fish liver
    # mammal liver 16000 (BAC), 64000 (EAC) ww
    # mammal hair 6100 (BAC), 24400 (EAC) dw
    # bird egg homongenate 110 (BAC), 470 (EAC) ww
    # bird liver 1400 (BAC) 7300 (EAC) ww
    # bird feather 1580 (BAC) 7920 (EAC) dw
    # bird blood 200 (BAC) 1000 (EAC) ww
    
    id <- out$determinand %in% "HG"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        
        BAC = case_when(
          .data$species_group %in% "Fish"                            ~ NA_real_,
          .data$sub.family %in% "Oyster"                      ~ NA_real_,
          .data$species_group %in% "Mammal" & .data$matrix %in% "LI" ~
            convert.basis(16000, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Mammal" & .data$matrix %in% "HA" ~
            convert.basis(6100, "D", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "EH"   ~
            convert.basis(110, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "LI"   ~
            convert.basis(1400, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "FE"   ~
            convert.basis(1580, "D", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "BL"   ~
            convert.basis(200, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          TRUE                                                ~ .data$BAC
        ),
        
        EQS.OSPAR = case_when(
          .data$species_group %in% "Mammal" & .data$matrix %in% "LI" ~
            convert.basis(64000, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Mammal" & .data$matrix %in% "HA" ~
            convert.basis(24400, "D", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "EH"   ~
            convert.basis(470, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "LI"   ~
            convert.basis(7300, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "FE"   ~
            convert.basis(7920, "D", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          .data$species_group %in% "Bird" & .data$matrix %in% "BL"   ~
            convert.basis(1000, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
          TRUE                                                ~ .data$EQS.OSPAR
        ),
        
        HQS = if_else(.data$species_group %in% "Fish" & .data$matrix %in% "LI", NA_real_, .data$HQS)
      )
    }
    
    
    # cadmium and lead
    # adjust HQS (MPC) for fish muscle
    # BACs in fish only apply to high lipid tissue
    
    out <- mutate(
      out,
      
      HQS = if_else(
        .data$determinand %in% "CD" & .data$species_group %in% "Fish" & .data$matrix %in% "MU",
        50,
        .data$HQS
      ),
      
      HQS = if_else(
        .data$determinand %in% "PB" & .data$species_group %in% "Fish" & .data$matrix %in% "MU",
        300,
        .data$HQS
      ),
      
      BAC = if_else(
        .data$determinand %in% c("CD", "PB") & .data$species_group %in% "Fish" &
          (is.na(.data$lipid_wt) | .data$lipid_wt < lipid_high),
        NA_real_,
        .data$BAC
      )
    )
    
    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))
    
    out
  }
  
  
  get.AC.biota.Chlorobiphenyls <- function(data, AC, lipid_high = 3) {
    
    out <- get.AC.biota.contaminant(data, AC, export_cf = TRUE)
    
    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "EAC", "HQS") %in% names(AC)
    )
    
    out <- bind_cols(out, data)
    
    out <- out %>%
      rownames_to_column() %>%
      mutate(
        species_group = ctsm_get_info("species", .data$species, "species_group"),
        sub.family = ctsm_get_info("species", .data$species, "species_subgroup")
      )
    
    # BACs in fish only apply to high lipid tissue
    # add in MPC (HQS) of 200 ww for fish liver
    
    out <- mutate(
      out,
      
      BAC = if_else(
        .data$species_group %in% "Fish" & (is.na(.data$lipid_wt) | .data$lipid_wt < lipid_high),
        NA_real_,
        .data$BAC
      ),
      
      HQS = if_else(
        .data$species_group %in% "Fish" & .data$matrix %in% "LI" & .data$determinand %in% "SCB6",
        convert.basis(200, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .data$HQS
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "SCB7",
        convert.basis(6.7, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .data$EAC
      )
    )
    
    
    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))
    
    out
  }
  
  
  get.AC.biota.Organochlorines <- function(data, AC, lipid_high = 3) {
    
    out <- get.AC.biota.contaminant(data, AC, export_cf = TRUE)
    
    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "EAC", "EQS.OSPAR", "HQS") %in% names(AC)
    )
    
    out <- bind_cols(out, data)
    
    out <- out %>%
      rownames_to_column() %>%
      mutate(
        species_group = ctsm_get_info("species", .data$species, "species_group"),
        species = as.character(species)
      )
    
    
    # BAC in fish only apply to high lipid tissue
    # EAC for birds fro HCB and HCH
    
    out <- mutate(
      out,
      
      BAC = if_else(
        .data$species_group %in% "Fish" & (is.na(.data$lipid_wt) | .data$lipid_wt < lipid_high),
        NA_real_,
        .data$BAC
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "HCB",
        convert.basis(2.0, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .data$EAC
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "HCH",
        convert.basis(2.0, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .data$EAC
      )
    )
    
    
    # HCHG HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "HCHG" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_LIPIDWT%"],
        .dry_mu = info.species[.data$species, "MU_DRYWT%"],
        HQS = convert.basis(61, "W", "L", .dry_mu, "", .lipid_mu, ""),
        HQS = convert.basis(.data$HQS, "L", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .lipid_mu = NULL,
        .dry_mu = NULL
      )
    }
    
    
    # HCB HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "HCB" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_LIPIDWT%"],
        .dry_mu = info.species[.data$species, "MU_DRYWT%"],
        HQS = convert.basis(10, "W", "L", .dry_mu, "", .lipid_mu, ""),
        HQS = convert.basis(.data$HQS, "L", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .lipid_mu = NULL,
        .dry_mu = NULL
      )
    }
    
    
    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))
    
    out
  }
  
  
  get.AC.biota.Organofluorines <- function(data, AC) {
    
    out <- get.AC.biota.contaminant(data, AC, export_cf = TRUE)
    
    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "EQS.OPAR", "HQS") %in% names(AC)
    )
    
    out <- bind_cols(out, data)
    
    out <- out %>%
      rownames_to_column() %>%
      mutate(species_group = ctsm_get_info("species", .data$species, "species_group"))
    
    
    # fish liver - multiply by 5
    
    out <- mutate(
      out,
      .id <- .data$species_group %in% "Fish" & .data$matrix %in% "LI" &
        .data$determinand %in% "PFOS",
      EQS.OSPAR = .data$EQS.OSPAR * if_else(.id, 5, 1),
      HQS = .data$HQS * if_else(.id, 5, 1),
      .id = NULL
    )
    
    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))
    
    out
  }
  
  
  get.AC.biota.PBDEs <- function(data, AC) {
    
    out <- get.AC.biota.contaminant(data, AC, export_cf = TRUE)
    
    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "FEQG", "HQS") %in% names(AC)
    )
    
    out <- bind_cols(out, data)
    
    out <- out %>%
      rownames_to_column() %>%
      mutate(
        species_group = ctsm_get_info("species", .data$species, "species_group"),
        species = as.character(species)
      )
    
    
    # HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "SBDE6" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_LIPIDWT%"],
        .dry_mu = info.species[.data$species, "MU_DRYWT%"],
        HQS = convert.basis(0.0085, "W", "L", .dry_mu, "", .lipid_mu, ""),
        HQS = convert.basis(.data$HQS, "L", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .lipid_mu = NULL,
        .dry_mu = NULL
      )
    }
    
    
    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))
    
    out
  }
  
  
  get.AC.biota.Dioxins <- function(data, AC) {
    
    out <- get.AC.biota.contaminant(data, AC, export_cf = TRUE)
    
    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "EQS.OPAR", "HQS") %in% names(AC)
    )
    
    out <- bind_cols(out, data)
    
    out <- out %>%
      rownames_to_column() %>%
      mutate(species_group = ctsm_get_info("species", .data$species, "species_group"))
    
    
    # add in HQS of 0.02 ww for fish liver
    
    out <- mutate(
      out,
      HQS = if_else(
        .data$species_group %in% "Fish" & .data$matrix %in% "LI" & .data$determinand %in% "TEQDFP",
        convert.basis(0.02, "W", .data$basis, .data$dry_wt, "", .data$lipid_wt, ""),
        .data$HQS
      )
    )
    
    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))
    
    out
  }
  
  
  get.AC.biota.Effects <- function(data, AC) {
    
    out <- as.data.frame(do.call("cbind", sapply(AC, function(i) rep(NA, nrow(data)), simplify = FALSE)))
    rownames(out) <- rownames(data)
    
    with(data, {
      
      if ("EROD" %in% data$determinand) {
        
        stopifnot("matrix" %in% names(data))
        
        id <- determinand %in% "EROD" & matrix %in% "LIMIC"
        if (any(id) & "BAC" %in% AC) {
          out$BAC[id & species %in% "Limanda limanda"] <- 680
          out$BAC[id & species %in% "Gadus morhua"] <- 145
          out$BAC[id & species %in% "Pleuronectes platessa"] <- 255
          out$BAC[id & species %in% "Lepidorhombus boscii"] <- 13
          out$BAC[id & species %in% "Callionymus lyra"] <- 202
        }
        
        if ("LIS9" %in% data$matrix) {
          stopifnot("sex" %in% names(data))
          
          id <- determinand %in% "EROD" & matrix %in% "LIS9"
          
          if ("BAC" %in% AC) {
            out$BAC[id & species %in% "Limanda limanda" & sex %in% "F"] <- 178
            out$BAC[id & species %in% "Limanda limanda" & sex %in% "M"] <- 147
            out$BAC[id & species %in% "Platichthys flesus" & sex %in% "M"] <- 24
            out$BAC[id & species %in% "Pleuronectes platessa" & sex %in% "M"] <- 9.5
            out$BAC[id & species %in% "Mullus barbatus" & sex %in% "M"] <- 208
          }
        }
      }
      
      if ("SFG" %in% data$determinand) {
        id <- ctsm_get_info("species", species, "sub.family") %in% "Mussel" &
          determinand %in% "SFG"
        if ("BAC" %in% AC) out$BAC[id] <- 25
        if ("EAC" %in% AC) out$EAC[id] <- 15
      }
      
      if ("SURVT" %in% data$determinand) {
        id <- ctsm_get_info("species", species, "sub.family") %in% "Mussel" &
          determinand %in% "SURVT"
        if ("BAC" %in% AC) out$BAC[id] <- 10
        if ("EAC" %in% AC) out$EAC[id] <- 5
      }
      
      if ("NRR" %in% data$determinand) {
        id <- determinand %in% "NRR"
        if ("BAC" %in% AC) out$BAC[id] <- 120
        if ("EAC" %in% AC) out$EAC[id] <- 50
      }
      
      if ("LP" %in% data$determinand) {
        id <- determinand %in% "LP"
        if ("BAC" %in% AC) out$BAC[id] <- 20
        if ("EAC" %in% AC) out$EAC[id] <- 10
      }
      
      if ("MNC" %in% data$determinand) {
        
        if (any(ctsm_get_info("species", species, "sub.family") %in% "Mussel" &
                determinand %in% "MNC"))
          stop("AC not coded for MNC in mussels")
        
        id <- determinand %in% "MNC"
        if ("BAC" %in% AC) {
          out$BAC[id & species %in% "Platichthys flesus"] <- 0.3
          out$BAC[id & species %in% "Limanda limanda"] <- 0.5
          out$BAC[id & species %in% "Zoarces viviparus"] <- 0.4
          out$BAC[id & species %in% "Gadus morhua"] <- 0.4
          out$BAC[id & species %in% "Mullus barbatus"] <- 0.3
        }
      }
      
      if ("%DNATAIL" %in% data$determinand) {
        id <- determinand %in% "%DNATAIL"
        if (any(id) & "BAC" %in% AC) {
          out$BAC[id & species %in% "Mytilus edulis"] <- 10
          out$BAC[id & species %in% "Gadus morhua"] <- 5
          out$BAC[id & species %in% "Limanda limanda"] <- 5
        }
      }
      
      out
    })
  }
  
  
  get.AC.biota.Metabolites <- function(data, AC) {
    
    out <- as.data.frame(do.call("cbind", sapply(AC, function(i) rep(NA, nrow(data)), simplify = FALSE)))
    rownames(out) <- rownames(data)
    
    with(data, {
      
      stopifnot("metoa" %in% names(data))
      
      if ("BAC" %in% AC) {
        id <- species %in% "Limanda limanda"
        out$BAC[id & determinand %in% "PYR1OH" & metoa %in% "HPLC-FD"] <- 16
        out$BAC[id & determinand %in% "PA1OH" & metoa %in% "HPLC-FD"] <- 3.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 0.15
        
        id <- species %in% "Gadus morhua"
        out$BAC[id & determinand %in% "PYR1OH" & metoa %in% "HPLC-FD"] <- 21
        out$BAC[id & determinand %in% "PA1OH" & metoa %in% "HPLC-FD"] <- 2.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 1.1
        
        id <- species %in% "Platichthys flesus"
        out$BAC[id & determinand %in% "PYR1OH" & metoa %in% "HPLC-FD"] <- 16
        out$BAC[id & determinand %in% "PA1OH" & metoa %in% "HPLC-FD"] <- 3.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 1.3
        
        id <- species %in% "Melanogrammus aeglefinus"
        out$BAC[id & determinand %in% "PYR1OH" & metoa %in% "HPLC-FD"] <- 13
        out$BAC[id & determinand %in% "PA1OH" & metoa %in% "HPLC-FD"] <- 0.8
        out$BAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 1.9
      }
      
      if ("EAC" %in% AC) {
        id <- species %in% "Limanda limanda"
        out$EAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 22
        
        id <- species %in% "Gadus morhua"
        out$EAC[id & determinand %in% "PYR1OH" & metoa %in% "GC-MS"] <- 483
        out$EAC[id & determinand %in% "PA1OH" & metoa %in% "GC-MS"] <- 528
        out$EAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 35
        
        id <- species %in% "Platichthys flesus"
        out$EAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 29
        
        id <- species %in% "Melanogrammus aeglefinus"
        out$EAC[id & determinand %in% "PYR1OHEQ" & metoa %in% "FLM-SS"] <- 35
      }
      
      out
    })
  }
  
  
  get.AC.biota.Imposex <- function(data, AC) {
    out <- as.data.frame(do.call("cbind", sapply(AC, function(i) rep(NA, nrow(data)), simplify = FALSE)))
    rownames(out) <- rownames(data)
    
    with(data, {

      if ("BAC" %in% AC)
      {
        out$BAC[determinand %in% "VDS" & species %in% "Nucella lapillus"] <- 0.3
        out$BAC[determinand %in% "VDS" & species %in% "Neptunea antiqua"] <- 0.3
      }
      
      if ("EAC" %in% AC)
      {
        out$EAC[determinand %in% "VDS" & species %in% "Nucella lapillus"] <- 2.0
        out$EAC[determinand %in% "VDS" & species %in% "Neptunea antiqua"] <- 2.0
        out$EAC[determinand %in% "VDS" & species %in% "Tritia nitida / reticulata"] <- 0.3
        out$EAC[determinand %in% "VDS" & species %in% "Buccinum undatum"] <- 0.3
      }
      
      out
    })
  }

}


if (info_AC_type == "HELCOM") {
  
  get.AC.biota.PAH_parent <- get.AC.biota.contaminant
  get.AC.biota.PAH_alkylated <- get.AC.biota.contaminant
  get.AC.biota.PBDEs <- get.AC.biota.contaminant
  get.AC.biota.Organotins <- get.AC.biota.contaminant
  get.AC.biota.Organobromines <- get.AC.biota.contaminant
  get.AC.biota.Chlorobiphenyls <- get.AC.biota.contaminant
  get.AC.biota.Organochlorines<- get.AC.biota.contaminant
  get.AC.biota.Dioxins <- get.AC.biota.contaminant
  
  get.AC.biota.Metals <- function(data, AC) {

    out <- get.AC.biota.contaminant(data, AC)

    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      ! c("BAC", "EQS") %in% names(AC)
    )

    out <- bind_cols(out, data)

    out <- rownames_to_column(out)


    # lead

    out <- mutate(
      out,
      BAC = if_else(
        .data$determinand %in% "PB" & .data$matrix %in% "MU",
        NA_real_,
        as.double(.data$BAC)
      ),
      EQS = if_else(
        .data$determinand %in% "PB" & .data$matrix %in% "LI",
        NA_real_,
        .data$EQS
      )
    )

    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))

    out
  }


  get.AC.biota.Organofluorines <- function(data, AC) {

    out <- get.AC.biota.contaminant(data, AC)

    stopifnot(
      length(intersect(names(data), names(out))) == 0,
      !c("EQS") %in% names(AC)
    )

    out <- bind_cols(out, data)

    out <- rownames_to_column(out)

    # fish liver - multiply by 5

    out <- mutate(
      out,
      EQS = if_else(
        .data$matrix %in% "LI" & .data$determinand %in% "PFOS",
        .data$EQS * 17.9,
        .data$EQS
      )
    )

    out <- out %>%
      column_to_rownames() %>%
      select(all_of(AC))

    out
  }

  
  get.AC.biota.Metabolites <- function(data, AC) {
    
    out <- as.data.frame(do.call("cbind", sapply(AC, function(i) rep(NA, nrow(data)), simplify = FALSE)))
    rownames(out) <- rownames(data)
    
    with(data, {
      
      stopifnot("metoa" %in% names(data))
      
      if ("EAC" %in% AC) {
        out$EAC[determinand %in% "PYR1OH" & metoa %in% "HPLC-FD"] <- 483
      }
      
      out
    })
  }
  
  
  get.AC.biota.Imposex <- function(data, AC) {
    out <- as.data.frame(do.call("cbind", sapply(AC, function(i) rep(NA, nrow(data)), simplify = FALSE)))
    rownames(out) <- rownames(data)
    
    with(data, {
      
      if ("EAC" %in% AC)
      {
        out$EAC[determinand %in% "VDS" & species %in% "Nucella lapillus"] <- 2.0
        out$EAC[determinand %in% "VDS" & species %in% "Neptunea antiqua"] <- 2.0
        out$EAC[determinand %in% "VDS" & species %in% "Tritia nitida / reticulata"] <- 0.3
        out$EAC[determinand %in% "VDS" & species %in% "Buccinum undatum"] <- 0.3
        out$EAC[determinand %in% "VDS" & species %in% "Peringia ulvae"] <- 0.1
        out$EAC[determinand %in% "INTS" & species %in% "Littorina littorea"] <- 0.3
      }
      
      out
    })
  }

}



get.AC.sediment.contaminant <- function(data, AC, export_all = FALSE) {   

  AC.data <- info.assessment.criteria$sediment
  stopifnot(AC %in% names(AC.data))
  
  out <- with(data, {
    determinand <- as.character(determinand)
    basis <- as.character(basis)
  
    if ("country" %in% names(data)) {
      country <- as.character(country)
      country[country != "Spain"] <- ""
    }
    else
      country <- rep("", length(determinand))
  
    wk.data <- data.frame(determinand, basis, country, id = 1:length(determinand))
  
    out <- merge(wk.data, AC.data, all.x = T)
    out[order(out$id),]
  })
  
  rownames(out) <- rownames(data)
  
  if (export_all) {
    return(out)
  } 
  
  out[AC]
}                           


get.AC.sediment.Metals <- function(data, AC) {
  
  out <- get.AC.sediment.contaminant(data, AC, export_all = TRUE)

  # manual adjustment when ERLs are less than BACs 
  
  # version 2_64 and before: use ERL for AS and NI for Spain, but only BAC for 
  # other countries
  
  # if (all(c("ERL", "BAC") %in% AC)) {
  #   out <- within(out, {
  #     id <- !is.na(ERL) & !is.na(BAC)
  #     ERL[id & ERL < BAC] <- NA
  #     BAC[id & ERL == BAC] <- NA
  #     rm(id)
  #   })
  # }
  
  # version 2_65 onwards: don't use ERL at all for AS and NI and 
  # ensure ERL and not BAC is used for CR

  out <- tibble::rownames_to_column(out)
  
  if ("BAC" %in% AC) {
    out <- dplyr::mutate(
      out, 
      BAC = dplyr::if_else(.data$determinand %in% "CR", NA_real_, .data$BAC)
    )
  }  
  
  if ("ERL" %in% AC) {
    out <- dplyr::mutate(
      out,
      ERL = dplyr::if_else(.data$determinand %in% "AS", NA_real_, .data$ERL),
      ERL = dplyr::if_else(.data$determinand %in% "NI", NA_real_, .data$ERL)
    )
  }
    
  out <- tibble::column_to_rownames(out)
  
  out[AC]
}                           

get.AC.sediment.PAH_parent <- get.AC.sediment.contaminant
get.AC.sediment.PAH_alkylated <- get.AC.sediment.contaminant
get.AC.sediment.Chlorobiphenyls <- get.AC.sediment.contaminant
get.AC.sediment.PBDEs <- get.AC.sediment.contaminant
get.AC.sediment.Organobromines <- get.AC.sediment.contaminant
get.AC.sediment.Organofluorines <- get.AC.sediment.contaminant
get.AC.sediment.Organochlorines <- get.AC.sediment.contaminant
get.AC.sediment.Organotins <- get.AC.sediment.contaminant
get.AC.sediment.Dioxins <- get.AC.sediment.contaminant


get.AC.water.contaminant <- function(data, AC) {   

  AC.data <- info.assessment.criteria$water
  stopifnot(AC %in% names(AC.data))
  
  out <- data[c("determinand", "basis")]

  out <- within(out, {
    determinand <- as.character(determinand)
    basis <- as.character(basis)
    id = 1:length(determinand)
  })
  
  out <- merge(out, AC.data, all.x = TRUE)
  out <- out[order(out$id),]

  rownames(out) <- rownames(data)
  out[AC]
}                           

get.AC.water.Metals <- get.AC.water.contaminant
get.AC.water.PAH_parent <- get.AC.water.contaminant
get.AC.water.PAH_alkylated <- get.AC.water.contaminant
get.AC.water.Chlorobiphenyls <- get.AC.water.contaminant
get.AC.water.PBDEs <- get.AC.water.contaminant
get.AC.water.Organobromines <- get.AC.water.contaminant
get.AC.water.Organochlorines <- get.AC.water.contaminant 
get.AC.water.Organofluorines <- get.AC.water.contaminant
get.AC.water.Organotins <- get.AC.water.contaminant
get.AC.water.Dioxins <- get.AC.water.contaminant
get.AC.water.Pesticides <- get.AC.water.contaminant


# Unit conversion ----

convert_units <- function(conc, from, to) {

  # information_functions.R
  # converts units; e.g. mg/kg to ug/kg
  
  # can supply non-standard units (e.g. for biological effects) provided that 
  #   from and to for these rows are identical (in which case no attempt is made 
  #   to convert
  # missing values are permitted provided consistent in from and to
  
  # error checking
    
  from = as.character(from)
  to = as.character(to)
  
  if (!(length(from) %in% c(1, length(conc)))) {
    stop("lengths of input data inconsistent: compare lengths of from and conc")
  }
  
  if (!(length(to) %in% c(1, length(conc)))) {
    stop("lengths of input data inconsistent: compare lengths of to and conc")
  }
  
  
  # set up working data frame
  
  data <- data.frame(conc, from, to)
  
  not_ok <- is.na(data$from) & !is.na(data$to)
  if (any(not_ok)) stop("cannot convert units since input unit is missing")

  not_ok <- !is.na(data$from) & is.na(data$to)
  if (any(not_ok)) stop("cannot convert units since output unit is missing")
  
  
  # identify the rows where a conversion is required

  convert <- !is.na(data$from) & (data$from != data$to)
  
  if (!any(convert)) return(data$conc)


  # setup output structure and then apply conversion to restricted rows 
  # where conversion required
  
  out <- data$conc
  
  
  # get unit type and magnitude for each unit
  
  to <- convert_units_workup(data$to[convert])
  from <- convert_units_workup(data$from[convert])
  
  
  # check that units in from and to are compatible
  
  if (!identical(to$group, from$group)) {
    stop('Cannot convert between units of different types')
  }
  
  
  # finally convert units
  
  out[convert] <- out[convert] * 10^(to$magnitude - from$magnitude)
  
  out
}

  
convert_units_workup <- function(units) {
  
  # information_functions.R
  # helper function for convert_units
  
  # valid unit measurements (for conversion) and their 'type'
  # can only convert between comparable types 
  # note the special case of percentage in the weight_weight category
    
  # could add in more units but those below cover all the ones needed to date
  
  unit_group <- list(
    "length" = c("km", "m", "cm", "mm"), 
    "weight" = c("kg", "g", "mg"), 
    "volume" = c("l", "ml"), 
    "weight_volume" = c(
      "kg/l", "g/l", "mg/l", "ug/l", "ng/l", "pg/l", 
      "g/ml", "mg/ml", "ug/ml", "ng/ml", "pg/ml"
    ), 
    "weight_weight" = c(
      "kg/kg", "g/kg", "mg/kg", "ug/kg", "ng/kg", "pg/kg", 
      "g/g", "mg/g", "ug/g", "ng/g", "pg/g", 
      "mg/mg", "ug/ug", "ng/ng", "pg/pg", 
      "%"
    ),
    "TEQ_weight_weight" = c("TEQ ug/kg", "TEQ pg/g"),
    "mol_min_weight" = c(
      "umol/min/mg protein", "nmol/min/mg protein", "pmol/min/mg protein"
    )
  )
  
  
  # the magnitude of each unit 
  
  unit_magnitude <- list(
    "-3" = c("km", "kg"),
    "0" = c(
      "m", "g", "l", "kg/l", "g/ml", 
      "kg/kg", "g/g", "mg/mg", "ug/ug", "ng/ng", "pg/pg"
    ),  
    "2" = c("%", "cm"), 
    "3" = c("mm", "mg", "ml", "g/l", "mg/ml", "g/kg", "mg/g"), 
    "6" = c("mg/l", "ug/ml", "mg/kg", "ug/g", "umol/min/mg protein"), 
    "9" = c("ug/l", "ng/ml", "ug/kg", "ng/g", "TEQ ug/kg", "nmol/min/mg protein"), 
    "12" = c("ng/l", "pg/ml", "ng/kg", "pg/g", "TEQ pg/g", "pmol/min/mg protein"), 
    "15" = c("pg/l", "pg/kg")
  )
  
  
  # ensure the two lists have identical units
  
  check1 <- sort(unname(unlist(unit_group)))
  check2 <- sort(unname(unlist(unit_magnitude)))
  stopifnot(identical(check1, check2))
  
  
  # check all units in the data have a valid type and magnitude 
  
  unit_levels <- unlist(unit_group)
  
  ok <- units %in% unit_levels
  if (any(!ok)) {
    stop('Unrecognised units: ', paste(unique(units[!ok]), collapse = ", "))
  }
  
  
  # create data frame with the group and mangitude of each unit in the data
  
  out <- data.frame(units)
  
  out$group <- out$magnitude <- factor(out$units, levels = unit_levels)
  
  levels(out$group) <- unit_group
  out$group <- as.character(out$group)
  
  levels(out$magnitude) <- unit_magnitude
  out$magnitude <- as.numeric(as.character(out$magnitude))
  
  out        
}



# Basis and matrix information and basis conversion ----

convert.basis <- function(
  conc, from, to, drywt, drywt.qflag, lipidwt = NA, lipidwt.qflag = NA, exclude) {

  library(dplyr)
  

  # converts between wet, dry and lipid basis of measurement

  data <- data.frame(
    conc, from, to, drywt, drywt.qflag, lipidwt, lipidwt.qflag, 
    stringsAsFactors = FALSE
  )
  
  if (missing(exclude))
    exclude <- rep(FALSE, length(data$conc))

  stopifnot(
    data$from[!exclude] %in% c("W", "D", "L", NA), 
    data$to[!exclude] %in% c("W", "D", "L")
  )


  # check to see if any conversions needed

  ok <- exclude | (is.na(data$from) | as.character(data$from) == as.character(data$to))
  if (all(ok)) return (data$conc)


  # do conversion on subset of data that needs it
  # ensure columns of data are of correct type

  data <- mutate_at(data, c("drywt", "lipidwt"), as.numeric)
  data <- mutate_at(data, c("drywt.qflag", "lipidwt.qflag"), as.character)
  
  data$conc[!ok] <- with(data[!ok, ], {

    from_value <- case_when(
      from == "W" ~ 100, 
      from == "D" ~ drywt, 
      from == "L" ~ lipidwt
    )
    
    from_qflag <- case_when(
      from == "W" ~ "", 
      from == "D" ~ drywt.qflag, 
      from == "L" ~ lipidwt.qflag
    )

    to_value <- case_when(
      to == "W" ~ 100, 
      to == "D" ~ drywt,
      to == "L" ~ lipidwt
    )
    
    to_qflag <- case_when(
      to == "W" ~ "", 
      to == "D" ~ drywt.qflag,
      to == "L" ~ lipidwt.qflag
    )

    qflag_ok <- from_qflag %in% "" & to_qflag %in% ""
    
    conc <- case_when(
      ! qflag_ok ~ NA_real_,
      TRUE ~ conc * from_value / to_value
    )
    conc
  })

  data$conc
}


get_basis <- function(purpose, ...) {
  do.call(paste("get_basis", purpose, sep = "_"), list(...))
}


get_basis_OSPAR <- function(compartment, group, matrix, determinand, species, lipid_high = 3.0) {
  
  switch(
    compartment, 
    sediment = if (missing(group)) "D" else rep("D", length(group)),
    water = if (missing(group)) "W" else rep("W", length(group)),
    biota = { 
           
      # combine input variables and get species family
           
      out <- data.frame(species, matrix, determinand, group)
           
      out <- mutate(
        out,
        across(everything(), as.character),
        species_group = ctsm_get_info("species", .data$species, "species_group")
      )
           
      # get typical lipid content by species and matrix 
           
      lipid_info <- info.species %>% 
        rownames_to_column("species") %>% 
        select(.data$species, contains("LIPIDWT%")) %>% 
        gather(key = "matrix", value = "lipid_wt", contains("LIPIDWT%"), na.rm = TRUE) %>% 
        separate(matrix, c("matrix", NA), sep = "_") 
      
      out <- left_join(out, lipid_info, by = c("species", "matrix"))
           
      # default basis W
      # bivalves and gastropods - D
      # fish and crustacea:
      #   organobromines and organochlorines (except chlorinated paraffins) L
      # mammals (based on data submissions): 
      #   hair D 
      #   metals W
      #   organics L (note organofluorines submitted on L, with no associated
      #     lw, so can't convert)
      # birds: 
      #   Alle alle (BL, FE) D
      #   Rissa tridactyla (ER) D
      #   remaining data (other than EH) (BL, FE, LI, MU) W
      #   EH metals W (apart from Larus argentatus D)
      #   EH organofluorines W
      #   EH organics L (Cepphus grylle, Haematopus ostralegus, Sterna hirundo)
      #               W (Larus argentatus, Somateria mollissima)

      lw_group <- c("PBDEs", "Organobromines", "Chlorobiphenyls", "Dioxins", "Organochlorines")

      out <- mutate(
        out,
        .lw = .data$group %in% lw_group & !(.data$determinand %in% c("MCCP", "SCCP")),
        new.basis = case_when(
          .data$group %in% c("Imposex", "Effects", "Metabolites")       ~ NA_character_,
          .data$species_group %in% c("Bivalvia", "Gastropoda")                 ~ "D",
          .data$species_group %in% c("Fish", "Crustacea") & 
            .lw &
            .data$lipid_wt >= lipid_high                                ~ "L",
          .data$species_group %in% c("Fish", "Crustacea")                      ~ "W",
          .data$species_group %in% "Mammal" &
            .data$matrix %in% "HA"                                      ~ "D",
          .data$species_group %in% "Mammal" &
            .data$group %in% "Metals"                                   ~ "W",
          .data$species_group %in% "Mammal"                                    ~ "L",
          .data$species %in% c("Alle alle", "Rissa tridactyla")         ~ "D",
          .data$matrix %in% "EH" &
            .data$group %in% "Metals" & 
            .data$species %in% "Larus argentatus"                       ~ "D",
          .data$matrix %in% "EH" &
            .data$group %in% "Metals"                                   ~ "W",
          .data$matrix %in% "EH" &
            .data$group %in% "Organofluorines"                          ~ "W",
          .data$matrix %in% "EH" & 
            .data$species %in% c(
              "Cepphus grylle", "Haematopus ostralegus", "Sterna hirundo"
            )                                                           ~ "L",
          .data$matrix %in% "EH" & 
            .data$species %in% c(
              "Larus argentatus", "Somateria mollissima"
            )                                                           ~ "W",
          .data$species_group %in% "Bird"                                      ~ "W"
        )
      )
      
      out$new.basis
    }
  )
}


get_basis_CSEMP <- get_basis_OSPAR


get_basis_HELCOM <- function(compartment, group, matrix, determinand, species) {
  
  switch(
    compartment, 
    sediment = if (missing(group)) "D" else rep("D", length(group)),
    water = if (missing(group)) "W" else rep("W", length(group)),
    biota = { 
           
      # combine input variables and get species family
           
      out <- data.frame(species, matrix, determinand, group)
      
      out <- mutate(
        out,
        across(everything(), as.character),
        species_group = ctsm_get_info("species", .data$species, "species_group")
      )
 
      # define new basis
      
      out <- mutate(
        out, 
        new.basis = case_when(
          .data$group %in% c("Imposex", "Metabolites")       ~ NA_character_,
          .data$species_group %in% "Bivalvia"                       ~ "W",
          .data$species_group %in% "Fish" & 
            .data$group %in% c("Metals", "Organofluorines")  ~ "W",
          .data$species_group %in% "Fish"                           ~ "L"
        )
      ) 

      out$new.basis
    }
  )
}

get_basis_AMAP <- function(data, compartment) {

  # chooses basis which is most reported (regardless of when or whether 
  # auxiliary variables are also present to enable conversion)
    
  if (compartment != "biota")
    stop("not coded")
  
  id <- do.call(paste, data[c("station", "species", "matrix", "group")])
  
  out <- by(data, id, function(x) {
    out <- with(x, table(basis))
    x$new.basis <- names(out)[which.max(out)]
    x
  })
  
  out <- do.call(rbind, out)
  
  out <- within(out, new.basis <- factor(new.basis))
  
  out
}


info.matrix <- read.csv(info.file("matrix.csv"), row.names = "matrix", stringsAsFactors = FALSE)


# Regions ----

info.regions <- sapply(
  c("AMAP", "HELCOM", "OSPAR"), 
  function(x) {
    infile <-info.file(paste(x, "regions.csv"))
    row_names_id <- switch(
      x,
      OSPAR = "OSPAR_subregion",
      HELCOM = "HELCOM_L4",
      AMAP = NULL
    )
    if (file.exists(infile))
      read.csv(infile, row.names = row_names_id)
    else 
      NULL
  },
  simplify = FALSE
)



# Method of extraction and pivot values ----

info.methodExtraction <- read.csv(
  info.file("method of extraction.csv"), 
  row.names = "METCX",  
  na.strings = ""
)

info.pivotValues <- read.csv(info.file("pivot values.csv"), na.strings = "")


# Html ----

# something has changed with this - need to investigate
# but there are better replacement functions anyway

# info.html <- read.csv(info.file("HTMLtranslate.csv"))


# Imposex ----

info.imposex <- read.csv(info.file("imposex.csv"))

get.info.imposex <- function(species, determinand, choice = c("min_value", "max_value"), 
                             na.action = c("fail", "ok")) {
  
  choice <- match.arg(choice)
  na.action <- match.arg(na.action)

  info <- info.imposex[[choice]]
  names(info) <- do.call("paste", info.imposex[c("species", "determinand")])
  
  id <- paste(species, determinand)
  ok <- id %in% names(info)
  
  if (!all(ok)) {
    message <- paste0("Species determinand combinations not recognised: ", 
                      paste(unique(id[!ok]), collapse = ", "))
    if (na.action == "fail") stop(message) else warning(message)
  }
  
  out <- rep(NA, length(id))
  out[ok] <- info[id[ok]]  
  names(out) <- NULL
  out
}


rm(info.path, info.file)


# ICES RECO codes ----

# reads in data from csv files exported from ICES RECO
#
# generally based on read.csv, but sometimes csv export compromised e.g when commas are present in 
# chemical parameter descriptions


get_RECO <- function(code, path = "information") {
  
  library(tidyverse)
  library(lubridate)
  
  # check code argument and convert to upper case
  
  if(!is.character(code) | length(code) != 1L)
    stop("code must be a length 1 character")
  
  code <- toupper(code)
  
  
  # get list of available RECO files
  
  files <- list.files(path)
  
  ok <-substring(files, 1, 5) == "RECO_"
  
  if (!any(ok)) 
    stop("no RECO files found")
  
  files <- files[ok]
  
  
  # get relevant file
  
  file_codes <- files %>% 
    strsplit("_") %>% 
    sapply("[[", 2)
  
  ok <- file_codes %in% code

  if (!any(ok)) 
    stop("RECO file for this code not found")
  
  if (sum(ok) >= 2L)
    stop("multiple REcO files found for this code")
  
  infile <- files[ok]
  
  infile <- file.path(path, infile)
  

  # read in file - need to do something different for PARAM
  
  if (!code %in% "PARAM") {
    out <- infile %>%
      read_csv(
        col_types = cols(
          tblCodeID = col_integer(),
          Code = col_character(),
          Description = col_character(),
          tblCodeTypeID = col_integer(),
          CodeType = col_character(),
          Created = col_date(format = ""),
          Modified = col_date(format = ""),
          Deprecated = col_logical()
        )
      ) %>%
      as.data.frame()
    
    names(out)[2] <- code
    return(out)
  }


  # special case for PARAM required because there are columns in the Description column
  
  # read in each row as a single character string and split by commas
  
  data <- infile %>% 
    read_delim(
      delim = "\n",
      col_types = cols(
        `tblCodeID,Code,Description,tblCodeTypeID,CodeType,Created,Modified,Deprecated` = col_character()
      )
    ) %>% 
    as.data.frame() %>%
    set_names("X") %>% 
    pull(.data$X) %>% 
    strsplit(",")
  
  out <- data.frame(
    tblCodeID = sapply(data, function(x) {n <- length(x); x[1]}),
    Code = sapply(data, function(x) {n <- length(x); x[2]}),
    Description = sapply(data, function(x) {n <- length(x); paste(x[3:(n-5)], collapse = ",")}),
    tblCodeTypeID = sapply(data, function(x) {n <- length(x); x[n-4]}),
    CodeType = sapply(data, function(x) {n <- length(x); x[n-3]}),
    Created = sapply(data, function(x) {n <- length(x); x[n-2]}),
    Modified = sapply(data, function(x) {n <- length(x); x[n-1]}),
    Deprecated = sapply(data, function(x) {n <- length(x); x[n]}),
    stringsAsFactors = FALSE
  )

  out <- out %>% 
    mutate_at(c("tblCodeID", "tblCodeTypeID"), as.integer) %>% 
    mutate_at(c("Created", "Modified"), as_date) %>% 
    mutate_at(c("Deprecated"), as.logical)

  names(out)[2] <- code
  out
}


# Station utility functions ----

# get station code from station name
  
get_station_code <- function(name, country, stations) {
  
  # gets the station code corresponding to the station name and country from the 
  # station dictionary
  
  # only works for one country at a time
  
  stopifnot(length(country) == 1)
  n <- length(name)
  
  out <- data.frame(station = name, country = country) 
  out <- mutate(out, across(.fns = as.character))
    
  
  stations <- stations[c("station", "country", "code")]
  stations <- mutate(stations, across(.fns = as.character)) 
  
  out <- left_join(out, stations, by = c("station", "country"))
  
  stopifnot(!is.na(out), n == nrow(out))
  
  out$code
}
