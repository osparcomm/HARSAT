# Information extractor function ----

ctsm_get_info <- function(
  ref_table, 
  input, 
  output, 
  compartment = NULL, 
  na_action = c("fail", "input_ok", "output_ok", "ok"),
  info_type = NULL,
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

  # if not supplied, pick up info_type from function call - designed for use
  # internally; add trailing space for printing
  
  if (is.null(info_type)) {
    wk <- substitute(ref_table)
    wk <- as.character(wk)
    if (length(wk) == 3L && identical(wk[1:2], c("$", "info"))) {
      info_type <- paste0(wk[3], " ")
    } else {
      info_type <- ""
    }   
  } else {
    info_type <- paste0(info_type, " ")
  }

  
  # ensure input is character
    
  input <- as.character(input)
  
  
  # check whether input is a combination of values - sometimes used when e.g. 
  # there are two methods used in the extraction of a chemical 
  
  split_input <- any(grepl("~", na.omit(input)))
  if (split_input) {
    input2 <- strsplit(input, "~", fixed = TRUE)
  }
  
  
  # check for failure due to missing values in data or in reference table 
  
  wk <- if(split_input) unlist(input2) else input

  if (na_action %in% c("fail", "output_ok") && any(is.na(wk))) {
    stop(
      "Missing input values to ctsm_get_info with na.action set to ", 
      na_action
    )
  } 
  
  if (na_action != "ok") {
    ctsm_check_reference_table(na.omit(wk), ref_table, info_type)
  }
  

  # construct output variables and check that all information is present 
  
  if (!is.null(compartment)) {
    output <- paste(compartment, output, sep = sep)
  }
  
  if (!(output %in% names(ref_table))) { 
    stop('Incorrect specification of output variable in function ctsm_get_info')
  }
  
  
  # check that if the input has multiple values (i.e. has had to be split) each 
  # element gives the same output - then simplify input to just one of the 
  # relevant values
  
  if (split_input) {
    ok <- sapply(input2, function(i) {
      wk <- ref_table[i, output]
      length(unique(wk)) == 1L
    })
    if (any(!ok)) {
      stop(
        '\nIncompatible data found in the ', info_type, 'reference table:\n', 
        paste(input[!ok], collapse = ", "), 
        call. = FALSE
      )
    }
    input <- sapply(input2, "[", 1)
  }
  
  out <- ref_table[input, output]
  
  ok <- switch(
    na_action,
    fail = !is.na(out),
    input_ok = is.na(input) | (!is.na(input) & !is.na(out)),
    TRUE
  )
  
  if (any(!ok)) { 
    stop (
      '\nMissing output for the following input values in the', info_type, 
      'reference table:\n', 
      paste(unique(input[!ok]), collapse = ", "), 
      call. = FALSE
    )
  }

  out
}  


ctsm_check_reference_table <- function(x, ref_table, info_type = "") {

  # information_functions.R
  # checks whether x is in the reference table
  
  id <- unique(x)
  
  ok <- id %in% row.names(ref_table)
  
  if (!all(ok)) {
    id <- id[!ok]
    id <- sort(id)
    stop(
      '\nThe following values are not in the ', info_type, 'reference table.\n',
      "Please add them to the reference table or edit your data to continue.\n",
      paste(id, collapse = ", "), 
      call. = FALSE
    )
  }
  
  invisible()
}


# Species information ----

ctsm_read_species <- function(file, path = "information") {

  var_id <- c(
    "reference_species" = "character",
    "submitted_species" = "character",
    "common_name" = "character",
    "species_group" = "character", 
    "species_subgroup" = "character", 
    "assess" = "logical"
  )
  
  required <- names(var_id)
  

  # check required variables are present in data
  
  data <- read.csv(
    file.path(path, file), 
    strip.white = TRUE, 
    nrows = 1
  )
  
  ok <- required %in% names(data)
  
  if (!all(ok)) {
    id <- required[!ok]
    id <- sort(id)
    stop(
      "\nThe following variables are missing from ", file, ".\n", 
      "Please update the file to continue, noting that variable names are case ", 
      "sensitive.\n",
      "Variables: ", paste(id, collapse = ", "),
      .call = FALSE
    )
  }
  
  
  # add drywt and lipidwt conversion factors to var_id

  cf_id <- grep("_drywt|_lipidwt", names(data), value = TRUE)
  is_cf <- length(cf_id > 0)
  
  if (is_cf) {
    cf_class <- rep("numeric", length(cf_id))
    names(cf_class) <- cf_id
    var_id <- c(var_id, cf_class)
    ok <- c(required, cf_id) %in% names(data)
  }
      
  
  data <- read.csv(
    file.path(path, file), 
    na.strings = c("", "NULL"),
    strip.white = TRUE,
    colClasses = var_id[ok]
  )
  
  
  # placeholder to create missing (non-required) variables 


  # check no missing values in required variables
  
  ok <- apply(data[required], 2, function(x) !any(is.na(x)))
  
  if (!all(ok)) {
    id <- required[!ok]
    id <- sort(id)
    stop(
      "\nThe following variables have missing values in ", file, ".\n", 
      "Please update the file to continue.\n",
      "Variables: ", paste(id, collapse = ", "),
      .call = FALSE
    )
  }
  
  
  # check no duplicate values in submitted_species

  ok <- !duplicated(data$submitted_species)
  
  if (!all(ok)) {
    id <- data$submitted_species[!ok]
    stop(
      "\nDuplicate 'submitted_species' not allowed in ", file, ".\n",
      "Please update the file to continue.\n",
      "Duplicated: ", paste(id, collapse = ", "),
      call. = FALSE
    )
  }
  
  
  # check recognised_species is a subset of submitted_species
  
  ok <- data$reference_species %in% data$submitted_species
  
  if (!all(ok)) {
    id <- data$recognised_species[!ok]
    stop(
      "\nSome 'reference_species' are not in 'submitted_species' in ", file, ".\n",
      "Please update the file to continue.\n",
      "Missing reference_species: ", paste(id, collapse = ", "),
      call. = FALSE
    )
  }
  
  
  # check species_group values are recognised
  
  group_id <- c(
    "Bird", "Bivalve", "Crustacean", "Fish", "Gastropod", "Macrophyte", 
    "Mammal", "Other"
  )
  
  ok <- data$species_group %in% group_id

  if (!all(ok)) {
    id <- data$species_group[!ok]
    id <- sort(id)
    stop(
      "\nUnrecongised 'species_group' values in ", file, ".\n",
      "Please update the file to continue or contact the HARSAT development ",
      "team to update\nthe list of recognised values:\n", 
      paste(group_id, collapse = ", "), "\n",
      "Unrecognised values: ", paste(id, collapse = ", "), 
      call. = FALSE
    )
  }
  

  # check no forward or backward slashes
  # will have to come back to backward slashes as these are hard to code for
  
  bad <- apply(data, 2, function(x) any(grepl("/", x, fixed = TRUE)))

  if (any(bad)) {
    id <- names(data)[bad]
    id <- sort(id)
    stop(
      "\nVariables in ", file, " have forward slashes which are not allowed.\n", 
      "Please update the file to continue.\n",
      "Variables: ", paste(id, collapse = ", "), 
      call. = FALSE
    )
  }
  
  # check conversion factors are between 0 and 100  
  
  if (is_cf) {
    values_range_check_species(data, 0, 100)
  }  
  
  
  # tidy up
  
  data <- tibble::column_to_rownames(data, "submitted_species")
  
  data
}


# function to check if species conversion factorsnumerical values are within pre-defined range 

values_range_check_species <- function(species_data, min_value, max_value) {
  library(tidyverse)
  
  # select conversion factor columns
  conversion_factors <- select(
    species_data, 
    ends_with("_drywt") | ends_with("_lipidwt")
  )
  
  # check the condition if conversion factors are in pre-defined range
  condition <- (conversion_factors >= min_value) & (conversion_factors <= max_value)
  
  # convert NA to TRUE to only have boolean values
  condition[is.na(condition)] <- TRUE
  
  # check if all values are TRUE
  check <- all(condition == TRUE)
  
  if(!check) {
    message(
      "Warning: not all conversion factors in species reference table are within ",
      "range [min_value, max_value] !!!"
    )
  }
  
  invisible()
}



# turn drywt or lipidwt data into a more usable dataframe

ctsm_get_species_cfs <- function(data, wt = c("drywt", "lipidwt")) {
  
  # location: information_functions.R
  # purpose: gets typical species drywt / lipidwt conversion factors and converts 
  #  to a more usable form

  wt <- match.arg(wt)
  wt_adj <- paste0("_", wt)

  data <- rownames_to_column(data, "species")
  data <- select(data, .data$species, ends_with(wt_adj)) 

  
  # deal with null case where there are no conversion factors 
  
  if (ncol(data) == 1) {
    data$matrix <- NA_character_
    data[[wt]] <- NA_real_
    return(data)
  }
    
  
  data <- pivot_longer(
    data, 
    ends_with(wt_adj), 
    names_to = "matrix", 
    values_to = wt, 
    values_drop_na = TRUE
  )
  
  data <- separate(data, .data$matrix, c("matrix", NA), sep = "_") 
  data <- as.data.frame(data)
  
  data
} 
                                   
                                   


# Determinand information and functions ----

ctsm_read_determinand <- function(
  file, 
  path = "information", 
  compartment = c("biota", "sediment", "water"), 
  simplify = TRUE) {
  
  # location: information_functions.r
  # purpose: reads determinand reference table
  # arguments:
  # - compartment only checks for variables relevant to those compartments
  # - simplify only keeps relevant variables  
  
  # argument validation
  
  ok <- is.character(compartment) && length(compartment) %in% 1:3 &&
    all(compartment %in% c("biota", "sediment", "water"))
  if (!ok) {
    stop(
      "\nInvalid argment 'compartment' in function ctsm_read_determinand.\n", 
      "Must be a character string with at least one of 'biota', 'sediment', ",
      "'water'.\n",
      call. = FALSE
    )
  }

  
  # initialise key variables
  
  var_id <- c(
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
    biota_sd_constant = "numeric",
    biota_sd_variable = "numeric",
    sediment_sd_constant = "numeric",
    sediment_sd_variable = "numeric",
    water_sd_constant = "numeric",
    water_sd_variable = "numeric",
    distribution = "character", 
    good_status = "character"
  )
  
  required <- c(
    "determinand", "common_name", "pargroup", "distribution", "good_status"
  )
  
  extra <- c("group", "assess", "unit", "auxiliary")
  extra <- paste(rep(compartment, each = 4), extra, sep = "_")
  required <- c(required, extra)
    
  optional <- c("sd_constant", "sd_variable")
  optional <- paste(rep(compartment, each = 2), optional, sep = "_")

  
  # check required variables are present in data and issue message if optional
  # variables are not
  
  data <- read.csv(
    file.path(path, file), 
    strip.white = TRUE, 
    nrows = 1
  )
  
  ok <- required %in% names(data)
  
  if (!all(ok)) {
    id <- required[!ok]
    id <- sort(id)
    stop(
      "\nThe following variables are missing from ", file, ".\n", 
      "Use the 'compartment' argument to limit the required variables or update\n", 
      "the reference table to continue. Note that variable names are case sensitive.\n",
      "Variables: ", paste(id, collapse = ", "),
      call. = FALSE
    )
  }
  
    
  ok <- optional %in% names(data)
  
  if (!all(ok)) {
    id <- optional[!ok]
    message(
      "The following variables are missing from ", file, ".\n", 
      "They will be created and populated with missing values: missing ", 
      "uncertainties\ncan not be imputed and corresponding measurements will be ", 
      "deleted.\n",
      "Variables: ", paste(id, collapse = ", ")
    )
  }


  ok <- names(var_id) %in% names(data)
  
  data <- read.csv(
    file.path(path, file), 
    na.strings = c("", "NULL"),
    strip.white = TRUE,
    colClasses = var_id[ok]
  )
  

  # fill in common name if missing and create optional variables if missing
  
  data$common_name <- ifelse(
    is.na(data$common_name), 
    data$determinand, 
    data$common_name 
  )

  
  ok <- optional %in% names(data)
  
  if (!all(ok)) {
    id <- optional[!ok]
    data[id] <- rep(NA_real_, nrow(data))
  }


  # check all auxiliary variables are determinands in their own right
  
  lapply(compartment, function(id) {
    
    auxiliary <- data[[paste0(id, "_auxiliary")]] 
    auxiliary <- strsplit(auxiliary, ", ")
    auxiliary <- unlist(auxiliary)

    auxiliary <- unique(na.omit(auxiliary))
    
    ok <- auxiliary %in% data$determinand
    if(!all(ok)) {
      stop(
        '\nNot found in determinand information file: ', 
        paste(auxiliary[!ok], collapse = ", ")
      )
    }
  })
  
  
  # retain relevant variables
  
  ok <- names(var_id) %in% c(required, optional)
  id <- names(var_id)[ok]
  data <- data[id]
  

  # tidy up for output
  
  data <- tibble::column_to_rownames(data, "determinand")
  
  data
}


# extractor functions

ctsm_get_determinands <- function(compartment = c("biota", "sediment", "water")) {
  
  # information_functions.R
  # gets determinands to be assessed from determinand reference table
  
  compartment <- match.arg(compartment)
  
  assess_id <- paste0(compartment, "_assess")
  ok <- info.determinand[[assess_id]]
  
  if (!any(ok)) {
    stop(
      "No determinands have been been selected for assessment:\n", 
      "  please update the determinand reference table or supply the\n",
      "  determinands directly via the determinands argument."
    )
  }
  
  row.names(info.determinand)[ok]
}  


ctsm_get_auxiliary <- function(determinands, info) {
  
  # information_functions.R
  # gets required auxiliary variables for determinands
  
  # in case determinands is a factor
  determinands <- as.character(determinands)
  
  determinands <- unique(determinands)
  
  # auxiliary_id <- paste0(compartment , "_auxiliary")
  # auxiliary <- info.determinand[determinands, auxiliary_id]
  
  auxiliary <- ctsm_get_info(
    info$determinand, 
    determinands, 
    "auxiliary", 
    info$compartment, 
    na_action = "output_ok", 
    sep = "_"
  )
  
  auxiliary <- strsplit(auxiliary, ", ")
  auxiliary <- unlist(auxiliary)
  
  unique(c(na.omit(auxiliary)))
}





# Toxic EQuivalents for WHO_DFP (health)

info_TEQ <- c(
  "CB77" = 0.0001, "CB81" = 0.0003, "CB105" = 0.00003, "CB118" = 0.00003, "CB126" = 0.1, 
  "CB156" = 0.00003, "CB157" = 0.00003, "CB167" = 0.00003, "CB169" = 0.03, 
  "CDD1N" = 1, "CDD4X" = 0.1, "CDD6P" = 0.01, "CDD6X" = 0.1, "CDD9X" = 0.1, "CDDO" = 0.0003,
  "CDF2N" = 0.3, "CDF2T" = 0.1, "CDF4X" = 0.1, "CDF6P" = 0.01, "CDF6X" = 0.1, "CDF9P" = 0.01,
  "CDF9X" = 0.1, "CDFO" = 0.00003, "CDFP2" = 0.03, "CDFX1" = 0.1, "TCDD" = 1
)



# Assessment criteria ----

ctsm_read_assessment_criteria <- function(files, path = "information")  {
  
  out <- sapply(names(files), function(compartment) {
    read.csv(
      file.path(path, files[[compartment]]), 
      na.strings = "",
      strip.white = TRUE
    )},
    simplify = FALSE, 
    USE.NAMES = TRUE
  )
    
  if ("sediment" %in% names(out)) {
    id <- is.na(out$sediment$country)
    out$sediment$country[id] <- ""
  }

  if ("biota" %in% names(out)) {
    wk <- strsplit(out$biota$sub.family, ",", fixed = TRUE)
    n <- sapply(wk, length)
    out$biota <- out$biota[rep(1:nrow(out$biota), times = n),]
    out$biota$sub.family <- unlist(wk)
  }

  out
}


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
    mutate(
      species_group = ctsm_get_info("species", species, "species_group"),
      sub.family = ctsm_get_info("species", species, "species_subgroup")
    ) 
  
  lipid_info <- ctsm_get_species_cfs(info.species, "lipidwt")
  
  drywt_info <- ctsm_get_species_cfs(info.species, "drywt")
  
  data <- left_join(data, lipid_info, by = c("species", "matrix"))
  data <- left_join(data, drywt_info, by = c("species", "matrix"))
  
  data <- data[c("rownames", "determinand", "sub.family", "basis", "drywt", "lipidwt")]
  
  data <- left_join(data, AC_data, by = c("determinand", "sub.family"))
  
  out <- sapply(AC, simplify = FALSE, FUN = function(i) {
    basis_AC <- paste("basis", i, sep = ".")
    ctsm_convert_basis(
      data[[i]], 
      from = data[[basis_AC]], 
      to = data$basis, 
      drywt = data$drywt, 
      lipidwt = data$lipidwt
    )
  })
  
  out <- data.frame(out)
  
  if (export_cf) {
    out <- bind_cols(out, data[c("drywt", "lipidwt")])
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
          .data$sub.family %in% "Oyster"                             ~ NA_real_,
          .data$species_group %in% "Mammal" & .data$matrix %in% "LI" ~
            ctsm_convert_basis(16000, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Mammal" & .data$matrix %in% "HA" ~
            ctsm_convert_basis(6100, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "EH"   ~
            ctsm_convert_basis(110, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "LI"   ~
            ctsm_convert_basis(1400, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "FE"   ~
            ctsm_convert_basis(1580, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "BL"   ~
            ctsm_convert_basis(200, "W", .data$basis, .data$drywt, .data$lipidwt),
          TRUE                                                       ~ .data$BAC
        ),
        
        EQS.OSPAR = case_when(
          .data$species_group %in% "Mammal" & .data$matrix %in% "LI" ~
            ctsm_convert_basis(64000, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Mammal" & .data$matrix %in% "HA" ~
            ctsm_convert_basis(24400, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "EH"   ~
            ctsm_convert_basis(470, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "LI"   ~
            ctsm_convert_basis(7300, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "FE"   ~
            ctsm_convert_basis(7920, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "BL"   ~
            ctsm_convert_basis(1000, "W", .data$basis, .data$drywt, .data$lipidwt),
          TRUE                                                       ~ .data$EQS.OSPAR
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
          (is.na(.data$lipidwt) | .data$lipidwt < lipid_high),
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
        .data$species_group %in% "Fish" & (is.na(.data$lipidwt) | .data$lipidwt < lipid_high),
        NA_real_,
        .data$BAC
      ),
      
      HQS = if_else(
        .data$species_group %in% "Fish" & .data$matrix %in% "LI" & .data$determinand %in% "SCB6",
        ctsm_convert_basis(200, "W", .data$basis, .data$drywt, .data$lipidwt),
        .data$HQS
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "SCB7",
        ctsm_convert_basis(6.7, "W", .data$basis, .data$drywt, .data$lipidwt),
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
        .data$species_group %in% "Fish" & (is.na(.data$lipidwt) | .data$lipidwt < lipid_high),
        NA_real_,
        .data$BAC
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "HCB",
        ctsm_convert_basis(2.0, "W", .data$basis, .data$drywt, .data$lipidwt),
        .data$EAC
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "HCH",
        ctsm_convert_basis(2.0, "W", .data$basis, .data$drywt, .data$lipidwt),
        .data$EAC
      )
    )
    
    
    # HCHG HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "HCHG" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_lipidwt"],
        .dry_mu = info.species[.data$species, "MU_drywt"],
        HQS = ctsm_convert_basis(61, "W", "L", .dry_mu, .lipid_mu),
        HQS = ctsm_convert_basis(.data$HQS, "L", .data$basis, .data$drywt, .data$lipidwt),
        .lipid_mu = NULL,
        .dry_mu = NULL
      )
    }
    
    
    # HCB HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "HCB" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_lipidwt"],
        .dry_mu = info.species[.data$species, "MU_drywt"],
        HQS = ctsm_convert_basis(10, "W", "L", .dry_mu, .lipid_mu),
        HQS = ctsm_convert_basis(.data$HQS, "L", .data$basis, .data$drywt, .data$lipidwt),
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
        .lipid_mu = info.species[.data$species, "MU_lipidwt"],
        .dry_mu = info.species[.data$species, "MU_drywt"],
        HQS = ctsm_convert_basis(0.0085, "W", "L", .dry_mu, .lipid_mu),
        HQS = ctsm_convert_basis(.data$HQS, "L", .data$basis, .data$drywt, .data$lipidwt),
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
        ctsm_convert_basis(0.02, "W", .data$basis, .data$drywt, .data$lipidwt),
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
      
      stopifnot("method_analysis" %in% names(data))
      
      if ("BAC" %in% AC) {
        id <- species %in% "Limanda limanda"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 16
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 3.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 0.15
        
        id <- species %in% "Gadus morhua"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 21
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 2.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 1.1
        
        id <- species %in% "Platichthys flesus"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 16
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 3.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 1.3
        
        id <- species %in% "Melanogrammus aeglefinus"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 13
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 0.8
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 1.9
      }
      
      if ("EAC" %in% AC) {
        id <- species %in% "Limanda limanda"
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 22
        
        id <- species %in% "Gadus morhua"
        out$EAC[id & determinand %in% "PYR1OH" & method_analysis %in% "GC-MS"] <- 483
        out$EAC[id & determinand %in% "PA1OH" & method_analysis %in% "GC-MS"] <- 528
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 35
        
        id <- species %in% "Platichthys flesus"
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 29
        
        id <- species %in% "Melanogrammus aeglefinus"
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 35
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
        out$EAC[determinand %in% "VDS" & species %in% "Tritia nitida (reticulata)"] <- 0.3
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
      
      stopifnot("method_analysis" %in% names(data))
      
      if ("EAC" %in% AC) {
        out$EAC[determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 483
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
        out$EAC[determinand %in% "VDS" & species %in% "Tritia nitida (reticulata)"] <- 0.3
        out$EAC[determinand %in% "VDS" & species %in% "Buccinum undatum"] <- 0.3
        out$EAC[determinand %in% "VDS" & species %in% "Peringia ulvae"] <- 0.1
        out$EAC[determinand %in% "INTS" & species %in% "Littorina littorea"] <- 0.3
      }
      
      out
    })
  }

}


if (info_AC_type == "EXTERNAL") {
  
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
            ctsm_convert_basis(16000, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Mammal" & .data$matrix %in% "HA" ~
            ctsm_convert_basis(6100, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "EH"   ~
            ctsm_convert_basis(110, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "LI"   ~
            ctsm_convert_basis(1400, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "FE"   ~
            ctsm_convert_basis(1580, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "BL"   ~
            ctsm_convert_basis(200, "W", .data$basis, .data$drywt, .data$lipidwt),
          TRUE                                                ~ .data$BAC
        ),
        
        EQS.OSPAR = case_when(
          .data$species_group %in% "Mammal" & .data$matrix %in% "LI" ~
            ctsm_convert_basis(64000, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Mammal" & .data$matrix %in% "HA" ~
            ctsm_convert_basis(24400, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "EH"   ~
            ctsm_convert_basis(470, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "LI"   ~
            ctsm_convert_basis(7300, "W", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "FE"   ~
            ctsm_convert_basis(7920, "D", .data$basis, .data$drywt, .data$lipidwt),
          .data$species_group %in% "Bird" & .data$matrix %in% "BL"   ~
            ctsm_convert_basis(1000, "W", .data$basis, .data$drywt, .data$lipidwt),
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
          (is.na(.data$lipidwt) | .data$lipidwt < lipid_high),
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
        .data$species_group %in% "Fish" & (is.na(.data$lipidwt) | .data$lipidwt < lipid_high),
        NA_real_,
        .data$BAC
      ),
      
      HQS = if_else(
        .data$species_group %in% "Fish" & .data$matrix %in% "LI" & .data$determinand %in% "SCB6",
        ctsm_convert_basis(200, "W", .data$basis, .data$drywt, .data$lipidwt),
        .data$HQS
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "SCB7",
        ctsm_convert_basis(6.7, "W", .data$basis, .data$drywt, .data$lipidwt),
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
        .data$species_group %in% "Fish" & (is.na(.data$lipidwt) | .data$lipidwt < lipid_high),
        NA_real_,
        .data$BAC
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "HCB",
        ctsm_convert_basis(2.0, "W", .data$basis, .data$drywt, .data$lipidwt),
        .data$EAC
      ),
      
      EAC = if_else(
        .data$species %in% c("Sterna hirundo", "Haematopus ostralegus") &
          .data$matrix %in% "EH" & .data$determinand %in% "HCH",
        ctsm_convert_basis(2.0, "W", .data$basis, .data$drywt, .data$lipidwt),
        .data$EAC
      )
    )
    
    
    # HCHG HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "HCHG" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_lipidwt"],
        .dry_mu = info.species[.data$species, "MU_drywt"],
        HQS = ctsm_convert_basis(61, "W", "L", .dry_mu, "", .lipid_mu, ""),
        HQS = ctsm_convert_basis(.data$HQS, "L", .data$basis, .data$drywt, .data$lipidwt),
        .lipid_mu = NULL,
        .dry_mu = NULL
      )
    }
    
    
    # HCB HQS in liver converted from muscle using muscle lipid content
    
    id <- out$determinand %in% "HCB" & out$species_group %in% "Fish" & out$matrix %in% "LI"
    
    if (any(id)) {
      
      out[id, ] <- mutate(
        out[id, ],
        .lipid_mu = info.species[.data$species, "MU_lipidwt"],
        .dry_mu = info.species[.data$species, "MU_drywt"],
        HQS = ctsm_convert_basis(10, "W", "L", .dry_mu, "", .lipid_mu, ""),
        HQS = ctsm_convert_basis(.data$HQS, "L", .data$basis, .data$drywt, .data$lipidwt),
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
        .lipid_mu = info.species[.data$species, "MU_lipidwt"],
        .dry_mu = info.species[.data$species, "MU_drywt"],
        HQS = ctsm_convert_basis(0.0085, "W", "L", .dry_mu, "", .lipid_mu, ""),
        HQS = ctsm_convert_basis(.data$HQS, "L", .data$basis, .data$drywt, .data$lipidwt),
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
        ctsm_convert_basis(0.02, "W", .data$basis, .data$drywt, .data$lipidwt),
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
      
      stopifnot("method_analysis" %in% names(data))
      
      if ("BAC" %in% AC) {
        id <- species %in% "Limanda limanda"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 16
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 3.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 0.15
        
        id <- species %in% "Gadus morhua"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 21
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 2.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 1.1
        
        id <- species %in% "Platichthys flesus"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 16
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 3.7
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 1.3
        
        id <- species %in% "Melanogrammus aeglefinus"
        out$BAC[id & determinand %in% "PYR1OH" & method_analysis %in% "HPLC-FD"] <- 13
        out$BAC[id & determinand %in% "PA1OH" & method_analysis %in% "HPLC-FD"] <- 0.8
        out$BAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 1.9
      }
      
      if ("EAC" %in% AC) {
        id <- species %in% "Limanda limanda"
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 22
        
        id <- species %in% "Gadus morhua"
        out$EAC[id & determinand %in% "PYR1OH" & method_analysis %in% "GC-MS"] <- 483
        out$EAC[id & determinand %in% "PA1OH" & method_analysis %in% "GC-MS"] <- 528
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 35
        
        id <- species %in% "Platichthys flesus"
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 29
        
        id <- species %in% "Melanogrammus aeglefinus"
        out$EAC[id & determinand %in% "PYR1OHEQ" & method_analysis %in% "FLM-SS"] <- 35
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
        out$EAC[determinand %in% "VDS" & species %in% "Tritia nitida (reticulata)"] <- 0.3
        out$EAC[determinand %in% "VDS" & species %in% "Buccinum undatum"] <- 0.3
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

ctsm_convert_basis <- function(
  conc, from, to, 
  drywt = NA_real_, 
  lipidwt = NA_real_, 
  drywt_censoring = ifelse(is.na(drywt), NA_character_, ""), 
  lipidwt_censoring = ifelse(is.na(lipidwt), NA_character_, ""), 
  exclude = FALSE, 
  print_warning = TRUE) {

  # source: information_functions.R
  # dependencies: dplyr
  
  # converts concentrations between wet, dry and lipid bases

  # from ~ current basis
  # to ~ target basis

  # drywt and lipidwt assumed to be percentages taking values in [0, 100]

  # lipidwt assumed to be on a wet weight basis - could generalise (issue raised)
  
  # exclude is a logical identifying records that do not need to be converted -
  # useful when the data contain biological effects measurements
  
  # print_warning gives the number of records lost during conversion
  
  
  # set up working data frame
  
  data <- data.frame(
    conc, from, to, drywt, drywt_censoring, lipidwt, lipidwt_censoring, exclude 
  )

  
  # ensure variables are stored as characters (rather than factors)
  
  var_id <- c("from", "to", "drywt_censoring", "lipidwt_censoring")
  data <- dplyr::mutate(data, dplyr::across(all_of(var_id), as.character))

  
  # check arguments have admissible values
  
  ctsm_check_convert_basis(data)
  

  # only need to convert records which are 
  # a) not excluded
  # b) concentration is not missing
  # c) from and to differ
  
  data$exclude <- data$exclude | 
    is.na(data$conc) |
    (!is.na(data$from) & !is.na(data$to) & data$from == data$to)
  
  if (all(data$exclude)) {
    return(data$conc)
  }


  # set up extra variables to do the conversion

  data <- dplyr::mutate(
    data, 
    from_value = dplyr::case_when(
      .data$from %in% "W" ~ 100, 
      .data$from %in% "D" ~ drywt, 
      .data$from %in% "L" ~ lipidwt,
      TRUE                ~ NA_real_
    ),
    from_censoring = dplyr::case_when(
      .data$from %in% "W" ~ "", 
      .data$from %in% "D" ~ drywt_censoring, 
      .data$from %in% "L" ~ lipidwt_censoring,
      TRUE                ~ NA_character_ 
    ),
    to_value = dplyr::case_when(
      .data$to %in% "W" ~ 100, 
      .data$to %in% "D" ~ drywt, 
      .data$to %in% "L" ~ lipidwt,
      TRUE              ~ NA_real_
    ),
    to_censoring = dplyr::case_when(
      .data$to %in% "W" ~ "", 
      .data$to %in% "D" ~ drywt_censoring, 
      .data$to %in% "L" ~ lipidwt_censoring,
      TRUE              ~ NA_character_ 
    )
  )
  
  
  # conversion not possible if missing or censored dry or lipid data or because 
  # the to_value is 0%
  
  data <- dplyr::mutate(
    data, 
    convert = !.data$exclude,
    convert_ok = .data$convert & 
      !is.na(.data$to_value) & !is.na(.data$from_value) & 
      data$from_censoring %in% "" & .data$to_censoring %in% "" &
      .data$to_value > 0
  )  
  
  if (sum(data$convert_ok) < sum(data$convert) && print_warning) {
    message(
      "   Losing ", sum(data$convert) - sum(data$convert_ok), " out of ", 
      nrow(data), " records in basis conversion due to missing, censored\n",
      "   or zero drywt or lipidwt values."
    )
  }
    
  
  # finally convert
  
  conc = dplyr::case_when(
    !data$convert    ~ data$conc,
    !data$convert_ok ~ NA_real_,
    TRUE             ~ data$conc * data$from_value / data$to_value
  )

  conc
}


ctsm_check_convert_basis <- function(data) {

  # source: information_functions.R
  # dependencies: dplyr

  # check inputs to ctsm_convert_basis are valid
  
  if (any(is.na(data$exclude))) {
    stop("Missing values not allowed in 'exclude'")
  }
  

  # only need to consider values in records which are 
  # a) not excluded
  # b) concentration is not missing
  # c) from and to differ

  data$exclude <- data$exclude | 
    is.na(data$conc) |
    (!is.na(data$from) & !is.na(data$to) & data$from == data$to)

  if (all(data$exclude)) {
    return(invisible())
  }
  
  data <- data[!data$exclude, ]
  
  
  # from and to 

  not_ok <- c(any(is.na(data$from)), any(is.na(data$to)))
  
  if (any(not_ok)) {
    txt <- c('from', 'to')[not_ok]
    txt <- paste(txt, collapse = ", and ")
    warning(
      "Missing basis values in ", txt,  
      "; associated concentrations will be set to missing.\n", 
      "The exclude argument can prevent basis conversion (and ", 
      "suppress this warning);\n",
      "this would be appropriate for biological effects or auxiliary variables.", 
      call. = FALSE, immediate. = TRUE
    )
  }
  
  if (!all(data$from %in% c("W", "D", "L", NA_character_))) {
    stop("Unexpected basis values in 'from': must be 'W', 'D', 'L' or NA")
  }  

  if (!all(data$to %in% c("W", "D", "L", NA_character_))) {
    stop("Unexpected basis values in 'to': must be 'W', 'D', 'L' or NA")
  }  
  

  # drywt and lipidwt
  
  if (!all(is.na(data$drywt) | dplyr::between(data$drywt, 0, 100))) {
    stop("Some drywt values are less than 0% or greater than 100%")
  }

  if (!all(is.na(data$lipidwt) | dplyr::between(data$lipidwt, 0, 100))) {
    stop("Some lipidwt values are less than 0% or greater than 100%")
  }
  

  # validate drywt_censoring and lipidwt_censoring

  if (any(!is.na(data$drywt) & is.na(data$drywt_censoring))) {
    warning(
      "Missing drywt_censoring values when drywt is present;\n",  
      "this might result in data getting lost during basis conversion.",
      call. = FALSE
    )
  }
    
  if (any(!is.na(data$lipidwt) & is.na(data$lipidwt_censoring))) {
    warning(
      "Missing lipidwt_censoring values when lipidwt is present;\n",  
      "this might result in data getting lost during basis conversion.",
      call. = FALSE
    )
  }
  
  if (!all(data$drywt_censoring %in% c("", "<", "D", "Q", NA_character_))) {
    stop(
      "Unexpected drywt_censoring values: must be '', '<', 'D', 'Q' or NA.\n",
      "If you have greater_than values, contact the HARSAT development team."
    )
  }  
  
  if (!all(data$lipidwt_censoring %in% c("", "<", "D", "Q", NA_character_))) {
    stop(
      "Unexpected lipidwt_censoring values: must be '', '<', 'D', 'Q' or NA.\n",
      "If you have greater_than values, contact the HARSAT development team."
    )
  }  

  invisible()
}




# four get_basis functions defined here

# default is very simplistic, but works in all cases for sediment and water
# most_common was used by AMAP in their mercury assessment and takes the most 
#   common basis reported in each station, species (biota), matrix and  
#   determinand_group combination
# biota_OSPAR is the current (2023) OSPAR biota configuration
# biota_HELCOM is the current (2023) HELCOM biota configuration


get_basis_default <- function(data, info) {
  
  # gets default target basis - information_functions.r
  
  # biota: W
  # sediment: D
  # water: W
  
  # the exceptions are biological effects measurements where it is assumed the 
  # data are submitted on the correct basis (or where basis isn't relevant)
  
  basis_id <- switch(
    info$compartment, 
    biota = "W", 
    sediment = "D",
    water = "W"
  )
  
  new_basis <- dplyr::if_else(
    data$group %in% c("Imposex", "Metabolites", "Effects"), 
    NA_character_, 
    basis_id
  )
  
  new_basis  
}


get_basis_most_common <- function(data, info) {

  # gets target basis defined as the most commonly reported basis 
  # information_function.r
  
  # the basis most reported within a particular station, species, matrix and 
  # determinand group (regardless of whether auxiliary variables are present to 
  # enable conversion)
  
  # get grouping identifier
  
  var_id <- c("station", "species", "matrix", "group")
  var_id <- intersect(var_id, names(data))
  
  data$.id <- do.call("paste", c(data[var_id], sep = "_"))
  
  
  # provide index to ensure output is in same order as original
  
  data$.order <- 1:nrow(data)
  
  
  # get modal basis within each group
  
  out <- by(data, data$.id, function(x) {
    
    # deal with e.g. biological effects which don't have a basis
    
    if (unique(x$group) %in% c("Metabolites", "Imposex", "Effects")) {
      x$new_basis <- rep(NA_character_, nrow(x))
      x <- x[c(".order", "new_basis")]
      return(x)
    }
    
    # check that have full basis information
    
    if (any(is.na(x$basis))) {
      stop("missing basis information for the following: ", unique(x$.id))
    }
    
    wk <- table(x$basis)
    x$new_basis = names(wk)[which.max(wk)]
    x[c(".order", "new_basis")]
  })
  
  out <- do.call(rbind, out)
  
  
  # return to the original ordering
  
  out <- out[order(out$.order), ]
  
  out$new_basis
}



get_basis_biota_OSPAR <- function(data, info) {
  
  # 2023 OSPAR biota target basis - information_functions.r
  
  # note hard-wiring of lipid_high which should be passed as a control variable
  
  
  if (info$compartment != "biota"){
    stop("Incorrect compartment specified")
  }
  
  
  # define cut-off for using lipid as a basis 
  
  lipid_high <- 3.0
  

  # combine input variables and get species family
  
  out <- data[c("species", "matrix", "determinand", "group")]
  
  out <- mutate(
    out,
    across(everything(), as.character),
    species_group = ctsm_get_info(info$species, .data$species, "species_group")
  )
  

  # get typical lipid content by species and matrix 
  
  lipid_info <- ctsm_get_species_cfs(info$species, "lipidwt")
  
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
    new_basis = case_when(
      .data$group %in% c("Imposex", "Effects", "Metabolites")            ~ NA_character_,
      .data$species_group %in% c("Bivalve", "Gastropod")               ~ "D",
      .data$species_group %in% c("Fish", "Crustacean") & 
        .lw &
        .data$lipidwt >= lipid_high                                     ~ "L",
      .data$species_group %in% c("Fish", "Crustacean")                    ~ "W",
      .data$species_group %in% "Mammal" & .data$matrix %in% "HA"         ~ "D",
      .data$species_group %in% "Mammal" & .data$group %in% "Metals"      ~ "W",
      .data$species_group %in% "Mammal"                                  ~ "L",
      .data$species %in% c("Alle alle", "Rissa tridactyla")              ~ "D",
      .data$matrix %in% "EH" & .data$group %in% "Metals" & 
        .data$species %in% "Larus argentatus"                            ~ "D",
      .data$matrix %in% "EH" & .data$group %in% "Metals"                 ~ "W",
      .data$matrix %in% "EH" & .data$group %in% "Organofluorines"        ~ "W",
      .data$matrix %in% "EH" & 
        .data$species %in% c(
          "Cepphus grylle", "Haematopus ostralegus", "Sterna hirundo"
        )                                                                ~ "L",
      .data$matrix %in% "EH" & 
        .data$species %in% c("Larus argentatus", "Somateria mollissima") ~ "W",
      .data$species_group %in% "Bird"                                    ~ "W"
    )
  )
  
  out$new_basis
}



# Matrix ----

ctsm_read_matrix <- function(file, path = "information") {
  read.csv(
    file.path(path, "matrix.csv"), 
    row.names = "matrix", 
    strip.white = TRUE
  )
}

info.matrix <- ctsm_read_matrix("matrix.csv")


# Regions ----

info.regions <- sapply(
  c("AMAP", "HELCOM", "OSPAR"), 
  function(x) {
    infile <- paste(x, "regions.csv")
    infile <- file.path("information", infile)
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

ctsm_read_method_extraction <- function(file, path = "information") {
  read.csv(
    file.path(path, file), 
    row.names = "METCX",  
    na.strings = "",
    strip.white = TRUE
  )
}

info.methodExtraction <- ctsm_read_method_extraction("method of extraction.csv") 


ctsm_read_pivot_values <- function(file, path = "information") {
  read.csv(
    file.path(path, file), 
    na.strings = "",
    strip.white = TRUE
  )
}

info.pivotValues <- ctsm_read_pivot_values("pivot values.csv")



# Imposex ----

info.imposex <- read.csv(file.path("information", "imposex.csv"))

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
  
get_station_code <- function(station_name, country, stations) {
  
  # gets the station code corresponding to the station name and country from the 
  # station dictionary
  
  # only works for one country at a time
  
  stopifnot(length(country) == 1)
  n <- length(station_name)
  
  out <- data.frame(station_name = station_name, country = country) 
  out <- mutate(out, across(.fns = as.character))
    
  
  stations <- stations[c("station_name", "country", "station_code")]
  stations <- mutate(stations, across(.fns = as.character)) 
  
  out <- left_join(out, stations, by = c("station_name", "country"))
  
  stopifnot(!is.na(out), n == nrow(out))
  
  out$station_code
}
