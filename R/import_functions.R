# main import and data processing functions ----


# read in data  ---- 

#' Read harsat data
#'
#' @param compartment a string, one of "biota", "sediment", and "water"
#' @param purpose a string, one of "OSPAR", "HELCOM", "AMAP" and "other"
#' @param contaminants a file reference for contaminant data
#' @param stations a file reference for station data
#' @param QA a file reference for QA data
#' @param data_dir The directory where the data files can be found (sometimes
#'   supplied using 'file.path'). Defaults to "."; i.e. the working directory.
#' @param info_files A list of files specifying reference tables which override
#'   the defaults. See examples.
#' @param info_dir The directory where the reference tables can be found
#'   (sometimes supplied using 'file.path'). Defaults to "."; i.e. the working
#'   directory
#' @param extraction A date saying when the extraction was made. Optional. This
#'   should be provided according to ISO 8601; for example, 29 February 2024
#'   should be supplied as "2024-02-29".
#' @param max_year
#' @param oddity_dir The directory where the 'oddities' will be written
#'   (sometimes supplied using 'file.path'). This directory (and subdirectories)
#'   will be created if it does not already exist.
#' @param control
#' @param data_format a string, one of "ICES_old", "ICES_new", and "external"
#'
#' @export
ctsm_read_data <- function(
  compartment = c("biota", "sediment", "water"), 
  purpose = c("OSPAR", "HELCOM", "AMAP", "other"), 
  contaminants, 
  stations, 
  QA, 
  data_dir = ".", 
  info_files = list(),
  info_dir = ".",
  extraction = NULL, 
  max_year = NULL,
  oddity_dir = "oddities",
  control = list(), 
  data_format = c("ICES_old", "ICES_new", "external")) {

  # import_functions.R
  # dependencies lubridate
  
  # reads in data from an ICES extraction or from an external data file
  
  
  # validate arguments
  
  compartment <- match.arg(compartment)
  purpose <- match.arg(purpose)
  data_format <- match.arg(data_format)

  if (!is.null(max_year)) {
    if (length(max_year) != 1) {
      stop("max_year must be a single integer valued year")
    } 
    if (!is.integer(max_year)) {
      if (!isTRUE(all.equal(max_year, as.integer(max_year)))) {
        stop("max_year must be an integer valued year")
      }
      max_year <- as.integer(max_year)
    }
  }    
    
  if (!is.null(extraction)) {
    if (length(extraction) > 1 | !is.character(extraction)) {
      stop(
        "extraction must be a single character string of the form ymd: e.g. \"",
        lubridate::today(), ""
      )
    }
  }  
  
  
  # turn extraction into date object

  if (!is.null(extraction)) {
    extraction <- suppressWarnings(lubridate::ymd(extraction))
    if (is.na(extraction)) {
      stop(
        "extraction not recognised: it must be a valid date of the form ", 
        "ymd: e.g. \"", lubridate::today(), "\""
      )
    }  
  }  

    
  # set up control (info) structure 
  # some of this should probably go in ctsm_create_timeSeries
    
  info <- list(
    compartment = compartment, 
    purpose = purpose, 
    extraction = extraction,
    data_format = data_format, 
    max_year = max_year, 
    oddity_dir = oddity_dir
  )
  
  control_default <- ctsm_control_default(purpose, compartment)
  
  control <- ctsm_control_modify(control_default, control)
  
  if (any(names(control) %in% names(info))) {
    warning("possible conflict between function arguments and control")
  }
  
  info <- append(info, control)
  
  
  # read in reference tables
  
  info <- ctsm_read_info(info, info_dir, info_files)
  

  # read in station dictionary, contaminant and biological effects data and QA data
  
  stations <- ctsm_read_stations(stations, data_dir, info)
  
  data <- ctsm_read_contaminants(contaminants, data_dir, info)
  
  if (data_format == "ICES_old") {
    QA <- ctsm_read_QA(QA, data_dir, purpose)
  }
  

  # populate max_year if not set in function call
  
  if (is.null(info$max_year)) {
    info$max_year <- max(data$year)
    cat(
      "\nArgument max_year taken to be the maximum year in the data:", 
      info$max_year
    )
  }


  # construct recent_years in which there must be some monitoring data
  
  info$recent_years <- seq(
    info$max_year - info$reporting_window + 1, 
    info$max_year
  )
  
    
  out <- list(
    call = match.call(), 
    info = info,
    data = data, 
    stations = stations
  )
  
  if (data_format == "ICES_old") out$QA <- QA
  
  out
}


ctsm_control_default <- function(purpose, compartment) {
  
  # import functions
  # sets up default values that control the assessment
  
  # reporting_window is set to 6 to match the MSFD reporting cycle
  # series with no data in the most recent_window are excluded
  
  # region_id are the names or the relevant regions in the ICES extraction
  # region_names are the names that will be used for reporting 
  
  # all_in_region is a logical that determines whether all data (and stations)
  # must be in a region

  # bivalve_spawning_season is a character vector of months when contaminant
  # data for bivalves and gastropds will be deleted because they are in the 
  # spawning season
  
  region_id <- switch(
    purpose, 
    OSPAR = c("OSPAR_region", "OSPAR_subregion"),
    HELCOM = c("HELCOM_subbasin", "HELCOM_L3", "HELCOM_L4"),
    NULL
  )
  
  region_names <- region_id
  
  all_in_region <- !is.null(region_id)
    
  bivalve_spawning_season <- switch(
    compartment, 
    biota = switch(
      purpose, 
      OSPAR = c("April", "May", "June", "July"),
      HELCOM = c("April", "May", "June", "July"),
      NULL
    ), 
    NULL
  )
  
  list(
    reporting_window = 6L, 
    region_id = region_id,
    region_names = region_names,
    all_in_region = all_in_region, 
    bivalve_spawning_season = bivalve_spawning_season
  )
}


ctsm_control_modify <- function(control_default, control) {

  # import functions
  # updates default control structure with user specification and does basic 
  # error checking
  
  # initialise region_names if region_id has been modified (and region_names is
  # not specified)
  
  if ("region_id" %in% names(control) & !("region_names" %in% names(control))) {
    control <- append(control, list(region_names = control$region_id))
  }
  
  control <- modifyList(control_default, control, keep.null = TRUE)

  if (length(control$region_id) != length(control$region_names)) {
    stop(
      "error in control argument: length of region_id and region_names must be ",
      "identical"
    )
  }
  
  if (!is.null(control$bivalve_spawning_season)) {
    months <- c(
      "January", "February", "March", "April", "May", "June", "July", "August",
      "September", "October", "November", "December"
    )
    if (!all(control$bivalve_spawning_season %in% months)) {
      stop(
        "error in control argument: invalid months in bivalve_spawning_season"
      )
    }
  }
  
  control
}



ctsm_read_info <- function(info, path, info_files) {
  
  # location: import_functions.R
  # purpose: reads in reference tables, using default files unless overruled 
  #  by info_control
  
  # default reference tables
  
  if (info$purpose == "OSPAR") {
    files <- list(
      determinand = "determinand_OSPAR_2022.csv", 
      species = "species_OSPAR_2022.csv",
      region_values = "regions_OSPAR.csv", 
      thresholds = paste0("thresholds_", info$compartment, "_OSPAR_2022.csv")
    )
  } else if (info$purpose == "HELCOM") {
    files <- list(
      determinand = "determinand_HELCOM_2023.csv", 
      species = "species_HELCOM_2023.csv",
      region_values = NULL,
      thresholds = paste0("thresholds_", info$compartment, "_HELCOM_2023.csv")
    )
  } else if (info$purpose == "AMAP") {
    files <- list(
      determinand = "determinand_AMAP_2022.csv", 
      species = "species_AMAP_2022.csv",
      region_values = NULL,
      thresholds = if (info$compartment == "biota") {
        "thresholds_biota_AMAP.csv"}
      else {
        NULL
      }
    )
  } else {
    files <- list(
      determinand = "determinand_default.csv", 
      species = "species_default.csv",
      region_values = NULL, 
      thresholds = NULL
    )
  }

  files$method_extraction <- "method_extraction.csv"
  files$pivot_values <- "pivot_values.csv"
  files$matrix <- "matrix.csv"
  files$imposex <- "imposex.csv"
  
  # modify with user supplied files
  
  files <- modifyList(files, info_files, keep.null = TRUE)
  

  # check determinand and species reference tables have not been made NULL
  # threshold values do not need to be specified (if only interested in trends)
  # other files unlikely to be altered by users at this stage

  if (is.null(files$determinand)) {
    stop(
      "\nDeterminand reference table not specified.\n",
      "Check for an error in argument 'info_files'",
      call. = FALSE
    )
  }
  
  if (info$compartment == "biota" && is.null(files$species)) {
    stop(
      "\nSpecies reference table not specified.\n",
      "Check for an error in argument 'info_files'",
      call. = FALSE
    )
  }
  

  info$determinand <- ctsm_read_determinand(
    files$determinand, path, info$compartment
  )
  
  info$matrix <- ctsm_read_matrix(files$matrix, path)
  
  if (info$compartment == "biota") {
    info$species <- ctsm_read_species(files$species, path)
    info$imposex <- ctsm_read_imposex(files$imposex, path)
  }

  if (info$compartment == "sediment") {
    info$method_extraction <- ctsm_read_method_extraction(
      files$method_extraction, path
    ) 
    info$pivot_values <- ctsm_read_pivot_values(files$pivot_values, path)
  }  
      
  if (!is.null(files$region_values)) {
    info$region_values <- ctsm_read_regions(
      files$region_values, path, info$purpose
    )
  } 
  
  if (!is.null(files$thresholds)) {
    info$thresholds <- ctsm_read_thresholds(
      files$thresholds, path, info$compartment
    )
  }
  
  info
}


ctsm_read_stations <- function(file, path, info) {

  # import functions
  # read in station dictionary
  
  infile <- file.path(path, file)
  cat("Reading station dictionary from:\n '", infile, "'\n", sep = "")

  
  if (info$data_format %in% c("ICES_old", "ICES_new")) {
    
    stations <- read.csv(
      infile, na.strings = c("", "NULL"), strip.white = TRUE
    )
    
    # rename columns to suit!
    
    names(stations) <- gsub("Station_", "", names(stations), fixed = TRUE)
    
    stations <- dplyr::rename(
      stations, 
      station_code = Code,
      country_ISO = Country,
      country = Country_CNTRY,
      station_name = Name,
      station_longname = LongName,
      station_latitude = Latitude,
      station_longitude = Longitude,
      startYear = ActiveFromDate,
      endYear = ActiveUntilDate,
      parent_code = AsmtMimeParent,
      replacedBy = ReplacedBy,
      programGovernance = ProgramGovernance,
      dataType = DataType,
      station_type = MSTAT,
      waterbody_type = WLTYP
    )
    
    if (info$purpose %in% c("OSPAR", "AMAP")) {
      stations <- dplyr::rename(stations, offshore = OSPAR_shore)
    }  
  
    # turn following variables into a character
    # code and parent code could be kept as integers, but safer to leave as 
    # characters until have decided how we are going to merge with non-ICES data
    
    id <- c(info$region_id, "station_code", "parent_code")
  
    stations <- dplyr::mutate(
      stations, 
      dplyr::across(any_of(id), as.character)
    )

  }    
        
  
  if (info$data_format == "external") {  
    
    var_id<- c(
      country = "character",
      station_name = "character",
      station_code = "character",
      station_longname = "character",
      station_latitude = "numeric", 
      station_longitude = "numeric",
      station_type = "character",
      waterbody_type = "character"
    )
    
    if (!is.null(info$region_id)) {
      extra <- rep("character", length(info$region_id))
      names(extra) <- info$region_id
      var_id <- c(var_id, extra)
    }
    
    required <- c(
      "country", "station_name", "station_code", "station_latitude", 
      "station_longitude", info$region_id
    )
    
    
    # check required variables are present in data
    
    stations <- read.csv(infile, strip.white = TRUE, nrows = 1)
    
    ok <- required %in% names(stations)
    
    if (!all(ok)) {
      id <- required[!ok]
      id <- sort(id)
      stop(
        "The following variables are not in the stations file. ", 
        "Please update the stations file to continue. ",
        "Note that the variable names are case sensitive.\n",
        "Variables: ", paste(id, collapse = ", ")
      )
    }
    
    
    # read data
    
    ok <- names(var_id) %in% names(stations)
    
    stations <- read.csv(
      infile, 
      na.strings = c("", "NULL"),
      strip.white = TRUE,
      colClasses = var_id[ok]
    )

    
    # create missing (non-required) variables 
    
    id <- c("station_longname", "station_type", "waterbody_type")

    for (i in id) {
      if (is.null(stations[[i]])) stations[[i]] <- NA_character_
    }
    
    
    if (!all(names(var_id) %in% names(stations))) {
      stop("coding error - seek help from HARSAT team")
    }

  }  
  
  
  stations
}


ctsm_read_contaminants <- function(file, path, info) {
  
  # import functions
  # read in contaminant (and biological effects) data
  
  infile <- file.path(path, file)
  cat("\nReading contaminant and biological effects data from:\n '", 
      infile, "'\n", sep = "")

  if (info$data_format == "ICES_old") {  

    data <- read.csv(
      infile, na.strings = c("", "NULL"), strip.white = TRUE
    )
      
    #data <- read.table(
    #  infile, strip.white = TRUE, sep = "\t", header = TRUE, quote = "\"" , 
    #  na.strings = c("", "NULL"), fileEncoding = "UTF-8-BOM", comment.char = ""  
    #)
    
  } 
  
  if (info$data_format == "ICES_new") {

    var_id <- c(
      "country" = "character",
      "mprog" = "character", 
      "helcom_subbasin" = "character",     
      "helcom_l3" = "character",
      "helcom_l4" = "character",
      "ices_ecoregion" = "character",
      "ospar_region" = "character", 
      "ospar_subregion" = "character", 
      "is_amap_monitoring" = "logical",
      "is_helcom_monitoring" = "logical", 
      "is_medpol_monitoring" = "logical", 
      "is_ospar_monitoring" = "logical",  
      "is_amap_area" = "logical", 
      "is_helcom_area" = "logical",
      "is_ospar_area"= "logical",
      "rlabo" = "character",
      "slabo" = "character",
      "alabo" = "character",               
      "statn" = "character",
      "sd_code_match" = "character",
      "sd_code_name" = "character", 
      "sd_code_replaced" = "character", 
      "sd_name_replaced" = "character", 
      "sd_code_final" = "character",
      "sd_name_final" = "character", 
      "myear" = "integer",
      "date" = "Date",                
      "latitude" = "numeric",
      "longitude" = "numeric",
      "dephu" = "numeric",               
      "dephl" = "numeric",
      "purpm" = "character",
      "finfl" = "character",               
      "param" = "character",
      "pargroup" = "character",
      "matrx" = "character",              
      "basis" = "character",
      "value" = "numeric",
      "munit" = "character",              
      "detli" = "numeric",
      "lmqnt" = "numeric",
      "uncrt" = "numeric",              
      "metcu" = "character",
      "qflag" = "character",
      "vflag" = "character",              
      "metoa" = "character",
      "metcx" = "character",
      "metpt" = "character",
      "metst" = "character",
      "metps" = "character",
      "metfp" = "character",               
      "smtyp" = "character",
      "smpno" = "character",
      "subno" = "character",              
      "dcflgs" = "character",
      "tblAnalysisid" = "integer",
      "tblparamid" = "integer",          
      "tblsampleid" = "character",
      "tblspotid" = "integer",
      "tbluploadid" = "integer"         
    )
    
    required <- names(var_id)
    
  } else if (info$data_format == "external") {

    var_id <- c(
      "country" = "character",
      "station_code" = "character",
      "station_name"= "character",
      "sample_latitude" = "numeric", 
      "sample_longitude" = "numeric", 
      "year" = "integer", 
      "date" = "Date",
      "depth" = "numeric",
      "subseries" = "character",
      "sample" = "character",
      "determinand" = "character",
      "matrix" = "character",
      "basis" = "character",
      "unit" = "character",
      "value" = "numeric",
      "censoring" = "character",
      "limit_detection" = "numeric",
      "limit_quantification" = "numeric", 
      "uncertainty" = "numeric", 
      "unit_uncertainty" = "character",
      "method_pretreatment" = "character",
      "method_analysis" = "character",
      "method_extraction" = "character"
    )

    if (info$compartment == "biota") {
      var_id = c(
        var_id, 
        "species" = "character",
        "sex" = "character",
        "n_individual" = "integer"
      )
    }    

    required <- c(
      "country", "station_code", "station_name", "year", "sample", "determinand",
      "matrix", "unit", "value"
    )

    if (info$compartment %in% c("biota", "sediment")) {
      required <- c(required, "basis")
    }

    if (info$compartment %in% c("biota")) {
      required <- c(required, "species")
    }

  }
  

  if (info$data_format %in% c("ICES_new", "external")) {  

    # check required variables are present in data
    
    data <- read.csv(infile, strip.white = TRUE, nrows = 1)
  
    ok <- required %in% names(data)
    
    if (!all(ok)) {
      id <- required[!ok]
      id <- sort(id)
      stop(
        "The following variables are not in the data file. ", 
        "Please update the data file to continue. ",
        "Note that the variable names are case sensitive.\n",
        "Variables: ", paste(id, collapse = ", ")
      )
    }
    
    
    # read data
    
    ok <- names(var_id) %in% names(data)

    data <- read.csv(
      infile, 
      na.strings = c("", "NULL"),
      strip.white = TRUE,
      colClasses = var_id[ok]
    )
    
  }
  

  # create missing (non-required) variables 
  
  if (info$data_format %in% "external") {
  
    # numeric (non-integer) variables
      
    id <- c(
      "sample_latitude", "sample_longitude", "depth", "limit_detection",
      "limit_quantification", "uncertainty"
    )

    for (i in id) {
      if (is.null(data[[i]])) data[[i]] <- NA_real_
    }
    
    
    # character variables
    
    id <- c(
      "subseries", "basis", "basis", "censoring", "unit_uncertainty", 
      "method_pretreatment", "method_analysis", "method_extraction"
    )
    
    for (i in id) {
      if (is.null(data[[i]])) data[[i]] <- NA_character_
    }
    
    
    # other variables
    
    if (is.null(data$date)) data$date <- as.Date(NA)
    
    if (info$compartment == "biota") {
      if (is.null(data$sex)) data$sex <- NA_character_
      if (is.null(data$n_individual)) data$n_individual <- NA_integer_
    }
    
    if (!all(names(var_id) %in% names(data))) {
      stop("coding error - seek help from HARSAT team")
    }
      
  } else if (info$data_format %in% c("ICES_new", "ICES_old")) {
    
    data$subseries <- NA_character_
    
  }  
  


  # check regional identifiers are in the extraction 

  if (info$data_format == "ICES_new" & info$purpose == "HELCOM") {

    data <- dplyr::rename(
      data,
      HELCOM_subbasin = helcom_subbasin,
      HELCOM_L3 = helcom_l3,
      HELCOM_L4 = helcom_l4,
    )

  }

  if (info$data_format %in% c("ICES_new", "ICES_old")) {
    
    pos <- names(data) %in% info$region_id
    if (sum(pos) != length(info$region_id)) {
      stop("not all regional identifiers are in the data extraction")
    }
  
    names(data)[!pos] <- tolower(names(data)[!pos])     # just in case!
    
  } else if (info$data_format == "external") {
    
    names(data) <- tolower(names(data))  
  }
  
  
  # create more useful names
  # for biota, tblsampleid is the species, tblbioid gives the subsample (individual) 
  # which with matrix gives the unique sample id
  # for sediment, tblsampleid and matrix give the unique sample id
  # to streamline, could make sample = tblbioid (biota) or tblsampleid (sediment)

  
  if (info$data_format == "ICES_new") {
    
    data <- dplyr::rename(
      data,
      submitted.station = statn, 
      sd_name = sd_code_name,
      sd_code = sd_code_match,
      station_name = sd_name_final,
      station_code = sd_code_final,
      method_analysis = metoa,
      method_extraction = metcx,
      method_pretreatment = metpt
    )
    
  } else if (info$data_format == "ICES_old" && info$purpose %in% c("OSPAR", "AMAP")) {
    
    data <- dplyr::rename(
      data,
      offshore = ospar_shore,
      submitted.station = stationname, 
      sd_name = sd_stationname,
      sd_code = sd_stationcode,
      station_name = sd_asmt_stationname,
      station_code = sd_asmt_stationcode,
    )
    
  } else if (info$data_format == "ICES_old" && info$purpose %in% "HELCOM") {
    
    data <- dplyr::rename(
      data,
      submitted.station = statn, 
      sd_name = sd_stationname,
      sd_code = sd_stationcode,
      station_name = sd_replacedby_stationname,
      station_code = sd_replacedby_stationcode,
      method_analysis = metoa
    )

  } 
  
  
  # variables common across purpose and compartments
  
  if(info$data_format %in% c("ICES_old", "ICES_new")) {

    data <- dplyr::rename(
      data,
      year = myear, 
      sample_latitude = latitude,
      sample_longitude = longitude,
      determinand = param, 
      matrix = matrx, 
      unit = munit,
      censoring = qflag,
      qalink = tblanalysisid,
      uncertainty = uncrt,
      unit_uncertainty = metcu, 
      replicate = tblparamid, 
      sample = tblsampleid, 
      limit_detection = detli,
      limit_quantification = lmqnt,
      upload = tbluploadid
    )
    
    # compartment specific variables
    
    var_id <- c("tblbioid", "sexco", "dephu", "noinp")
    replacement <- c("sub.sample", "sex", "depth", "n_individual")
    
    pos <- match(var_id, names(data), nomatch = 0)
    
    if (any(pos > 0)) {
      ok <- pos > 0
      pos <- pos[ok]
      replacement <- replacement[ok]
      names(data)[pos] <- replacement
    }
    
  }
  

  
  # ensure further consistency
  
  data <- dplyr::mutate(
    data, 
    determinand = toupper(.data$determinand),
    unit = tolower(.data$unit)
  )
  
  
  if (info$data_format %in% "ICES_old") {

    id <- c(
      info$region_id, "censoring", "sample", "sub.sample", "sd_code", "station_code", 
      "sd_name", "station_name"
    )
    
    data <- dplyr::mutate(
      data, 
      dplyr::across(any_of(id), as.character)
    )
  
  }
      
  data
}


ctsm_read_QA <- function(file, path, purpose) {
  
  # import functions
  # read in method data (and additional crm information)

  # read in data
  
  infile <- file.path(path, file)
  cat("\nReading QA data from '", infile, "'\n", sep = "")
  
  crm <- read.csv(
    infile, na.strings = c("", "NULL"), strip.white = TRUE
  )

  names(crm) <- tolower(names(crm))

  crm <- dplyr::rename(
    crm,
    crm = crmco,
    basis = crmmb,
    value = crmmv,
    year = myear,
    unit = munit,
    data_type = dtype,
    determinand = param,
    qalink = tblanalysisid, 
    method_analysis = metoa,
    method_extraction = metcx
  )

  crm <- dplyr::mutate(crm, determinand = toupper(.data$determinand))
  
  crm
}

#' Tidy the data
#' 
#' Reduces the size of the extraction by removing redundant variables.
#' Any ad-hoc changes will usually be made between `read_data` and `simplify_data`.
#' The output is in the correct format for `create_timeSeries`.
#'
#' @param ctsm_obj 
#' @export
ctsm_tidy_data <- function(ctsm_obj) {
  
  # import_functions.R
  
  library(tidyverse)
  
  
  info <- ctsm_obj$info
  data <- ctsm_obj$data
  stations <- ctsm_obj$stations
  
  
  # set up oddity directory and back up any previous oddity files
  
  ctsm_initialise_oddities(info$oddity_dir, info$compartment)
  
  
  # tidy station dictionary and contaminant data
  
  stations <- ctsm_tidy_stations(stations, info)

  data <- ctsm_tidy_contaminants(data, info)

  if (info$data_format == "ICES_old") {
    
    QA <- ctsm_obj$QA
    
    wk <- ctsm_tidy_QA(QA, data, info)
    
    data <- wk$data
    QA <- wk$QA
    
    id <- c("method_analysis", "method_extraction") 
    data[id] <- QA[as.character(data$qaID), id]
    
    data$qaID <- NULL
    data$qalink <- NULL

  } 

  
  ctsm_obj$stations <- stations
  ctsm_obj$data <- data
  
  if (info$data_format %in% "ICES_old") ctsm_obj$QA <- NULL
  
  ctsm_obj
}


ctsm_tidy_stations <- function(stations, info) {
  
  cat("\nCleaning station dictionary\n")
  
  # purpose-specific changes for ICES data
  
  if (info$data_format %in% c("ICES_old", "ICES_new")) {
    stations <- do.call(
      paste("ctsm_tidy_stations", info$purpose, sep = "_"),
      list(stations = stations, info = info)
    )
  }
  
  
  # ensure country is consistent
  
  stations <- dplyr::mutate(stations, country = str_to_title(.data$country))
  
  
 # replace backward slash with forward slash in station (long) name

  stations <- dplyr::mutate(
    stations, 
    station_longname = gsub("\\", "/", .data$station_longname, fixed = TRUE)
  )
    

  if (info$data_format %in% c("ICES_old", "ICES_new")) {
      
    # remove stations that have been replaced
    
    stations <- dplyr::filter(stations, is.na(.data$replacedBy))
    

    # check whether any remaining duplicated stations 
    # if present, select most recent of these
    # crude method to sort by station, startYear and endYear (since NAs are sorted last)
    # and then reverse the ordering of the whole data frame so it works with ctsm_check
    
    stations <- tidyr::unite(
      stations, 
      "station_id", 
      all_of(c("country", "station_name")), 
      remove = FALSE
    )
    
    stations <- dplyr::arrange(stations, .data$station_id, .data$startYear, .data$endYear)

    stations <- dplyr::arrange(stations, desc(row_number()))
    
    stations <- ctsm_check(
      stations, 
      station_id, 
      action = "delete.dups",
      message = "Duplicate stations - first one selected",
      file_name = "duplicate_stations", 
      info = info
    )
    
    stations$station_id <- NULL
  }    
  
  # select useful columns
  
  
  col_id <- c(
    info$region_id, "country", "station_name", "station_code", "station_longname", 
    "station_latitude", "station_longitude", "station_type", "waterbody_type"
  )
  
  
  stations <- stations[col_id]
  
  
  # order data
  
  col_id <- c(info$region_id, "country", "station_name")
  
  stations <- dplyr::arrange(stations, dplyr::across(all_of(col_id))) 
  
  stations
}


ctsm_tidy_stations_OSPAR <- function(stations, info) {
  
  # restrict to OSPAR stations used for temporal monitoring 
  # also includes grouped stations: Assessment Grouping for OSPAR MIME
  
  stations <- dplyr::filter(
    stations, 
    grepl("OSPAR", .data$programGovernance) & grepl("T", .data$PURPM)
  )
  
  
  # and stations that are used for contaminants or biological effects
  
  stations <- dplyr::filter(
    stations, 
    switch(
      info$compartment, 
      biota = grepl("CF|EF", .data$dataType), 
      sediment = grepl("CS", .data$dataType), 
      water = grepl("CW", .data$dataType)
    )
  )

  # delete stations outside the OSPAR region
  
  ctsm_check(
    stations, 
    is.na(OSPAR_region) | is.na(OSPAR_subregion), 
    action = "delete", 
    message = "Stations outside OSPAR region", 
    file_name = "stations_outside_area",
    info = info
  )
  
  
  # check all regions are within the appropriate OSPAR_region 
  # can sometimes go wrong due to local shape file errors
  
  if (all(c("OSPAR_region", "OSPAR_subregion") %in% info$region_id)) {
    
    ok <- local({
      id1 <- stations$OSPAR_region
      id2 <- info$region_values[stations$OSPAR_subregion, "OSPAR_region"]
      !is.na(id2) & id1 == id2
    })
    
    ctsm_check(
      stations, !ok, action = "warning", 
      message = "OSPAR subregion in wrong OSPAR region", 
      file_name = "region_information_incorrect", 
      info = info
    )
    
  }    
  
  
  stations
}


ctsm_tidy_stations_HELCOM <- function(stations, info) {
  
  # restrict to stations that have a HELCOM region
  
  stations <- tidyr::drop_na(stations, "HELCOM_subbasin")
  
  
  stations
}


ctsm_tidy_stations_AMAP <- function(stations, info) {
  
  # restrict to AMAP stations used for temporal monitoring 
  # also includes grouped stations: Assessment Grouping for OSPAR MIME
  
  stations <- dplyr::filter(
    stations,
    grepl("AMAP", .data$programGovernance) & grepl("T", .data$PURPM)
  )
  
  # and stations that are used for contaminants or biological effects
  
  stations <- dplyr::filter(
    stations,
    switch(
      info$compartment,
      biota = grepl("CF|EF", .data$dataType),
      sediment = grepl("CS", .data$dataType),
      water = grepl("CW", .data$dataType)
    )
  )
  
  # assume missing OSPAR region information corresponds to AMAP stations in 'region 0'
  
  # add in regional information for East AMAP area
  
  stations <- dplyr::mutate(
    stations, 
    across(c("OSPAR_region", "OSPAR_subregion"), as.character)
  )
  
  stations <- within(stations, {
    OSPAR_region[is.na(OSPAR_region)] <- "0"
    OSPAR_subregion[OSPAR_region %in% "0"] <- "East AMAP area"
  })
  
  
  stations
}


ctsm_tidy_contaminants <- function(data, info) {
  
  cat("\nCleaning contaminant and biological effects data\n")

  
  # check for compatibility between data, max_year, and reporting_window
  
  if (any(data$year > info$max_year)) { 

    message_txt <- paste0(
      "years greater than ", info$max_year, " will be excluded; to change this use", 
      "\n   the max_year argument of ctsm_read_data"
    )

    not_ok <- data$year > info$max_year
      
    data <- ctsm_check(
      data, 
      not_ok, 
      action = "delete", 
      message = message_txt, 
      file_name = "data_too_recent", 
      info = info
    )
    
  }
  
  if (!any(data$year %in% info$recent_years)) {
    
    stop(
      "no data in the period ", min(info$recent_years), " to ", 
      max(info$recent_years), " inclusive; nothing to assess!\n",
      "Consider changing the max_year argument or the reporting_window control\n", 
      "option in ctsm_read_data.",
      call. = FALSE
    )
    
  }
  

  if (info$data_format %in% c("ICES_old", "ICES_new")) {
    
    # check whether submitted station has matched to station correctly for 
    # those countries where extraction is by station
    
    if (info$purpose %in% c("OSPAR", "AMAP")) {
      
      # Denmark, France, Ireland, Norway, Spain (2005 onwards), Sweden, UK
      # need to update this list
      
      odd <- 
        data$country %in% c(
          "Denmark", "France", "Ireland", "Norway", "Sweden", "United Kingdom"
        ) | 
        (data$country %in% "Spain" & data$year > 2005)
      
      # have to use sd_name because station_name also incorporates 
      # replacements and groupings
      
      odd <- odd & !is.na(data$submitted.station) & is.na(data$sd_name)
      
    } else if (info$purpose %in% "HELCOM") {
      
      odd <- data$country %in% c("Denmark", "Sweden") 
      
      odd <- odd & !is.na(data$submitted.station) & is.na(data$station_name)
      
    }    
    
    ctsm_check(
      data, 
      odd, 
      action = "warning", 
      message = "Submitted.station unrecognised by dictionary", 
      file_name = "stations_unrecognised",
      info = info
    )
  }    

  
  # add required variables 'replicate' and 'pargroup' (not present in external data)
  
  if (info$data_format == "external") {
    data$replicate <- seq(from = 1, to = nrow(data), by = 1)
    data$pargroup <- ctsm_get_info(info$determinand, data$determinand, "pargroup")
  }
  
      
  # drop data with no stations  
  
  ok <- !is.na(data$station_name)
  if (!all(ok)) {
    cat("   Dropping data with no stations\n")  
    data <- data[ok, ]
  }
  
  
  # ICES biota data: sample is the species identifier (within a haul say) and 
  # sub.sample is what we usually think of as the sample 
  # swap over and delete sub.sample

  if (info$data_format %in% c("ICES_old", "ICES_new") && info$compartment == "biota") {
    data$sample <- data$sub.sample
    data$sub.sample <- NULL
  }
  
  data$sample <- as.character(data$sample)

  
  # retain useful variables
  
  var_id <- c(
    "station_code", "sample_latitude", "sample_longitude", 
    "year", "date", "time", "depth", 
    "species", "sex", "n_individual", "subseries", "sample", "replicate", 
    "determinand", "pargroup", "matrix", "basis", "filtered", 
    "method_analysis", "method_extraction", "method_pretreatment",
    "unit", "value", "censoring", "limit_detection", "limit_quantification", 
    "uncertainty", "unit_uncertainty", "alabo", "qalink"
  )
  
  data <- data[intersect(var_id, names(data))]
  
  
  data
}



ctsm_tidy_QA <- function(QA, data, info) {
  
  # filter QA by data_type 
  
  id <- switch(info$compartment, biota = "CF", sediment = "CS", water = "CW")
  
  QA <- dplyr::filter(QA, .data$data_type %in% id)
  
  
  # retain useful variables  
  
  QA <- QA[c("alabo", "year", "determinand", "method_extraction", "method_analysis", "qalink")]
  
  
  # relabel organotins (which have been relabelled in adjustment file)
  
  QA <- ctsm_TBT_convert(QA, action = "relabel")
  
  
  # drop crm data with no useful information
  
  QA <- dplyr::filter(QA, !(is.na(.data$method_extraction) & is.na(.data$method_analysis)))
  
  
  # retain unique information
  
  QA <- unique(QA)
  
  
  # qalink and determinand should be unique combinations 
  
  QA <- dplyr::arrange(
    QA, 
    dplyr::across(c("qalink", "determinand", "method_extraction", "method_analysis"))
  )
  
  QA <- ctsm_check(
    QA, 
    paste(qalink, determinand), 
    action = "delete.dups",
    message = "Conflicting QA information - first one selected",
    file_name = "conflicting_QA", 
    info = info
  )
  
  ctsm_link_QA(QA, data, info$compartment)
}


ctsm_link_QA <- function(QA, data, compartment) {
  
  # create unique qaID to link data and QA files 
  # should be able to use qalink and determinand, but alabo and year give
  # useful extra information
  
  var_id <- c("determinand", "qalink", "alabo", "year")
  
  data <- tidyr::unite(
    data, 
    "qaID", 
    all_of(var_id), 
    remove = FALSE
  )
  
  QA <- tidyr::unite(
    QA, 
    "qaID", 
    all_of(var_id), 
    remove = FALSE
  )
  
  
  # restrict QA to those values found in data
  
  QA <- dplyr::filter(QA, .data$qaID %in% data$qaID)
  
  
  # get digestion method based on method_extraction - only needed for sediment
  
  # if (compartment == "sediment") {
  #   QA <- ctsm.link.QA.digestion(QA, compartment)
  # }
  
  
  # tidy up
  
  QA <- tibble::column_to_rownames(QA, "qaID")
  
  list(QA = QA, data = data)
}


#' Create a time series
#' 
#' Cleans the data and turns it into time series structures ready for assessment
#' 
#' @export
ctsm_create_timeSeries <- function(
  ctsm.obj, 
  determinands = ctsm_get_determinands(ctsm.obj$info), 
  determinands.control = NULL, 
  oddity_path = "oddities", 
  return_early = FALSE, 
  print_code_warnings = FALSE, 
  get_basis = get_basis_default,
  normalise = FALSE, 
  normalise.control = list()) {

  # import_functions.R
  
  library(tidyverse)

  # arguments
  
  # determinands: character string of determinands to be assessed; default is to 
  #   take values from determinand reference table
  
  determinands <- force(determinands)
  
  if (!is.character(determinands)) {
    stop("determinands must be a character string")
  }
  
  if (length(determinands) == 0L) {
    stop(
      "No determinands have been been selected for assessment:\n", 
      "  please update the determinand reference table or supply the\n",
      "  determinands directly via the determinands argument."
    )
  }

    
  # normalisation can either be a logical (TRUE uses default normalisation function)
  # or a function
  
  if (length(normalise) != 1L) {
    stop("normalise should be a length 1 logical or a function")
  }
  
  if (!(is.logical(normalise) | is.function(normalise))) {
    stop("normalise should be a length 1 logical or a function")
  }
  

  # get key data structures, and initialise output

  out <- list(call = match.call(), call.data = ctsm.obj$call, info = ctsm.obj$info)
  
  info <- ctsm.obj$info
  station_dictionary <- ctsm.obj$stations
  data <- ctsm.obj$data
  
  rm(ctsm.obj)
 

  is.recent <- function(year) year %in% info$recent_years


  # lots of data cleansing - first ensure oddity directory exists and back up
  # any previous oddity files

  oddity_path <- ctsm_initialise_oddities(info$oddity_dir, info$compartment)


  # checks station dictionary:
  # - no backward slashes in station_name or station_longname
  # - no duplicated or missing station_code 
  
  ctsm_check_stations(station_dictionary)
  
  
  # retains determinands of interest, including auxiliary determinands and those
  # required by determinands.control$variables 
  # checks all determinands of interest are recognised by info$determinand

  wk <- ctsm_check_determinands(info, data, determinands, determinands.control)
  
  data <- wk$data
  determinands <- wk$determinands
  determinands.control <- wk$control

  
  cat("\nCleaning data\n")
  
  id <- c("station_code", "date", "filtered", "species", "sample", "matrix", "determinand")
  id <- intersect(id, names(data))
  ord <- do.call("order", data[id])
  data <- data[ord, ]

  # merge with country and station_name from station dictionary to give more 
  # meaningful output in oddity files

  var_id <- c("country", "station_code", "station_name")
  var_id <- intersect(var_id, names(station_dictionary))

  data <- dplyr::left_join(data, station_dictionary[var_id], by = "station_code")
  data <- dplyr::relocate(data, all_of(var_id))


  # check all stations are in station dictionary 

  ok <- data$station_code %in% station_dictionary$station_code
  data <- ctsm_check(
    data, 
    !ok, 
    action = "delete", 
    message = "Stations in data not in station dictionary", 
    file_name = "unidentified_stations", 
    info = info
  )

  
  # replace synonyms for species by reference_species 
  
  if (info$compartment == "biota") {
    data$species <- ifelse(
      data$species %in% row.names(info$species), 
      info$species[data$species, "reference_species"], 
      data$species
    )
  }  

  
  # drop stations (sediment) or station / species combinations (biota) with no
  # data in the relatively recent period 
  # recent years reduces size of data and removes legacy species not in info files
  
  cat(
    "   Dropping stations with no data between", min(info$recent_years), "and", 
    max(info$recent_years), "\n"
  )

  id.names <- intersect(c("station_code", "species"), names(data))
  id <- do.call("paste", data[id.names])
  ok <- id %in% id[is.recent(data$year)]
  data <- droplevels(data[ok, ])
  

  # drop species that aren't going to be assessed
  # do this after we drop station species combinations because legacy species
  # then don't need to be in the reference tables
  
  if (info$compartment == "biota") {
  
    ok <- ctsm_get_info(info$species, data$species, "assess")
    
    if (!all(ok)) { 
      id <- data$species[!ok]
      id <- unique(id)
      id <- sort(id)
      cat(
        "   Dropping following species - see species reference table:\n", 
        paste(id, collapse = ", "), "\n"
      )
    }
    
    data <- data[ok, ]
  }

  
  # drop data corresponding to stations outside the region (mainly OSPAR or HELCOM requirement)
  
  if (info$all_in_region) {
    
    data <- left_join(
      data, 
      station_dictionary[c("station_code", info$region_id)],
      by = "station_code"
    )
    
    if (any(is.na(data[info$region.id]))) {
      cat("   Dropping stations with missing region information in station dictionary\n")
      data <- tidyr::drop_na(data, all_of(info$region_id))
    }
    
    data <- data[setdiff(names(data), info$region_id)]
  }


  # add variables that are going to be useful throughout
  # pargroup in ICES extraction, but can also be got from info$determinand
  
  data$group <- ctsm_get_info(
    info$determinand, data$determinand, "group", info$compartment, sep = "_"
  )

  if (info$compartment == "biota") {
    data$species_group <- ctsm_get_info(info$species, data$species, "species_group")
  }
  
  if (!"pargroup" %in% names(data)) {
    data$pargroup <- ctsm_get_info(info$determinand, data$determinand, "pargroup")
  }
  

  # drop samples which only have auxiliary data
  
  ok <- with(data, sample %in% sample[group != "Auxiliary"])
  if (any(!ok)) {
    cat("   Dropping samples with only auxiliary variables\n")
    data <- data[ok, ]
  }  

  
  # check variables in the data file have valid values
  # - biota: removes species that are not to be assessed
  # - biota: check species_group, sex and n_individual are appropriate
  # - check all determinands have a valid matrix, basis and unit
  # - check method_analysis (only relevant for bile metabolites)
  # - check value is valid (e.g. > 0 for concentrations)
  
  wk <- c(
    "species_group", "sex", "no_individual", "matrix", "basis", "unit", 
    "method_analysis", "value"
  )
  
  for (var_id in wk) {
    data <- ctsm_check_variable(data, var_id, info)
  }


  # get digestion method for metals in sediments and check that these are valid
  
  if (info$compartment == "sediment") {
    data$digestion <- ctsm_get_digestion(data, info)
  }

  
  # check no replicate measurements (some are genuine replicates, mostly from
  # early years) but will just delete these for simplicity (taking the average
  # has been tried, but gets really messy when there are less thans, different
  # detection limits (sediment), uncertainties etc.)
  
  # not sure of the best place for this
  # should come before the merging with auxiliary variables, otherwise get multiple matches
  # should also come before the linking of e.g. CHR with CHRTR, otherwise the duplicates aren't 
  # necessarily the correct ones
  
  data <- ctsm_check(
    data, paste(sample, determinand, matrix), 
    action = "delete.dups", 
    message = "Replicate measurements, only first retained", 
    file_name = "replicate_measurements", 
    info = info
  )


  # ensure censoring, limit of detection and limit of quantification are consistent

  data <- ctsm_check_censoring(data, info, print_code_warnings)
  

  # convert uncertainty into standard deviations, and remove any associated variables
  
  data <- ctsm_check(
    data, 
    !is.na(uncertainty) & uncertainty <= 0, 
    action = "make.NA", 
    message = "Non-positive uncertainties", 
    file_name = "non_positive_uncertainties", 
    missing_id = "uncertainty",
    info = info
  )
  
  data <- dplyr::mutate(
    data, 
    uncertainty_sd = dplyr::case_when(
      unit_uncertainty %in% "U2" ~ uncertainty / 2, 
      unit_uncertainty %in% "%" ~ value * uncertainty / 100, 
      TRUE ~ uncertainty
    ),
    uncertainty_rel = 100 * (uncertainty_sd / value)
  )                 

  wk_id <- match("unit_uncertainty", names(data))
  wk_n <- ncol(data)
  data <- data[c(
    names(data)[1:wk_id], 
    "uncertainty_sd", "uncertainty_rel", 
    names(data)[(wk_id+1):(wk_n-2)])]
  
  ctsm_check(
    data, 
    !is.na(uncertainty) & uncertainty_rel >= 100, 
    action = "warning", 
    message = "Large uncertainties", 
    file_name = "large uncertainties",
    info = info
  )
  
  # delete data with large relative uncertainties
  
  data <- dplyr::mutate(
    data, 
    uncertainty_sd = dplyr::if_else(
      .data$uncertainty_rel < 100, 
      .data$uncertainty_sd, 
      NA_real_
    ),
    uncertainty = .data$uncertainty_sd, 
    unit_uncertainty = NULL,
    uncertainty_sd = NULL, 
    uncertainty_rel = NULL
  )
  

  # sort out determinands where several determinands represent the same variable of interest
  # three types of behaviour: replace, sum and bespoke

  for (i in names(determinands.control)) {
    
    wk <- determinands.control[[i]]
    linkFunction <- switch(
      wk$action, 
      replace = determinand.link.replace,
      sum = determinand.link.sum,
      bespoke = get(paste("determinand.link", i, sep = "."), mode = "function")
    )
    
    data <- do.call(
      linkFunction, 
      list(data = data, keep = i, drop = wk$det, info = info)
    )
  }  

  # drop any remaining unwanted determinands (from sum and perhaps bespoke functions);
  # could make this more elegant!
  
  id <- c(determinands, ctsm_get_auxiliary(determinands, info))
  
  data <- dplyr::filter(data, .data$determinand %in% id)


  # remove 'extra' data from time series which have a subseries classification 
  # for some but not all records

  data <- ctsm_check_subseries(data)
    

  # set rownames to NULL(ie. back to auto numeric)
  
  rownames(data) <- NULL

  cat("\nCreating time series data\n")  

  data <- data[setdiff(names(data), c("qalink", "alabo"))]


  # create new.unit and concentration columns comprising the details from the
  # determinand file in the information folder, required to get correct unit details
  
  data$new.unit <- ctsm_get_info(
    info$determinand, data$determinand, "unit", info$compartment, sep = "_"
  )
  
  data$concentration <- data$value
  

  # convert data to conventional units - see "information determinands.csv"
  # NB value is in original units, concentration is in converted units
  # this applies to all determinands (including auxiliary variables)
  # can't change bases here because need to merge with auxiliary variables

  id <- c("concentration", "uncertainty", "limit_detection", "limit_quantification")

  data[id] <- lapply(data[id], convert_units, from = data$unit, to = data$new.unit) 


  # merge auxiliary data with determinand data

  data <- ctsm_merge_auxiliary(data, info)
  

  # impute %femalepop when missing and sex = 1 - write out remaining
  # missing values for correction
  
  if (info$compartment == "biota") {
    data <- ctsm.imposex.check.femalepop(data)
  }
  

  # convert data to basis of assessment
  
  data <- ctsm_convert_to_target_basis(data, info, get_basis)
  

  if (return_early) {
    out  = c(
      out, 
      ctsm.import.value(
        data, 
        station_dictionary, 
        info$compartment, 
        info$purpose, 
        print_code_warnings)
    )
    
    return(out)
  }


  # estimate missing uncertainties

  data$uncertainty <- ctsm_estimate_uncertainty(data, "concentration", info)

  if (info$compartment == "sediment") {
    
    for (norm_id in c("AL", "LI", "CORG", "LOIGN")) {
      
      if (norm_id %in% names(data)) {
        norm_uncrt <- paste0(norm_id, ".uncertainty")
        data[[norm_uncrt]] <- ctsm_estimate_uncertainty(data, norm_id, info)
      }    
    }
    
  }

  
  # filter contaminant data to remove bivalve and gastropod records in the 
  # spawning season when they are elevated / more variable
  
  # this is placed here for convenience at present
  # it should be in tidy_contaminants, where other purpose specific filtering 
  # operations will be done (when code for matching new ICES data with stations
  # is available) 
  # however, first need to make pargroup information available for all 
  # determinands that might be downloaded, not just those that are of interest

  if (info$compartment == "biota" && !is.null(info$bivalve_spawning_season)) {

    data <- dplyr::mutate(
      data, 
      .month = months(as.Date(.data$date)),
      .not_ok = species_group %in% c("Bivalve", "Gastropod") & 
        (is.na(.month) | .month %in% info$bivalve_spawning_season) & 
        !group %in% c("Effects", "Imposex", "Metabolites") 
    )
    
    if (any(data$.not_ok)) {
      cat(
        "   Dropping bivalve and gastropod contaminant data collected during the\n", 
        "   spawning season, which is taken to be the following months:\n   ",
        paste(info$bivalve_spawning_season, collapse = ", "), 
        "\n"
      )
      data <- dplyr::filter(data, !.not_ok)
    }
    
    data[c(".month", ".not_ok")] <- NULL
  }


  # normalise data
  # e.g. sediment data to 5% aluminium and 2.5% organic carbon (OSPAR)
  # biota data to 5% lipid (HELCOM)

  if (is.logical(normalise) && normalise) {
    data <- switch(
      info$compartment,
      biota = 
        ctsm_normalise_biota(data, station_dictionary, info, normalise.control),
      sediment = 
        ctsm_normalise_sediment(data, station_dictionary, info, normalise.control),
      water = stop("there is no default normalisation function for water")
    )
  } else if (is.function(normalise)) {
    data <- normalise(data, station_dictionary, info, normalise.control)
  }
    

  # remove concentrations where cv of uncertainty > 100%
  # tidy up missing data for uncertainties and censoring
  
  data <- within(data, {
    concentration[uncertainty > concentration] <- NA
    uncertainty[is.na(concentration)] <- NA
    censoring[is.na(concentration)] <- NA
  })
    

  # drop groups of data at stations with no data in recent years

  cat("   Dropping groups of compounds / stations with no data between", 
      min(info$recent_years), "and", max(info$recent_years), "\n")
  id_names <- intersect(c("station_code", "filtered", "species", "group"), names(data))
  id <- do.call("paste", data[id_names])
  ok <- id %in% id[is.recent(data$year) & !is.na(data$concentration)]
  data <- data[ok, ]
  

  # create seriesID and timeSeries object
  
  out <- c(
    out, 
    ctsm_import_value(data, station_dictionary, info)
  )
  
  out
}


# default routine to check (and clean) data attributes when creating timeSeries

ctsm_check <- function(
  data, id, 
  action = c("delete", "warning", "make.NA", "delete.dups"), 
  message, file_name, info, missing_id) {

  # location: import_functions.R
  # purpose: general checking routine 
  
  # action options: 
  # delete - delete oddities (in id) 
  # warning - print out oddities, but don't change data 
  # make.NA - replace values of missingID corresponding to oddities with NA 
  # delete.dups - delete any duplicated values identified in id, but print out 
  #   all instances of duplication to see what is going on
  
  action <- match.arg(action)

  # evaluate logical identifying oddities:
  # if action = delete.dup find all duplicated values in id
  # else evaluate id to obtain logical of oddities

  if (action == "delete.dups") {
    duplicateVar <- eval(substitute(id), data, parent.frame())
    odd_id <- duplicated(duplicateVar) | duplicated(duplicateVar, fromLast = TRUE)
  } else {
    odd_id <- eval(substitute(id), data, parent.frame())
    if (!is.logical(odd_id)) stop("'id' must be logical")
    odd_id <- odd_id & !is.na(odd_id)
  }

  
  # do nothing if no oddities
  
  if (!any(odd_id)) return(data) 


  # get oddities and output results
  
  oddities <- data[odd_id, ]


  # output results 
  
  outfile_name <- paste0(file_name, ".csv")
  outfile_name <- gsub(" ", "_", outfile_name)
  outfile <- file.path(info$oddity_dir, info$compartment, outfile_name)
  
  txt <- switch(action, delete = ": deleted data in '", ": see ")
  message("   ", message, txt, outfile_name) 

  readr::write_excel_csv(oddities, outfile, na = "")
  
  if (action == "warning") return()
  
  data <- switch(
    action, 
    delete = data[!odd_id, ],
    make.NA = {data[odd_id, missing_id] <- NA ; data},
    delete.dups = data[!duplicated(duplicateVar), ]
  )

  data <- droplevels(data)
    
  data        
}
  

ctsm_import_value <- function(data, station_dictionary, info) {
  
  # import_functions.R
  
  # identifies timeseries in the data and creates timeSeries metadata object


  # order data 

  id = c(
    "station_code", "species", "filtered", "year", "sex", "sample", "group", 
    "determinand"
  )
  
  data <- dplyr::arrange(data, dplyr::across(any_of(id)))


  # select variables of interest

  id <- c(
    "station_code", "sample_latitude", "sample_longitude", "filtered", 
    "species", "sex", "depth",
    "year", "date", "time", "sample",   
    "matrix", "subseries", "group", "determinand", "basis", "unit", "value", 
    "method_analysis", "n_individual", 
    "concOriginal", "censoringOriginal", "uncrtOriginal", 
    "concentration", "new.basis", "new.unit", "censoring",  
    "limit_detection", "limit_quantification", "uncertainty"
  )
  
  auxiliary <- ctsm_get_auxiliary(data$determinand, info)
  auxiliary_id <- paste0(
    rep(auxiliary, each = 5), 
    c("", ".censoring", ".limit_detection", ".limit_quantification", ".uncertainty") 
  )
    
  id <- c(id, auxiliary_id)

  data <- dplyr::select(data, any_of(id))

  row.names(data) <- NULL  
  data <- droplevels(data)
  
  
  # get seriesID and timeSeries structure 
  
  id <- c("station_code", "determinand")
  
  if (info$compartment %in% c("biota")) {
    id <- c(id, "species", "matrix", "subseries", "sex", "method_analysis")
  }

  if (info$compartment %in% "sediment") {
    id <- c(id, "matrix", "subseries")
  }
  
  if (info$compartment %in% "water") {
    id <- c(id, "filtered", "subseries") 
  }
  
  timeSeries <- data[id]
  

  # only retain sex and method_analysis information in timeSeries where 
  # necessary to distinguish different series
  
  # currently hard-wired for EROD and Metabolites in biota: issue raised as 
  # future enhancement 

  if (info$compartment == "biota") {
    timeSeries <- dplyr::mutate(
      timeSeries,
      sex = dplyr::if_else(.data$determinand %in% "EROD", .data$sex, NA_character_),
      .group = ctsm_get_info(
        info$determinand, .data$determinand, "group", "biota", sep = "_"
      ),
      method_analysis = dplyr::if_else(
        .group %in% "Metabolites", 
        .data$method_analysis, 
        NA_character_
      ),
      .group = NULL,
    )
  }
    
 
  # create seriesID column in data, filled with the concatenated values of timeSeries and 
  # reorder variables

  timeSeries <- tidyr::unite(
    timeSeries, 
    "seriesID", 
    all_of(names(timeSeries)), 
    sep = " ", 
    remove = FALSE, 
    na.rm = TRUE
  )
 
  data <- cbind(seriesID = timeSeries$seriesID, data)

  
  # pick up new.basis and new.unit for each timeSeries 
  
  timeSeries$basis <- data$new.basis
  timeSeries$unit <- data$new.unit
  
  data$new.basis <- data$new.unit <- NULL
  

  # timeSeries is now the unique rows of timeSeries
  
  timeSeries <- droplevels(dplyr::distinct(timeSeries))
  
  # id <- setdiff(names(timeSeries), c("basis", "unit"))
  
  timeSeries <- column_to_rownames(timeSeries, "seriesID")

  
  # change timeSeries output columns to fit the levels of the xml requirements
  # this is a legacy requirement and an issue has been raised to fix this
  
  # timeSeries <- changeToLevelsForXML(timeSeries, info)

  
  # check no replicate measurements within time series
  
  data <- ctsm_check(
    data, 
    paste(seriesID, sample), 
    action = "delete.dups",
    message = "Measurements still replicated, only first retained",
    file_name = "replicate_measurements_extra",
    info = info
  )
  

  # remove any series that don't have data in recent years
  
  is_recent_data <- (data$year %in% info$recent_years) & !is.na(data$concentration)
  
  id <- data$seriesID[is_recent_data]
  id <- unique(id)
  
  data <- data[data$seriesID %in% id, ]
  data <- droplevels(data)
  
  ok <- row.names(timeSeries) %in% id
  timeSeries <- timeSeries[ok, ]
  timeSeries <- droplevels(timeSeries)
  

  # remove unused stations and rownames from the station dictionary

  ok <- station_dictionary$station_code %in% timeSeries$station_code 
  
  station_dictionary <- station_dictionary[ok, ]
  
  row.names(station_dictionary) <- NULL
  
  
  
  list(data = data, stations = station_dictionary, timeSeries = timeSeries)
}




ctsm_check_stations <- function(stations) {

  # import_functions.R
    
  # ensure no backward slashes in station dictionary
  
  not_ok1 <- grepl("\\", stations$station_name, fixed = TRUE) 
  not_ok2 <- grepl("\\", stations$station_longname, fixed = TRUE)
  
  not_ok1 <- any(not_ok1)
  not_ok2 <- any(not_ok2)
  
  if (not_ok1 | not_ok2) {
    var_id <- c("station_name", "station_longname")
    var_id <- var_id[c(not_ok1, not_ok2)]
    stop(
      "backward slashes not allowed in station dictionary variables: ", 
      toString(var_id)
    )
  }

  
  # ensure no duplicated station codes
  
  if (any(duplicated(stations$station_code))) {
    stop("duplicated values of station_code not allowed in station dictionary")
  }


  invisible()  
}



ctsm_get_digestion <- function(data, info) {

  # import_functions.R
  # get digestion method for metals based on method_extraction
  
  data$digestion <- rep(NA_character_)


  # only apply ctsm_get_info on metals data to avoid all the potential issues
  # associated with organics submissions
  
  is_metal <- data$pargroup %in% "I-MET"  
  
  data[is_metal, "digestion"] <- ctsm_get_info(
    info$method_extraction,
    data[is_metal, "method_extraction"], 
    "digestion", 
    na_action = "input_ok"
  )
  

  # trap mercury digestion wth method_extraction of NON  
  
  id <- data$method_analysis %in% "AAS-CV" & data$method_extraction %in% "NON"
  if (any(id)) {
    message("   Warning: Ad-hoc way of dealing with mercury digestion that has method_extraction of NON")
    data[id, "digestion"] <- "Tot"
  }
  

  # convert Pe digestion to Ps due to lack of pivot values
  
  id <- data$digestion %in% "Pe"
  if (any(id)) {
    message("   Warning: Pe digestion converted to Ps due to lack of pivot values\n")
    data[id, "digestion"] <- "Ps"
  }
  

  # check for missing values, organics analysed by a digestion method, and CORG
  # analysed by a digestion method other than HCL

  not_ok <- is_metal & is.na(data$digestion) 

  ctsm_check(
    data, 
    not_ok, 
    action = "warning", 
    message = "Missing digestion methods", 
    file_name = "digestion_errors", 
    info = info
  )

  data$digestion  
}


ctsm_initialise_oddities <- function(path, compartment) {

  # location: import_functions.R
  # purpose: sets up oddity folder and backs up previous runs 
  
  # create oddity directories and backup if necessary
  
  if (!dir.exists(path)) dir.create(path)

  output <- file.path(path, compartment)

  if (!dir.exists(output)) {
    dir.create(output)
    cat("\nOddities will be written to '", output, "'\n", sep = "")
    return(output)
  } 

  backup <- file.path(path, paste(compartment, "backup", sep = "_"))
  
  if (!dir.exists(backup)) dir.create(backup) 
      
  cat("\nOddities will be written to '", output, "' with previous oddities ", 
      "backed up to\n '", backup, "'\n", sep = "")
  
  old.files <- dir(output, full.names = TRUE)
  file.copy(from = old.files, to = backup, overwrite = TRUE)
  file.remove(old.files)
    
  invisible()  
}


ctsm.imposex.check.femalepop <- function(data, info) {

  # makes missing values == 100% if sex == F
  
  if (!any(data$group %in% "Imposex"))
    return(data)

  impID <- data$group %in% "Imposex"
  
  replaceID <- impID & is.na(data[["%FEMALEPOP"]]) & data$sex == "F"
  data[replaceID, "%FEMALEPOP"] <- 100

  missingID <- impID & is.na(data[["%FEMALEPOP"]])
  ctsm_check(
    data[impID, ], 
    missingID, 
    action = "warning", 
    message = "Missing FEMALEPOP values", 
    file_name = "imposex_missing_FEMALEPOP",
    info = info
  )
  
  return (data)
}


ctsm_check_determinands <- function(info, data, determinands, control = NULL) {

  # checks all determinands are recognised in info files
  # checks determinands are not also in control (if they are to be replaced)
  # reduces data file and determinand structures so that they only contain required values
  
  # utility function to get all determinand names from control structure
  
  get_control_dets <- function(control, .names = TRUE) {  
    if (is.null(control)) 
      return(NULL)
    
    out <- lapply(control, "[[", "det")
    out <- unlist(out)
    out <- unname(out)
    
    if (.names) {
      return(c(names(control), out))
    } else {
      return(out)
    }
  }
  
  
  # check control values to be replaced are not also in determinands
  
  if (!is.null(control)) {
    lapply(control, function(ls) {
      not_ok <- ls$det %in% determinands & ls$action %in% "replace"
      if (any(not_ok)) {
        stop ("Replaced determinands in 'determinands.control' are also in 'determinands': ", 
              paste(ls$det[not_ok], collapse = ", ")
        )
      }
    })  
  }
  

  # check all determinands are recognised by info$determinands
  # note all auxiliaries are guaranteed to be info$determinands
  
  id <- c(determinands, get_control_dets(control))

  ctsm_check_reference_table(id, info$determinand, "determinand")  
  
    
  # simplify determinands and determinand.control so they only contain values that are
  # required
  
  id <- unique(data$determinand)
  
  if (!is.null(control))
    id <- c(id, names(control))
  
  determinands <- intersect(determinands, id)

  
  auxiliary <- ctsm_get_auxiliary(determinands, info)
  
  if (!is.null(control)) {
    ok <- names(control) %in% c(determinands, auxiliary, get_control_dets(control, .names = FALSE))
    
    if (!any(ok)) {
      control <- NULL
    } else {
      control <- control[ok]
    }
  }
    

  # only retain required determinands and auxiliary variables 

  id <- c(determinands, get_control_dets(control), auxiliary)

  ok <- data$determinand %in% id
  
  data <- data[ok, ]
  
  list(data = data, determinands = determinands, control = control)
}



determinand.link.check <- function(data, keep, drop, printDuplicates = TRUE, ...) {

  # check whether any drop and keep are both submitted for the same sample and 
  # matrix and, if so, delete drop - note that ctsm_check doesn't do the
  # deleting because it isn't necessarily the first instance that is retained
  
  ID <- with(data, paste(sample, matrix))
  
  dropID <- data$determinand %in% drop 
  keepID <- data$determinand %in% keep
  
  dups <- ID %in% intersect(ID[dropID], ID[keepID])

  if (printDuplicates) {
    
    dropTxt <- paste(drop, collapse = ", ")
    
    dupsID <- dups & (dropID | keepID)
    
    ctsm_check(
      data, 
      dupsID, 
      action = "warning",  
      message = paste(
        keep, "and", dropTxt, "submitted in same sample - deleting", dropTxt, 
        "data"
      ), 
      file_name = paste("determinand_link", keep, sep = "_"), 
      ...
    )
  }
  
  data[!(dups & dropID), ]
}  
  

determinand.link.replace <- function(data, keep, drop, ...) {

  # core function for relabelling determinand 'drop' as determinand 'keep'
  # most of the work is checking that there aren't data submitted as both for the same
  # sample
  
  stopifnot(length(keep) == 1, length(drop) == 1)
  
  if (any(data$determinand %in% drop)) 
    cat("   Data submitted as", drop, "relabelled as", keep, "\n")
  else return(data)
  
  
  # check for samples with both drop and keep and, if they exist, delete drop

  data <- determinand.link.check(data, keep, drop, ...)
  
  
  # relabel the levels so that drop becomes keep
  
  id <- data$determinand %in% drop
  data$determinand[id] <- keep
  
  data
}  


determinand.link.imposex <- function(data, keep, drop, ...) {
  
  stopifnot(length(keep) == 1, length(drop) == 1)

  detID <- c(keep, drop)
  
  # for imposex, indices and stages aren't linked by sample, but by visit
  # will assume, for simplicity, that only one visit per year
  
  visitID <- with(data, paste(station_code, year))
  
  
  # find visits when both indivuduals and stages reported and check consistent
  
  ok <- by(data, visitID, function(x) {
    with(x, {
      if (! all(detID %in% determinand)) return(TRUE)
      sumKeep <- sum(value[determinand == keep]) 
      sumDrop <- mean(value[determinand == drop]) * length(value[determinand == keep])
      abs(sumKeep - sumDrop) < 0.5
    })
  })

  dups <- visitID %in% names(ok)[!ok] & data$determinand %in% detID
  
  ctsm_check(
    data, 
    dups, 
    action = "warning",  
    message = paste("inconsistent", keep, "and", drop, "submitted in same year"), 
    file_name = paste("determinand_link", keep, sep = "_"), 
    ...
  )
  
  # delete indices submitted in same visit as individual data (whether consistent or not)

  dups <- tapply(data$determinand, visitID, function(x) all(detID %in% x))
  dups <- visitID %in% names(dups)[dups]
  
  data <- data[!(dups & data$determinand %in% drop), ]
    
  
  # relabel any remaining indices as stage data
    
  if (any(data$determinand %in% drop)) 
    message("   Data submitted as ", drop, "relabelled as ", keep)
  else return(data)
  
  id <- data$determinand %in% drop
  data$determinand[id] <- keep

  data
}  

determinand.link.VDS <- determinand.link.IMPS <- determinand.link.INTS <- determinand.link.imposex

determinand.link.BBKF <- function(data, keep, drop, ...) {
  
  stopifnot(
    identical(keep, "BBKF"), 
    identical(sort(drop), c("BBF", "BBJF", "BBJKF", "BKF"))
  )
  
  # first sum samples with both BBF and BKF
  
  data <- determinand.link.sum(data, "BBKF", c("BBF", "BKF"))
  
  # now sum samples with both BBJF and BKF to give BBJKF
  
  data <- determinand.link.sum(data, "BBJKF", c("BBJF", "BKF"))
  
  # now replace BBJKF with BBKF
  
  data <- determinand.link.replace(data, "BBKF", "BBJKF")
  
  data
}



assign("determinand.link.LIPIDWT%", function(data, keep, drop, ...) {

  stopifnot(identical(keep, "LIPIDWT%"), identical(sort(drop), c("EXLIP%", "FATWT%")))

  # if multiple values present, choose FATWT%, then LIPIDWT%, then EXLIP% (from Foppe)
  
  data <- determinand.link.check(data, keep = "LIPIDWT%", drop = "EXLIP%", printDuplicates = FALSE, ...)
  data <- determinand.link.check(data, keep = "FATWT%", drop = "EXLIP%", printDuplicates = FALSE, ...)
  data <- determinand.link.check(data, keep = "FATWT%", drop = "LIPIDWT%", printDuplicates = FALSE, ...)

  if (!any(data$determinand %in% drop)) return(data)

  # relabel the levels so that drop becomes keep
  
  cat("   Data submitted as EXLIP% or FATWT% relabelled as LIPIDWT%", "\n")

  id <- data$determinand %in% drop
  data$determinand[id] <- keep
  
  data
})  


determinand.link.sum <- function(data, keep, drop, ...) {
  
  stopifnot(length(keep) == 1, length(drop) > 1)
  
  if (!any(data$determinand %in% drop)) 
    return(data)


  # identify samples with drop and not keep, which are the ones that will be summed
  # if keep already exists, then don't need to do anything
  # don't delete drop data because might want to assess them individually
  
  ID <- with(data, paste(sample, matrix))
  
  dropID <- data$determinand %in% drop 
  keepID <- data$determinand %in% keep
  
  sum_ID <- ID %in% setdiff(ID[dropID], ID[keepID])

  if (length(sum_ID) == 0)
    return(data)
  
  
  dropTxt <- paste(drop, collapse = ", ")
  cat("   Data submitted as", dropTxt, "summed to give", keep, fill = TRUE)


  # get relevant sample matrix combinations
  
  data <- split(data, with(data, determinand %in% drop & sum_ID))
  
  ID <- with(data[["TRUE"]], paste(sample, matrix))

  summed_data <- by(data[["TRUE"]], ID, function(x) {
    
    # check all bases are the same 

    stopifnot(n_distinct(x$basis) == 1)
    
    if (!all(drop %in% x$determinand)) return(NULL)

    
    # adjust values if units vary
    # ideally use unit in info$determinand, but makes it more awkward because
    # have to pass in compartment
    
    if (n_distinct(x$unit) > 1) {

      # get modal value of unit
      
      unit_values <- unique(x$unit)
      target_unit <- unit_values[which.max(tabulate(match(x$unit, unit_values)))]
           
      id <- c("value", "uncertainty", "limit_detection", "limit_quantification")
      
      x[id] <- lapply(x[id], convert_units, from = x$unit, to = target_unit)

      x$unit <- target_unit
    }      
    
    
    # make output row have all the information from the largest determinand (ad-hoc) 
    # ensures a sensible qaID, method_analysis, etc.
    
    out <- x[which.max(x$value), ]
    
    out$determinand <- keep
    
    # sum value and limit_detection, make it a less-than if all are less-thans, and take 
    # proportional uncertainty from maximum value (for which uncertainty is reported)
    
    out$value <- sum(x$value)
    out$limit_detection <- sum(x$limit_detection)
    out$limit_quantification <- sum(x$limit_quantification)
    
    if ("" %in% x$censoring)
      out$censoring <- ""
    else if (n_distinct(x$censoring) == 1) 
      out$censoring <- unique(x$censoring) 
    else 
      out$censoring <- "<"

    if (all(is.na(x$uncertainty))) 
      out$uncertainty <- NA
    else {
      wk <- x[!is.na(x$uncertainty), ]
      pos <- which.max(wk$value)
      upct <- with(wk, uncertainty / value)[pos]
      out$uncertainty <- out$value * upct
    }

    out
    
  })

  summed_data <- do.call(rbind, summed_data)

  
  # see how many samples have been lost due to incomplete submissions
  # need a trap for no summed_data

  nTotal <- length(unique(ID))
  nSummed <- if (is.null(summed_data)) 0 else nrow(summed_data)
  nLost <- nTotal - nSummed
  
  if (nLost > 0) 
    message("     ", nLost, " of ", nTotal, " samples lost due to incomplete submissions")

  
  # combine data for both drop and keep and then add back into main data set
  
  data[["TRUE"]] <- rbind(data[["TRUE"]], summed_data)
    
  data <- do.call(rbind, data)
  
  data
}  




determinand.link.TEQDFP <- function(data, keep, drop, ...) {
  
  stopifnot(length(keep) == 1, length(drop) > 1)
  
  if (!any(data$determinand %in% drop)) 
    return(data)
  
  
  # identify samples with drop and not keep, which are the ones that will be summed
  # if keep already exists, then don't need to do anything
  # don't delete drop data because might want to assess them individually
  
  ID <- with(data, paste(sample, matrix))
  
  dropID <- data$determinand %in% drop 
  keepID <- data$determinand %in% keep
  
  sum_ID <- ID %in% setdiff(ID[dropID], ID[keepID])
  
  if (length(sum_ID) == 0)
    return(data)
  
  dropTxt <- paste(drop, collapse = ", ")
  cat("   Data submitted as", dropTxt, "summed to give", keep, fill = TRUE)
  
  
  # get relevant sample matrix combinations
  
  data <- split(data, with(data, determinand %in% drop))
  
  ID <- with(data[["TRUE"]], paste(sample, matrix))
  
  summed_data <- by(data[["TRUE"]], ID, function(x) {
    
    # only some of the determinands are mandatory - otherwise we woudld lose everything 
    # mandatory determinands contribute at least 1% to the total TEQ based on a quick look-see!
    # the order below is based on % contribution

    # mandatory <- c(
    #   "CB126", "CDF2N", "CDD1N", "CDF2T", "TCDD", "CB169", "CB118", "CDFP2", "CDD6X", "CDF4X", 
    #   "CDF6X")
    # 
    # if (!all(mandatory %in% x$determinand)) 
    #   return(NULL)
    

    # check all bases are the same 
    
    if (!all(drop %in% x$determinand)) return(NULL)

    stopifnot(n_distinct(x$basis) == 1)
    
    
    # convert to ug/kg and then to TEQ
    
    id <- c("value", "uncertainty", "limit_detection", "limit_quantification")
    
    x[id] <- lapply(x[id], convert_units, from = x$unit, to = "ug/kg")
    
    TEQ <- info_TEQ[as.character(x$determinand)]
    
    x[id] <- lapply(x[id], "*", TEQ)
    
    
    # make output row have all the information from the largest determinand (ad-hoc) 

    # ensures a sensible qaID, method_analysis, etc.

    
    out <- x[which.max(x$value), ]
    
    out$determinand <- keep
    out$unit <- "TEQ ug/kg"
    out$group <- "Dioxins"
    out$pargroup <- "OC-DX"
    
    # sum value and limit_detection, make it a less-than if all are less-thans, and take 
    # proportional uncertainty from maximum value (for which uncertainty is reported)
    # if no uncertainties reported at all, then have provided value of CB126 in info.unertainty
    # with sdConstant multiplied by 0.1 to reflect TEQ effect on detection limit
    
    out$value <- sum(x$value)
    out$limit_detection <- sum(x$limit_detection)
    out$limit_quantification <- sum(x$limit_quantification)
    
    if ("" %in% x$censoring)
      out$censoring <- ""
    else if (n_distinct(x$censoring) == 1) 
      out$censoring <- unique(x$censoring) 
    else 
      out$censoring <- "<"
    out$censoring <- if(all(x$censoring %in% "<")) "<" else ""
    
    if (all(is.na(x$uncertainty))) 
      out$uncertainty <- NA
    else {
      wk <- x[!is.na(x$uncertainty), ]
      pos <- which.max(wk$value)
      upct <- with(wk, uncertainty / value)[pos]
      out$uncertainty <- out$value * upct
    }
    
    out
    
  })
  
  summed_data <- do.call(rbind, summed_data)
  
  
  # see how many samples have been lost due to incomplete submissions
  
  nTotal <- length(unique(ID))
  nLost <- length(unique(ID)) - nrow(summed_data)
  if (nLost > 0) 
    message("     ", nLost, " of ", nTotal, " samples lost due to incomplete submissions")

  
  # combine data for both drop and keep and then add back into main data set
  
  data[["TRUE"]] <- rbind(data[["TRUE"]], summed_data)

  data <- do.call(rbind, data)
  
  data
}  


ctsm_check_censoring <- function(data, info, print_code_warnings) {
  
  # location: import_functions.R
  # purpose: checks that censoring, limit_detection and limit_quantification are 
  #   consistent and makes sensible adjustments to censoring
  
  # define function for testing equality - also need a relative component to cope with some tiny
  # value submitted with units g/g
  
  my_near <- function(x, y) {
    abs <- dplyr::near(x, y) 
    rel <- dplyr::case_when(
      x > 0   ~ abs((x - y) / y) < 1e-5,
      TRUE    ~ TRUE
    )
    abs & rel
  }
  
  
  # check detection limits associated with data are positive 
  
  if (print_code_warnings) {
    warning(
      'Need to make checking of detection limits determinand specific', 
      call. = FALSE
    )
  }  
  
  data <- ctsm_check(
    data, 
    limit_detection <= 0, 
    action = "make.NA", 
    message = "Non-positive detection limits", 
    file_name = "non_positive_det_limits", 
    missing_id = "limit_detection",
    info = info
  )
  
  data <- ctsm_check(
    data, 
    limit_quantification <= 0, 
    action = "make.NA", 
    message = "Non-positive quantification limits", 
    file_name = "non_positive_quant_limits", 
    missing_id = "limit_quantification",
    info = info
  )


  # limit_quantification must be greater than limit_detection - 
  # otherwise set both to missing
  
  data <- ctsm_check(
    data,  
    limit_quantification <= limit_detection, 
    action = "make.NA", 
    message = "Limit of quantification less than limit of detection", 
    file_name = "limits_inconsistent", 
    missing_id = c("limit_detection", "limit_quantification"), 
    info = info
  )


  # check for valid values of censoring
  # censoring cannot contain a > (although possible for some biological effects - 
  # need to revisit)
  # also check for other unrecognised characters
  
  if (print_code_warnings) {
    warning(
      "Need to make censoring determinand specific when biological effects are introduced", 
      call. = FALSE
    )
  }

  data <- within(data, {
    levels(censoring) <- c(levels(censoring), "")
    censoring[is.na(censoring)] <- ""
    censoring <- dplyr::recode(censoring, "<~D" = "D", "D~<" = "D", "<~Q" = "Q", "Q~<" = "Q")
  })

  data <- ctsm_check(
    data, 
    !censoring %in% c("", "D", "Q", "<"), 
    action = "delete", 
    message = "Unrecognised censoring values", 
    file_name = "censoring_codes_unrecognised", 
    info = info
  )


  # if censoring = D, then value must equal detection_limit
  # if censoring = Q, then value must equal limit_quantification
  
  id <- with(data, {
    out1 <- censoring == "D" & 
      (is.na(limit_detection) | !my_near(value, limit_detection))
    out2 <- censoring == "Q" & 
      (is.na(limit_quantification) | !my_near(value, limit_quantification))
    out1 | out2
  })

  ctsm_check(
    data, 
    id, 
    action = "warning",
    message = "Censoring codes D and Q inconsistent with respective limits",
    file_name = "censoring_codes_inconsistent",
    info = info
  )


  # resolve these inconsistencies
  
  data <- within(data, {
    censoring <- dplyr::case_when(
      censoring %in% "" ~ "",
      !is.na(limit_detection) & my_near(value, limit_detection) ~ "D",
      !is.na(limit_quantification) & my_near(value, limit_quantification) ~ "Q",
      TRUE ~ "<"
    )
  })    


  # check limit_detection is less than (or equal to) concentration
  # NB would need to be revised if limits are submitted in different units to concentrations
  
  data <- ctsm_check(
    data, 
    id = censoring %in% c("", "<") & limit_detection > value, 
    action = "make.NA", 
    message = "Detection limit higher than data", 
    file_name = "detection_limit_high", 
    missing_id = c("limit_detection", "limit_quantification"),
    info = info
  )

  data
}


ctsm_check_subseries <- function(data, info) {

  # import_functions.R
  
  # if a time series has a subseries classification for some, but not all, 
  # records, delete the records that have no classification

  # avoids the situation where e.g. there is a subgroup for medium sized fish 
  # and the remaining samples are the small and large fish
  
  # note that this should probably be done other after 'grouping' variables
  # have been identified (e.g. sex for EROD) - future enhancement

  data <- tidyr::unite(
    data, 
    ".series", 
    any_of(c("station_code", "species", "determinand", "matrix")), 
    remove = FALSE
  )
  
  has_subgroup <- tapply(
    data$subseries, 
    data$.series, 
    function(x) any(!is.na(x))
  )
  
  has_subgroup <- names(has_subgroup)[has_subgroup]
  
  # remove extra data in these series that do not have a subgroup 
  
  not_ok <- data$.series %in% has_subgroup & is.na(data$subseries)
  
  message_txt <- paste0(
    "Some times series only have partial subseries classifications\n", 
    "            unclassified records are deleted"
  )
  
  data <- ctsm_check(
    data, 
    not_ok, 
    action = "delete", 
    message = message_txt,  
    file_name = "subseries_unclassified_data", 
    info = info
  )
  
  # replace remaining missing values with Not_applicable
  
  # data <- dplyr::mutate(data, subseries = replace_na(subseries, "Not_applicable"))

  data
}


ctsm_merge_auxiliary <- function(data, info) {

  # import_functions.R
  # merge auxiliary variables with data

  # identify auxiliary variables and split data set accordingly
    
  auxiliary_var <- ctsm_get_auxiliary(data$determinand, info)

  id <- data$determinand %in% auxiliary_var
    
  auxiliary <- data[id, ]
  data <- data[!id, ]
  
  
  # ensure all auxiliary variables are present in output, by creating a 
  # factor with levels given by auxiliary_var, and then splitting by this factor
  
  auxiliary$determinand <- factor(
    auxiliary$determinand, 
    levels = auxiliary_var
  )
  
  auxiliary <- split(auxiliary, auxiliary$determinand)
  
  
  # catch for LNMEA measured in WO and ES for birds - probably shouldn't happen because the 
  # sample will differ?  need to check
  
  
  for (aux_id in names(auxiliary)) {

    auxiliary_data <- auxiliary[[aux_id]]
    
    # standard merging variables, new variables and their names
    
    merge_id <- "sample"
    new_id <- "concentration"
    new_names <- aux_id
    
    # additional variables for some auxiliaries
    # need to make this more flexible - issue raised
    
    if (aux_id %in% c("DRYWT%", "LIPIDWT%", "CORG", "LOIGN", "AL", "LI")) {
      
      merge_id <- c(merge_id, "matrix")
      
      extra_id <- c(
        "censoring", "basis", "limit_detection", "limit_quantification", 
        "uncertainty"
      )
      
      if (aux_id %in% c("AL", "LI")) {
        extra_id <- c(extra_id, "digestion")
      }
      
      extra_names <- paste(aux_id, extra_id, sep = ".")
      
      new_id <- c(new_id, extra_id)
      new_names <- c(new_names, extra_names)
    } 
    
    auxiliary_data <- auxiliary_data[c(merge_id, new_id)]
    names(auxiliary_data) <- c(merge_id, new_names)
    
    
    # check for no replicate measurements when merging (merge should be many to one)
    
    dup_check <- do.call("paste", auxiliary_data[merge_id]) 
    
    if (any(duplicated(dup_check))) {

      message_txt <- paste0(
        "Cannot merge data with ", aux_id, " measurements due to\n",
        "            ambiguous matches; this might be because there are\n", 
        "            ", aux_id, " measurements for the same sample in more\n",
        "            than one matrix"
      )
      
      auxiliary_data <- ctsm_check(
        auxiliary_data, 
        dup_check %in% dup_check[duplicated(dup_check)], 
        action = "warning",
        message = message_txt,
        file_name = "replicated_auxiliary_measurements", 
        info = info
      )
      
      stop(
        "To proceed, edit the data to remove any ambiguous matches.\n", 
        "  Otherwise, contact the HARSAT development team to resolve the issue."
      )
      
    }
    
    
    # finally merge data
    
    data <- merge(data, auxiliary_data, all.x = TRUE)
  }
  
  data <- droplevels(data)
}

ctsm_convert_to_target_basis <- function(data, info, get_basis) {

  # location: import_functions.R
  # purpose:  convert data and auxiliary variables to their target basis as 
  #           determined by get_basis function

  # no conversion currently required for water, because all measurements
  # already on a wet weight basis, and no auxiliaries need converting
  
  # routines are more complicated than they need to be because of the need
  # to convert lipid data to a wet weight basis before further basis conversion
   
  # arguably can be done better by first merging dry and lipid data, converting 
  # lipid data, converting remaining data, and then merging with remaining 
  # auxiliaries - however, this requires more care with the exclude function and the 
  # merge_auxiliary routine, which should strictly match determinand with auxiliary
  
  # getting rid of extra 'basis' variables at the end is something that needs to 
  # be tidied up
  
  
  cat("   Converting data to appropriate basis for statistical analysis", fill = TRUE)
  
  data$new.basis <- get_basis(data, info)
  

  if (info$compartment == "biota") {
    
    # ensure lipidwt and drywt columns are present (if they haven't been supplied)
    
    is_lipid <- "LIPIDWT%" %in% names(data)
    is_dry <- "DRYWT%" %in% names(data)
     
    if (!is_lipid) {
      data[["LIPIDWT%"]] <- NA_real_
      data[["LIPIDWT%.censoring"]] <- NA_character_
    }
    
    if (!is_dry) {
      data[["DRYWT%"]] <- NA_real_
      data[["DRYWT%.censoring"]] <- NA_character_
    }
    
        
    # convert lipid weight data to wet weight basis (if present)
    
    if (is_lipid) {
    
      id <- c("", ".uncertainty", ".limit_detection", ".limit_quantification") 
      id <- paste0("LIPIDWT%", id)
      
      data[id] <- lapply(
        data[id], 
        ctsm_convert_basis,  
        from = data[["LIPIDWT%.basis"]], 
        to = "W",
        drywt = data[["DRYWT%"]], 
        drywt_censoring = data[["DRYWT%.censoring"]], 
        print_warning = FALSE
      )

    }
          
    
    # convert measurement data
    # print_warning gives the number of failures, which is the same for all id
    # so only print first time round
    
    id <- c("concentration", "uncertainty", "limit_detection", "limit_quantification")
    print_warning = c(TRUE, rep(FALSE, 3))
    
    data[id] <- mapply(
      FUN = ctsm_convert_basis,  
      conc = data[id], 
      print_warning = print_warning, 
      MoreArgs = list(
        from = data$basis, 
        to = data$new.basis,
        drywt = data[["DRYWT%"]], 
        drywt_censoring = data[["DRYWT%.censoring"]], 
        lipidwt = data[["LIPIDWT%"]], 
        lipidwt_censoring = data[["LIPIDWT%.censoring"]], 
        exclude = data$group %in% c("Imposex", "Metabolites", "Effects")
      ), 
      SIMPLIFY = FALSE
    )
    

    # remove dummy lipidwt and drywt variables
    
    if (!is_lipid) {
      data[["LIPIDWT%"]] <- NULL
      data[["LIPIDWT%.censoring"]] <- NULL
    }
    
    if (!is_dry) {
      data[["DRYWT%"]] <- NULL
      data[["DRYWT%.censoring"]] <- NULL
    }
    
  }
 
  
  if (info$compartment == "sediment") { 

    # ensure drywt columns are present (if they haven't been supplied)
    
    is_dry <- "DRYWT%" %in% names(data)
    
    if (!is_dry) {
      data[["DRYWT%"]] <- NA_real_
      data[["DRYWT%.censoring"]] <- NA_character_
    }
    
    
    # convert measurement data
    
    id <- c("concentration", "uncertainty", "limit_detection", "limit_quantification")
    print_warning = c(TRUE, rep(FALSE, 3))
    
    data[id] <- mapply(
      FUN = ctsm_convert_basis,  
      conc = data[id], 
      print_warning = print_warning, 
      MoreArgs = list(
        from = data$basis, 
        to = data$new.basis,
        drywt = data[["DRYWT%"]], 
        drywt_censoring = data[["DRYWT%.censoring"]]
      ), 
      SIMPLIFY = FALSE
    )
    
    # convert auxiliary (normalisation) data
    
    norm_id <- c("AL", "LI", "CORG", "LOIGN")
    ok <- norm_id %in% names(data)
    norm_id <- norm_id[ok]
    
    id <- c("", ".uncertainty", ".limit_detection", ".limit_quantification")
    
    conc_id <- paste0(rep(norm_id, each = length(id)), id)
    
    from_id <- paste0(norm_id, ".basis")
    from_id <- match(from_id, names(data))
    from_id <- rep(from_id, each = length(id))
    
    data[conc_id] <- mapply(
      FUN = ctsm_convert_basis, 
      conc = data[conc_id], 
      from = data[from_id], 
      MoreArgs = list(
        to = data$new.basis,
        drywt = data[["DRYWT%"]], 
        drywt_censoring = data[["DRYWT%.censoring"]], 
        print_warning = FALSE
      ), 
      SIMPLIFY = FALSE
    )


    # remove dummy drywt variables
    
    if (!is_dry) {
      data[["DRYWT%"]] <- NULL
      data[["DRYWT%.censoring"]] <- NULL
    }

  }


  # drop unwanted basis variables (apart from C13D and N15D which haven't been converted)
  # retain new.basis for time series structure
  
  id <- grep(".basis", names(data), value = TRUE)
  id <- setdiff(id, c("new.basis", "C13D.basis", "N15D.basis"))
  id <- setdiff(names(data), id)
  data <- droplevels(data[id])
  
  data
}


#' @export
ctsm_normalise_sediment <- function(data, station_dictionary, info, control) {
  
  # normalises sediment concentrations
  
  # method supplied by control
  
  ctsm_normalise_default <- list(
    metals = list(method = "pivot", normaliser = "AL", extra = NULL), 
    organics = list(method = "simple", normaliser = "CORG", value = 2.5), 
    exclude = NULL
  )
  
  control <- modifyList(ctsm_normalise_default, control)
  
  
  # save non-normalised concentrations and censoring codes for plotting purposes 
  # later on
  # also save uncertainties just in case
  
  data <- dplyr::mutate(
    data, 
    concOriginal = .data$concentration,     
    censoringOriginal = .data$censoring,
    uncrtOriginal = .data$uncertainty
  )
  
  
  # exclude any data that do not need to be normalised 
  # can do this globally with method = "none", but useful e.g. in the OSPAR 
  #   assessment where sediments in the Iberian Sea and Gulf of Cadiz are not
  #   normalised
  
  # excluded stations are evaluated in the station dictionary
  # corresponding rows in the data are then identified
  
  
  if (!is.null(control$exclude)) {
    exclude_id <- eval(control$exclude, station_dictionary, parent.frame())
    exclude_id <- station_dictionary[exclude_id, "station_code"]
    exclude_id <- data$station_code %in% exclude_id
  } else {
    exclude_id <- FALSE
  }
  
  if (any(exclude_id)) {
    excluded_data <- data[exclude_id, ]
    data <- data[!exclude_id, ]
  }
  
  
  # split into metals and organics and then normalise each with e.g. AL and 
  # CORG respectively
  
  groupID <- factor(
    data$group == "Metals", 
    levels = c(TRUE, FALSE), 
    labels = c("metals", "organics")
  )
  
  data <- split(data, groupID)
  
  
  data <- mapply(
    names(data), data, control[names(data)],
    SIMPLIFY = FALSE, 
    FUN = function(group, data, control) {
      
      # check normalisation method fully specified by control
      
      if (! control$method %in% c("none", "simple", "pivot")) {
        stop("uncoded normalisation method specified: current methods are none, simple, pivot")
      }
      
      
      # exit if nothing to be done
      
      if (control$method == "none") {
        message("   No normalisation for ", group)
        return(data)
      }
      
      
      # extract normaliser and print summary information
      
      normaliser <- control$normaliser
      
      if (! normaliser %in% names(data)) {
        stop("Normaliser ", normaliser, " not found in data")
      }
      
      switch(
        control$method,
        simple = {
          unit <- ctsm_get_info(
            info$determinand, normaliser, "unit", "sediment", sep = "_"
          )
          message("   Normalising ", group, " to ", control$value, unit, " ", normaliser)
        },
        pivot = message("   Normalising ", group, " to ", normaliser, " using pivot values")
      )
      
      
      # function to get normaliser variables
      
      getNdata <- function(x) {
        id <- paste(normaliser, x, sep = ".")
        data[[id]]
      }      
      
      
      # get concentration and normaliser 
      
      Cm <- data$concentration
      Nm <- data[[normaliser]]
      
      var_Cm <- data$uncertainty ^ 2
      var_Nm <- getNdata("uncertainty") ^ 2
      
      
      # get pivot values
      
      if (control$method == "simple") {
        Cx <- 0
        Nx <- 0
        Nss <- control$value
      }    
      
      if (control$method == "pivot") {
        
        # get pivot data and make row.names the appropriate combination of 
        #determinand and digestion
        
        pivot <- info$pivot_values
        pivot <- pivot[pivot$determinand %in% c(as.character(data$determinand), normaliser), ]
        rownames(pivot) <- with(pivot, paste(determinand, digestion))
        pivot <- droplevels(pivot)
        
        
        # get digestion for both contaminant and normaliser (nn indicates not required)
        # if missing digestion for normaliser, assume the strongest method for the 
        #   same sample matrix combination
        # might both be missing if QA does not contain the method_extraction code 
        
        Cdigestion <- data$digestion
        Ndigestion <- getNdata("digestion")
        
        if ("PNL" %in% Cdigestion) {
          warning(
            "ad-hoc fix to deal with new NL method - must resolve for next assessment", 
            call. = FALSE
          )
        }
        
        
        notOK <- !is.na(data[[normaliser]]) & is.na(Ndigestion)
        
        if (any(notOK)) {
          message(
            "   Inferring missing normaliser digestion from corresponding ",
            "contaminant digestion"
          )
          
          # order factor levels in increasing strength of digestion
          
          Cdigestion <- factor(Cdigestion, levels = c("nn", "PNL", "Pw", "Ps", "Tot"))
          Ndigestion <- factor(Ndigestion, levels = c("nn", "PNL", "Pw", "Ps", "Tot"))
          
          
          # get strongest digestion for same sample matrix combination
          
          linkID <- with(data, paste(sample, matrix))
          
          replacement <- tapply(Cdigestion, linkID, function(x) {
            if (all(is.na(x))) return(NA)
            id <- max(as.numeric(x), na.rm = TRUE)
            levels(x)[id]
          })
          
          replacement <- factor(replacement[linkID], levels(Cdigestion))
          
          Ndigestion[notOK] <- replacement[notOK]
        }
        
        
        # check that all pivot values are there (provided digestion information present)
        
        CpivotID <- paste(data$determinand, Cdigestion)
        ok <- CpivotID %in% rownames(pivot) | is.na(Cdigestion)
        if (!all(ok)) {
          stop(
            'Not found in pivot information file: ', 
            paste(CpivotID[notOK], collapse = ", ")
          )
        }
        
        NpivotID <- paste(normaliser, Ndigestion)
        NpivotID[Ndigestion %in% "PNL" & data$determinand %in% "CR"] <- "AL Ps"
        ok <- NpivotID %in% rownames(pivot) | is.na(Ndigestion)
        if (!all(ok)) {
          stop(
            'Not found in pivot information file: ', 
            paste(NpivotID[notOK], collapse = ", ")
          )
        }        
        
        
        # get pivot values and normalise
        # note natural variability of pivot values varCx and varNx forced to be zero
        
        Cx <- pivot[CpivotID, "Nx"]          # pivot values
        Nx <- pivot[NpivotID, "Nx"]     
        
        Nss <- pivot[NpivotID, "Nss"]        # pivot conc
        
        
        # special treatment for CR Nss values - see Foppe communication
        
        if ("CR" %in% data$determinand) {
          
          if (! normaliser %in% c("AL", "LI")) {
            stop("Nss values for CR not coded for ", normaliser)
          }
          
          message(
            "   Warning: Nss values for CR hard-wired to 5.0 (AL) or ", 
            "52 (LI) for all digestions"
          )
          
          Nss[data$determinand %in% "CR"] <- switch(normaliser, AL = 5.0, LI = 52, NA)
        }
        
        
        # manual adjustment for French data from Region II
        # Table 7 of https://archimer.ifremer.fr/doc/00461/57294/59374.pdf
        
        if ("country" %in% names(station_dictionary) && 
            "France" %in% station_dictionary$country) {
          
          station_id <- station_dictionary$country %in% "France" & 
            station_dictionary$OSPAR_region %in% "2"
          station_id <- station_dictionary[station_id, "station_code"]
          id <- data$station_code %in% station_id
          
          if (any(id)) {
            message(
              "   Warning: French data from Region II normalised with ", 
              "user-supplied pivot values"
            )
            
            if (!normaliser %in% "AL") {
              stop(
                "French metal data in Region II with no pivot values as ",
                "normaliser is not AL"
              )
            }
            
            if (!all(Ndigestion[id] %in% "Tot")) {
              stop(
                "French metal data in Region II with partial digestion so no ",
                "pivot values"
              )
            }
            
            if (any(data$determinand[id] %in% "AS")) {
              stop("French arsenic data in Region II with no pivot values")
            } 
            
            Cx[id & data$determinand %in% "CD"] <-  0.05
            Cx[id & data$determinand %in% "HG"] <-  0.002
            Cx[id & data$determinand %in% "PB"] <-  8.5
            Cx[id & data$determinand %in% "CR"] <-  6.6
            Cx[id & data$determinand %in% "CU"] <-  0.29
            Cx[id & data$determinand %in% "NI"] <-  1.8
            Cx[id & data$determinand %in% "ZN"] <- 10.0
            
            Nx[id]  <- 1.02
            Nss[id] <- 5.0
          }
          
        }
        
      }      
      
      
      # normalise concentrations and put back into standard variables
      
      out <- ctsm_normalise_calculate(Cm, Nm, Nss, var_Cm, var_Nm, Cx, Nx)
      
      concentration <- out$Css 
      uncertainty <- sqrt(out$var_Css)
      
      
      # check no normalisers are less-thans; 
      # if so, normalised concentration should be a greater than, but
      # haven't coded this yet
      
      notOK <- getNdata("censoring") %in% c("<", "D", "Q")
      if (any(notOK)) {
        message('   Removing sediment data where normaliser is a less than')
        concentration[notOK] <- NA
      }
      
      data$concentration <- concentration
      data$uncertainty <- uncertainty
      
      data
    })
  
  data <- unsplit(data, groupID)
  
  if (any(exclude_id)) {
    data <- dplyr::bind_rows(data, excluded_data)
  }
  
  data
}

#' @export
ctsm_normalise_sediment_HELCOM <- function(data, station_dictionary, info, control) {
  
  # normalises sediment concentrations
  
  # method supplied by control
  
  ctsm_normalise_default <- list(
    metals = list(method = "pivot", normaliser = "AL", extra = NULL), 
    copper = list(method = "hybrid", normaliser = "CORG", value = 5),
    organics = list(method = "simple", normaliser = "CORG", value = 5), 
    exclude = NULL
  )
  
  control <- modifyList(ctsm_normalise_default, control)
  
  
  # save non-normalised concentrations and censoring codes for plotting purposes 
  # later on
  # also save uncertainties just in case
  
  data <- dplyr::mutate(
    data, 
    concOriginal = .data$concentration,     
    censoringOriginal = .data$censoring,
    uncrtOriginal = .data$uncertainty
  )
  
  
  # exclude any data that do not need to be normalised 
  # can do this globally with method = "none", but useful e.g. in the OSPAR 
  #   assessment where sediments in the Iberian Sea and Gulf of Cadiz are not
  #   normalised
  
  # excluded stations are evaluated in the station dictionary
  # corresponding rows in the data are then identified
  
  
  if (!is.null(control$exclude)) {
    exclude_id <- eval(control$exclude, station_dictionary, parent.frame())
    exclude_id <- station_dictionary[exclude_id, "station_code"]
    exclude_id <- data$station_code %in% exclude_id
  } else {
    exclude_id <- FALSE
  }
  
  if (any(exclude_id)) {
    excluded_data <- data[exclude_id, ]
    data <- data[!exclude_id, ]
  }
  
  
  # make ad-hoc change to deal with LOIGN
  # must undo at the end of the code

  data <- dplyr::mutate(
    data, 
    .tmp = CORG,
    .tmp.censoring = CORG.censoring,
    .tmp.uncertainty = CORG.uncertainty,
    CORG = dplyr::if_else(is.na(.tmp), 0.35 * LOIGN, CORG),
    CORG.censoring = dplyr::if_else(
      is.na(.tmp), 
      as.character(LOIGN.censoring), 
      as.character(CORG.censoring)
    ),
    CORG.censoring = factor(CORG.censoring),
    CORG.uncertainty = dplyr::if_else(is.na(.tmp), 0.35 * LOIGN.uncertainty, CORG.uncertainty)
  )

  
  # split into metals (CD, PB), copper and organics and then normalise each 
  # with AL, CORG (LOIGN) and CORG (LOIGN) respectively#
  
  groupID <- dplyr::case_when(
    data$determinand %in% c("CD", "PB") ~ "metals",
    data$determinand %in% "CU"          ~ "copper",
    TRUE                                ~ "organics"
  ) 

  groupID <- factor(groupID)
    
  data <- split(data, groupID)
  
  
  data <- mapply(
    names(data), data, control[names(data)],
    SIMPLIFY = FALSE, 
    FUN = function(group, data, control) {
      
      # check normalisation method fully specified by control
      
      if (! control$method %in% c("none", "simple", "pivot", "hybrid")) {
        stop(
          "uncoded normalisation method specified: ", 
          "current methods are none, simple, pivot, hybrid"
        )
      }
      
      
      # exit if nothing to be done
      
      if (control$method == "none") {
        message("   No normalisation for ", group)
        return(data)
      }
      
      
      # extract normaliser and print summary information
      
      normaliser <- control$normaliser
      
      if (! normaliser %in% names(data)) {
        stop("Normaliser ", normaliser, " not found in data")
      }
      
      switch(
        control$method,
        simple = {
          unit <- ctsm_get_info(
            info$determinand, normaliser, "unit", "sediment", sep = "_"
          )
          message("   Normalising ", group, " to ", control$value, unit, " ", normaliser)
        },
        pivot = message("   Normalising ", group, " to ", normaliser, " using pivot values"),
        hybrid = message("   Normalising ", group, " to ", normaliser, " using pivot values")
      )
      
      
      # function to get normaliser variables
      
      getNdata <- function(x) {
        id <- paste(normaliser, x, sep = ".")
        data[[id]]
      }      
      
      
      # get concentration and normaliser 
      
      Cm <- data$concentration
      Nm <- data[[normaliser]]
      
      var_Cm <- data$uncertainty ^ 2
      var_Nm <- getNdata("uncertainty") ^ 2
      
      
      # get pivot values
      
      if (control$method == "simple") {
        Cx <- 0
        Nx <- 0
        Nss <- control$value
      }    
      
      if (control$method == "pivot") {
        
        # get pivot data and make row.names the appropriate combination of 
        #determinand and digestion
        
        pivot <- info$pivot_values
        pivot <- pivot[pivot$determinand %in% c(as.character(data$determinand), normaliser), ]
        rownames(pivot) <- with(pivot, paste(determinand, digestion))
        pivot <- droplevels(pivot)
        
        
        # get digestion for both contaminant and normaliser (nn indicates not required)
        # if missing digestion for normaliser, assume the strongest method for the 
        #   same sample matrix combination
        # might both be missing if QA does not contain the method_extraction code 
        
        Cdigestion <- data$digestion
        Ndigestion <- getNdata("digestion")
        

        notOK <- !is.na(data[[normaliser]]) & is.na(Ndigestion)
        
        if (any(notOK)) {
          message(
            "   Inferring missing normaliser digestion from corresponding ",
            "contaminant digestion"
          )
          
          # order factor levels in increasing strength of digestion
          
          Cdigestion <- factor(Cdigestion, levels = c("nn", "PNL", "Pw", "Ps", "Tot"))
          Ndigestion <- factor(Ndigestion, levels = c("nn", "PNL", "Pw", "Ps", "Tot"))
          
          
          # get strongest digestion for same sample matrix combination
          
          linkID <- with(data, paste(sample, matrix))
          
          replacement <- tapply(Cdigestion, linkID, function(x) {
            if (all(is.na(x))) return(NA)
            id <- max(as.numeric(x), na.rm = TRUE)
            levels(x)[id]
          })
          
          replacement <- factor(replacement[linkID], levels(Cdigestion))
          
          Ndigestion[notOK] <- replacement[notOK]
        }
        
        
        # check that all pivot values are there (provided digestion information present)
        
        CpivotID <- paste(data$determinand, Cdigestion)
        ok <- CpivotID %in% rownames(pivot) | is.na(Cdigestion)
        if (!all(ok)) {
          stop(
            'Not found in pivot information file: ', 
            paste(CpivotID[notOK], collapse = ", ")
          )
        }
        
        NpivotID <- paste(normaliser, Ndigestion)
        ok <- NpivotID %in% rownames(pivot) | is.na(Ndigestion)
        if (!all(ok)) {
          stop(
            'Not found in pivot information file: ', 
            paste(NpivotID[notOK], collapse = ", ")
          )
        }        
        
        
        # get pivot values and normalise
        # note natural variability of pivot values varCx and varNx forced to be zero
        
        Cx <- pivot[CpivotID, "Nx"]          # pivot values
        Nx <- pivot[NpivotID, "Nx"]     
        
        Nss <- pivot[NpivotID, "Nss"]        # pivot conc
        
      }      
      

      if (control$method == "hybrid") {
        
        # get pivot data and make row.names the appropriate combination of 
        # determinand and digestion
        
        pivot <- info$pivot_values
        pivot <- pivot[pivot$determinand %in% c(as.character(data$determinand), normaliser), ]
        rownames(pivot) <- with(pivot, paste(determinand, digestion))
        pivot <- droplevels(pivot)
        
        
        # get digestion for contaminant 
        # not needed for normaliser (CORG)

        Cdigestion <- data$digestion


        # check that all pivot values are there (provided digestion information present)
        
        CpivotID <- paste(data$determinand, Cdigestion)
        ok <- CpivotID %in% rownames(pivot) | is.na(Cdigestion)
        if (!all(ok)) {
          stop(
            'Not found in pivot information file: ', 
            paste(CpivotID[notOK], collapse = ", ")
          )
        }
        

        # get pivot values and normalise
        # note natural variability of pivot values varCx and varNx forced to be zero
        
        Cx <- pivot[CpivotID, "Nx"]          # pivot values
        Nx <- 0     
        
        Nss <- control$value                 # pivot conc
        
      }      
      
            
      # normalise concentrations and put back into standard variables
      
      out <- ctsm_normalise_calculate(Cm, Nm, Nss, var_Cm, var_Nm, Cx, Nx)
      
      concentration <- out$Css 
      uncertainty <- sqrt(out$var_Css)
      
      
      # check no normalisers are less-thans; 
      # if so, normalised concentration should be a greater than, but
      # haven't coded this yet
      
      notOK <- getNdata("censoring") %in% c("<", "D", "Q")
      if (any(notOK)) {
        message('   Removing sediment data where normaliser is a less than')
        concentration[notOK] <- NA
      }
      
      data$concentration <- concentration
      data$uncertainty <- uncertainty
      
      data
    })
  
  data <- unsplit(data, groupID)


  data <- dplyr::mutate(
    data, 
    CORG = .tmp,
    CORG.censoring = .tmp.censoring,
    CORG.uncertainty = .tmp.uncertainty,
    .tmp = NULL,
    .tmp.censoring = NULL,
    .tmp.uncertainty = NULL
  )

  if (any(exclude_id)) {
    data <- dplyr::bind_rows(data, excluded_data)
  }
  
  data
}

#' @export
ctsm_normalise_biota_HELCOM <- function(data, station_dictionary, info, control) {
  
  # normalises fish concentrations in contaminants other than metals or 
  # organofluorines to 5% lipid
  # these concentrations are expressed on a wet weight basis
  
  # method supplied by control
  
  ctsm_normalise_default <- list(
    lipid = list(method = "simple", value = 5), 
    other = list(method = "none"), 
    exclude = NULL
  )
  
  control <- modifyList(ctsm_normalise_default, control)
  
  
  # save non-normalised concentrations and censoring codes for plotting purposes 
  # later on
  # also save uncertainties just in case
  
  data <- dplyr::mutate(
    data, 
    concOriginal = .data$concentration,     
    censoringOriginal = .data$censoring,
    uncrtOriginal = .data$uncertainty
  )
  
  
  # exclude any data that do not need to be normalised 
  # can do this globally with method = "none", but useful e.g. in the OSPAR 
  #   assessment where sediments in the Iberian Sea and Gulf of Cadiz are not
  #   normalised
  
  # excluded stations are evaluated in the station dictionary
  # corresponding rows in the data are then identified
  
  
  if (!is.null(control$exclude)) {
    exclude_id <- eval(control$exclude, station_dictionary, parent.frame())
    exclude_id <- station_dictionary[exclude_id, "station_code"]
    exclude_id <- data$station_code %in% exclude_id
  } else {
    exclude_id <- FALSE
  }
  
  if (any(exclude_id)) {
    excluded_data <- data[exclude_id, ]
    data <- data[!exclude_id, ]
  }
  
  
  # split into contaminants that are to be normalised to 5% lipid and 
  # others
  
  groupID <- dplyr::if_else(
    data$species_group %in% "Fish" & 
      !(data$group  %in% c("Metals", "Organofluorines", "Metabolites")), 
    "lipid", 
    "other"
  )
  
  groupID <- factor(groupID)
  
  data <- split(data, groupID)
  
  
  data <- mapply(
    names(data), data, control[names(data)],
    SIMPLIFY = FALSE, 
    FUN = function(group, data, control) {
      
      # check normalisation method fully specified by control
      
      if (! control$method %in% c("none", "simple")) {
        stop("uncoded normalisation method specified: current methods are none, simple")
      }
      
      
      # exit if nothing to be done
      
      if (control$method == "none") {
        message("   No normalisation for ", group)
        return(data)
      }
      
      
      # normalise to a specified value of lipid content
      # data are already on a lipid basis (i.e. 100% lipid)
      
      message("   Normalising ", group, " to ", control$value, "%")

      conversion <- dplyr::if_else(
        data[["LIPIDWT%.censoring"]] %in% "", 
        control$value / data[["LIPIDWT%"]],
        NA_real_
      )
      
      data$concentration <- data$concentration * conversion
      data$uncertainty <- data$uncertainty * conversion
      
      data
    })
  
  data <- unsplit(data, groupID)
  

  if (any(exclude_id)) {
    data <- dplyr::bind_rows(data, excluded_data)
  }
  
  data
}




ctsm_normalise_calculate <- function(Cm, Nm, Nss, var_Cm, var_Nm, Cx, Nx, var_Cx = 0, var_Nx = 0) {
  
  # calculates normalised concentrations and their uncertainties based on pivot relationship
  # setting Nx = 0, Cx = 0 equates to a directly proportional relationship (e.g. organics with CORG)
  
  # Cm = measured concentration of contaminant
  # Nm = measured concentration of normaliser
  # Nss = target concentration of normaliser(e.g. 5% aluminium, 2.5% organic carbon)
  # var_Cm = variance of Cm
  # var_Nm = variance of Nm
  # Cx, Nx = pivot values of contaminant and normaliser
  # var_Cx, var_Nx = variances of Cx, Nx - usually set to zero in a time series analysis
                                              
                                        
  # allows for measured concentrations that are below the pivot concentration 
  # (provided the normalised concentration is still positive)
  
  ok <- Nm > Nx  &  Cm * (Nss - Nx) > Cx * (Nss - Nm)     
  
  Css <- dplyr::if_else(ok, (Cm - Cx) * (Nss - Nx) / (Nm - Nx) + Cx, NA_real_)
  
  var_Css <- dplyr::if_else(ok, ((Cm - Cx) / (Nm - Nx))^2, NA_real_)
  var_Css <- (Nss - Nx)^2 * (var_Cm + var_Nm * var_Css) + (Nss - Nm)^2 * (var_Cx + var_Nx * var_Css)
  var_Css <- dplyr::if_else(ok, var_Css / ((Nm - Nx)^2), NA_real_)
  
  data.frame(Cm, Nm, Css, var_Css)
}



ctsm_estimate_uncertainty <- function(data, response_id, info) {

  # import_functions.R
  
  # estimating missing uncertainties
  
  # only returns uncertainty so can modify data object at will
  
  # check response_id is valid (could be "concentration" or an auxiliary variable)
  # check compartment is valid and not a variable in data
  
  stopifnot(
    is.character(response_id),
    length(response_id) == 1L,
    response_id %in% names(data),
    
    is.character(info$compartment),
    length(info$compartment) == 1L
  )
  

  # standardise variable names when response_id is an auxiliary variable

  if (response_id != "concentration") {
    data <- dplyr::mutate(data, determinand = rep(response_id, nrow(data)))

    var1 <- c("uncertainty", "censoring", "limit_detection", "limit_quantification")  
    var2  <- paste(response_id, var1, sep = ".")
    
    data[c("concentration", var1)] <- data[c(response_id, var2)]
  }
    

  # is there anything to do?

  missing_id <- is.na(data$uncertainty)
  
  if (!any(missing_id)) { 
    return(data$uncertainty)
  }

  
  # get two components of variance
  
  data$sd_constant = ctsm_get_info(
    info$determinand, 
    data$determinand, 
    "sd_constant", 
    info$compartment, 
    na_action = "output_ok", 
    sep = "_"
  )
  
  data$sd_variable = ctsm_get_info(
    info$determinand, 
    data$determinand, 
    "sd_variable", 
    info$compartment, 
    na_action = "output_ok", 
    sep = "_"
  )
  
  
  # adjust sd_constant to correct basis (
  # for biota sd_constant is on a wet weight basis but the sample might be on 
  # a dry or lipid weight

  if (info$compartment == "biota") {
    data$sd_constant <- ctsm_convert_basis(
      data$sd_constant, 
      "W", 
      data$new.basis, 
      data[["DRYWT%"]], 
      data[["LIPIDWT%"]], 
      drywt_censoring = data[["DRYWT%.censoring"]], 
      lipidwt_censoring = data[["LIPIDWT%.censoring"]], 
      exclude = data$group %in% c("Imposex", "Effects", "Metabolites"),
      print_warning = FALSE
    )
  }

  
  # adjust sd_constant to use limit and censoring information where possible 
  # however censoring < is difficult to work with as it is inconsistent with both lod and loq so
  #   take lower of concentration / 3 and sd_constant
  # censoring = "" is also difficult to work with because we can't trust the lod and loq 
  #   option - again take lower of many options
  
  data$sd_constant <- dplyr::case_when(
    data$censoring %in% "D" ~ data$limit_detection / 3,
    data$censoring %in% "Q" ~ data$limit_quantification / 10,
    data$censoring %in% "<" ~ pmin(data$sd_constant, 
                               data$concentration / 3, 
                               na.rm = TRUE),
    data$censoring %in% ""  ~ pmin(data$sd_constant, 
                               data$limit_detection / 3, 
                               data$limit_quantification / 10, 
                               data$concentration / 3, 
                               na.rm = TRUE),
    TRUE                ~ NA_real_
  )
  
  
  data <- dplyr::mutate(
    data, 
    estimated_uncrt = .data$sd_constant ^ 2 + 
      (.data$sd_variable * .data$concentration) ^ 2,
    estimated_uncrt = sqrt(.data$estimated_uncrt),
    uncertainty = dplyr::if_else(
      is.na(.data$uncertainty), .data$estimated_uncrt, .data$uncertainty
    )
  )

  data$uncertainty
}



#' Convert tin concentrations
#'
#' Convert tin concentrations to cation concentrations.
#' Also `change_unit` moves units from tin units to conventional units
#' @export 
ctsm_TBT_convert <- function(
  data, subset, action, from = c("tin", "cation"), 
  convert_var = c("value", "limit_detection", "limit_quantification", "uncertainty")) {
  
  # relabel turns tin labels to cation labels
  # convert moves tin concentrations to cation concentrations
  # change_unit moves unit from tin unit to conventional unit

  from <- match.arg(from)
  
  stopifnot(action %in% c("relabel", "convert", "change_unit"))
  
  if (from == "cation" & "relabel" %in% action) 
    stop("shouldn't relabel cation determinands")
  
  if ("convert" %in% action & "uncertainty" %in% convert_var & 
      ! "unit_uncertainty" %in% names(data))
    stop("unit_uncertainty not provided to allow conversion of uncertainty")
  
  
  id <- data$determinand %in% switch(
    from, 
    tin = c("TBTIN", "TPTIN", "DBTIN", "DPTIN", "MBTIN", "MPTIN"), 
    cation = c("TBSN+", "TPSN+", "DBSN+", "DPSN+", "MBSN+", "MPSN+")
  )
  
  if (!missing(subset)) 
    id <- id & eval(substitute(subset), data, parent.frame()) 
  
  if (sum(id) == 0)
    return(data)
  
  if ("relabel" %in% action) {
    data[id, ] <- within(data[id, ], {
      determinand <- gsub("TIN", "SN+", determinand, fixed = TRUE)
    })
  }
  
  if ("convert" %in% action) {
    if ("uncertainty" %in% convert_var) 
      with(data[id, ], stopifnot(is.na(unit_uncertainty) | unit_uncertainty %in% "U2"))
    
    tin_id <- with(data[id, ], substring(determinand, 1, 2))
    conversion <- dplyr::recode(
      tin_id, DB = "1.96", MB = "1.48", TB = "2.44", TP = "2.95", MP = "1.65", DP = "2.30")
    conversion <- as.numeric(conversion)
    
    data[id, convert_var] <- data[id, convert_var] * conversion
  }
  
  if ("change_unit" %in% action) {
    data[id, ] <- within(data[id, ], {
      stopifnot(unit %in% c("ug sn/kg", "gm/g"))
      unit <- dplyr::recode(unit, "ug sn/kg" = "ug/kg", "gm/g" = "g/g")
    })
  }  
  
  data
}  

#' @export
ctsm_TBT_remove_early_data <- function(ctsm_obj, recent_year = min(ctsm_obj$info$recent_years)) {
  
  tin <- c("TBTIN", "TPTIN", "DBTIN", "DPTIN", "MBTIN", "MPTIN")
  cation <- c("TBSN+", "TPSN+", "DBSN+", "DPSN+", "MBSN+", "MPSN+")
  
  id <- with(ctsm_obj$data, determinand %in% c(tin, cation))
  id <- factor(id, levels = c("TRUE", "FALSE"))
  
  data <- split(ctsm_obj$data, id)
  
  tin_data <- data[["TRUE"]]
  
  # drop data with no stations
  
  ok <- with(tin_data, !is.na(station_code))
  tin_data <- tin_data[ok, ]
  
  # drop time series with no data in recent years - avoids sorting out old data that won't get used
  
  wk <- with(tin_data, tapply(year, station_code, max))
  ok <- with(tin_data, station_code %in% names(wk)[wk >= recent_year])
  tin_data <- tin_data[ok, ]
  
  ctsm_obj$data <- rbind(data[["FALSE"]], tin_data)
  
  ctsm_obj
} 


   
