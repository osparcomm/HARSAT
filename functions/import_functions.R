# main import and data processing functions ----


# read in data  ---- 

ctsm_read_data <- function(
  compartment = c("sediment", "biota", "water"), 
  purpose = c("OSPAR", "HELCOM", "AMAP"), 
  contaminants, 
  stations, 
  QA, 
  path = ".", 
  extraction, 
  max_year, 
  control = list(), 
  data_format = c("old", "new")) {

  # import_functions.R
  # reads in data from an ICES extraction
  
  
  data_format <- match.arg(data_format)
  
  
  # validate arguments
  
  compartment <- match.arg(compartment)
  purpose <- match.arg(purpose)
  
  if (!is.character(path) | length(path) > 1) {
    stop("path should be a single character string")
  }
  
  # max_year is the maximum data year allowed 
  # print warning if there are more recent data which have been submitted 'early'
  
  if (length(max_year) != 1 |
      !isTRUE(all.equal(max_year, as.integer(max_year)))) {
    stop("max_year must be a single (integer valued) year")
  }


  # sort out extraction date
    
  if (length(extraction) > 1) {
    stop("extraction must be a single date")
  }
  extraction <- as.POSIXlt(extraction)
  
  extraction_year <- as.numeric(format(extraction, "%Y"))
  extraction_month <- as.numeric(format(extraction, "%m"))
  
  if (extraction_year < max_year) {
    stop("extraction predates some monitoring data")
  }
  
  extraction_text <- format(extraction, "%d %B %Y")
  if (substring(extraction_text, 1, 1) == "0") { 
    extraction_text <- substring(extraction_text, 2)
  }
  
  
  # set up control (info) structure 
  # some of this should probably go in ctsm_create_timeSeries
    
  info <- list(
    compartment = compartment, 
    purpose = purpose, 
    extraction = extraction_text, 
    max_year = max_year 
  )
  
  control_default <- ctsm_control_default(purpose)
  
  control <- ctsm_control_modify(control_default, control)
  
  if (any(names(control) %in% names(info))) {
    warning("possible conflict between function arguments and control")
  }
  
  info <- append(info, control)
  

  # construct recent_years in which there must be some monitoring data
  
  info$recent_years <- seq(max_year - control$reporting_window + 1, max_year)
  

    
  # read in station dictionary, contaminant and biological effects data and QA data
  
  stations <- ctsm_read_stations(stations, path, info$purpose, info$region_id)
  
  data <- ctsm_read_contaminants(
    contaminants, path, info$purpose, info$region_id, data_format
  )
  
  if (data_format == "old") 
    QA <- ctsm_read_QA(QA, path, purpose)
  

  # check no data after max_year
  
  if (any(data$year > max_year)) { 
    warning("data submitted after max_year", call. = FALSE)
  }
  
  if (!any(data$year == max_year)) { 
    warning("no data in max_year - possible error in function call", call. = FALSE)
  }
  
    
  out <- list(
    call = match.call(), 
    info = info,
    data = data, 
    stations = stations
  )
  
  if (data_format == "old") out$QA <- QA
  
  out
}


ctsm_control_default <- function(purpose) {
  
  # import functions
  # sets up default values that control the assessment
  
  # reporting_window is set to 6 to match the MSFD reporting cycle
  # series with no data in the most recent_window are excluded
  
  # region_id are the names or the relevant regions in the ICES extraction
  # region_names are the names that will be used for reporting 
  
  # all_in_region is a logical that determines whether all data (and stations)
  # must be in a region

  region_id <- switch(
    purpose, 
    OSPAR = c("OSPAR_region", "OSPAR_subregion"),
    HELCOM = c("HELCOM_subbasin", "HELCOM_L3", "HELCOM_L4"),
    NULL
  )
  
  region_names <- region_id
  
  all_in_region <- !is.null(region_id)
    
  list(
    reporting_window = 6L, 
    region_id = region_id,
    region_names = region_names,
    all_in_region = all_in_region
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

  control
}


ctsm_read_stations <- function(infile, path, purpose, region_id) {

  # import functions
  # read in station dictionary
  
  infile <- file.path(path, infile)
  cat("Reading station dictionary from '", infile, "'\n", sep = "")

  stations <- read.table(
    infile, strip.white = TRUE, sep = "\t",  header = TRUE, quote = "\"", 
    na.strings = c("", "NULL"), fileEncoding = "UTF-8-BOM", comment.char = ""
  )
  
  
  # rename columns to suit!
  
  names(stations) <- gsub("Station_", "", names(stations), fixed = TRUE)

  stations <- dplyr::rename(
    stations, 
    code = Code,
    country_ISO = Country,
    country = Country_CNTRY,
    station = Name,
    name = LongName,
    latitude = Latitude,
    longitude = Longitude,
    startYear = ActiveFromDate,
    endYear = ActiveUntilDate,
    parent_code = AsmtMimeParent,
    replacedBy = ReplacedBy,
    programGovernance = ProgramGovernance,
    dataType = DataType
  )

  if (purpose %in% c("OSPAR", "AMAP")) {
    stations <- dplyr::rename(stations, offshore = OSPAR_shore)
  }  

  # turn following variables into a character
  # code and parent code could be kept as integers, but safer to leave as 
  # characters until have decided how we are going to merge with non-ICES data

  id <- c(region_id, "code", "parent_code")
  
  stations <- dplyr::mutate(
    stations, 
    dplyr::across(any_of(id), as.character)
  )
  
  stations
}


ctsm_read_contaminants <- function(infile, path, purpose, region_id, data_format) {
  
  # import functions
  # read in contaminant (and biological effects) data
  
  infile <- file.path(path, infile)
  cat("\nReading contaminant and biological effects data from '", infile, "'\n", 
      sep = "")

  if (data_format == "old") {  

    data <- read.table(
      infile, strip.white = TRUE, sep = "\t", header = TRUE, quote = "\"" , 
      na.strings = c("", "NULL"), fileEncoding = "UTF-8-BOM", comment.char = ""
    )
    
  } else {

    # specify type of each column
    # station codes and tblsampleid are converted to characters
    #   to make it easier to combine with external data - this could be revisited
    # could maybe filter out unwanted variables at this stage to save space 
    #   e.g. don't need Ospar variables for a Helcom assessment
    
    data <- read.table(
      infile, 
      strip.white = TRUE, 
      sep = "\t", 
      header = TRUE, 
      quote = "\"" , 
      na.strings = c("", "NULL"), 
      fileEncoding = "UTF-8-BOM", 
      comment.char = "",
      colClasses = c(
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
    )
    
  }
    
    
  # check regional identifiers are in the extraction 

  if (data_format == "new" & purpose == "HELCOM") {

    data <- dplyr::rename(
      data,
      HELCOM_subbasin = helcom_subbasin,
      HELCOM_L3 = helcom_l3,
      HELCOM_L4 = helcom_l4,
    )

  }
    
  pos <- names(data) %in% region_id
  if (sum(pos) != length(region_id)) {
    stop("not all regional identifiers are in the data extraction")
  }
  
  names(data)[!pos] <- tolower(names(data)[!pos])     # just in case!

  
  # create more useful names
  # for biota, tblsampleid is the species, tblbioid gives the subsample (individual) 
  # which with matrix gives the unique sample id
  # for sediment, tblsampleid and matrix give the unique sample id
  # to streamline, could make sample = tblbioid (biota) or tblsampleid (sediment)

  
  if (data_format == "new") {
    
    data <- dplyr::rename(
      data,
      submitted.station = statn, 
      sd_name = sd_code_name,
      sd_code = sd_code_match,
      station_name = sd_name_final,
      station_code = sd_code_final
    )
    
  } else if (purpose %in% c("OSPAR", "AMAP")) {
    
    data <- dplyr::rename(
      data,
      offshore = ospar_shore,
      submitted.station = stationname, 
      sd_name = sd_stationname,
      sd_code = sd_stationcode,
      station_name = sd_asmt_stationname,
      station_code = sd_asmt_stationcode,
    )
    
  } else if (purpose %in% "HELCOM") {
    
    data <- dplyr::rename(
      data,
      submitted.station = statn, 
      sd_name = sd_stationname,
      sd_code = sd_stationcode,
      station_name = sd_replacedby_stationname,
      station_code = sd_replacedby_stationcode,
    )
    
  } 
  
  
  # variables common across purpose and compartments

  data <- dplyr::rename(
    data,
    year = myear, 
    determinand = param, 
    matrix = matrx, 
    unit = munit, 
    qalink = tblanalysisid,
    uncertainty = uncrt,
    methodUncertainty = metcu, 
    replicate = tblparamid, 
    sample = tblsampleid, 
    limit_detection = detli,
    limit_quantification = lmqnt,
    upload = tbluploadid
  )
  

  # compartment specific variables
    
  var_id <- c("tblbioid", "sexco", "dephu")
  replacement <- c("sub.sample", "sex", "depth")
  
  pos <- match(var_id, names(data), nomatch = 0)
  
  if (any(pos > 0)) {
    ok <- pos > 0
    pos <- pos[ok]
    replacement <- replacement[ok]
    names(data)[pos] <- replacement
  }
  
  
  # ensure further consistency
  
  data <- dplyr::mutate(
    data, 
    determinand = toupper(.data$determinand),
    unit = tolower(.data$unit)
  )
  
  
  if (data_format == "old") {

    id <- c(
      region_id, "qflag", "sample", "sub.sample", "sd_code", "station_code", 
      "sd_name", "station_name"
    )
    
    data <- dplyr::mutate(
      data, 
      dplyr::across(any_of(id), as.character)
    )
  
  }
      
  data
}


ctsm_read_QA <- function(QA, path, purpose) {
  
  # import functions
  # read in method data (and additional crm information)

  # read in data
  
  infile <- file.path(path, QA)
  cat("\nReading QA data from '", infile, "'\n", sep = "")
  
  crm <- read.table(
    infile, strip.white = TRUE, header = TRUE, quote = "\"" , sep = "\t",
    na.strings = c("", "NULL"), fileEncoding = "UTF-8-BOM", comment.char = ""
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
    qalink = tblanalysisid
  )

  crm <- dplyr::mutate(crm, determinand = toupper(.data$determinand))
  
  crm
}



ctsm_create_timeSeries <- function(
  ctsm.obj, 
  determinands = ctsm_get_determinands(ctsm.obj$info$compartment), 
  determinands.control = NULL, 
  oddity.path = "oddities", 
  return_early = FALSE, 
  print_code_warnings = FALSE, 
  output = c("time_series", "uncertainties"), 
  normalise = FALSE, 
  normalise.control = list()) {

  # import_functions.R
  # cleans data and turns into time series structures ready for assessment
  
  library(tidyverse)


  # arguments
  
  # determinands: character string of determinands to be assessed; default is to 
  #   take values from determinand reference table
  
  determinands <- force(determinands)
  
  # output: 
  # time_series is default
  # uncertainties provides data for estimating missing uncertainties - converts all 
  #   concentrations to a wet weight basis

  output = match.arg(output)
  
  # normalisation can either be a logical (TRUE uses default normalisation function)
  # or a function
  
  if (length(normalise) != 1L) {
    stop("normalise should be a length 1 logical or a function")
  }
  
  if (!(is.logical(normalise) | is.function(normalise))) {
    stop("normalise should be a length 1 logical or a function")
  }
  

  # identify data format (old or new) - temporary requirement during development
  
  if (is.null(ctsm.obj$call$data_format)) {
    data_format <- "old"
  } else {
    data_format <- ctsm.obj$call$data_format
  }
  
  
  # get key data structures, and initialise output

  out <- list(call = match.call(), call.data = ctsm.obj$call, info = ctsm.obj$info)
  
  info <- ctsm.obj$info
  station.dictionary <- ctsm.obj$stations
  data <- ctsm.obj$data
  
  if (data_format == "old") QA <- ctsm.obj$QA
  
  rm(ctsm.obj)
 

  is.recent <- function(year) year %in% info$recent_years


  # lots of data cleansing - first ensure oddity directory exists and back up
  # any previous oddity files

  oddity.dir <- ctsm.initialise.oddities(oddity.path)


  # clean station dictionary, and provide rownames that will link to stationID 
  # in the data file

  station.dictionary <- ctsm.clean.stations(
    info$purpose, station.dictionary, info$compartment)


  # clean data - get rid of data with no station, construct stationID 
  # and sampleID, and retain variables that are going to be used 

  data <- ctsm.clean.contaminants(info$purpose, data, info$compartment)


  # retains determinands of interest, including auxiliary determinands and those
  # required by determinands.control$variables 
  # checks all determinands of interest are recognised by info.determinand

  wk <- ctsm_check_determinands(
    info$compartment, 
    data, 
    determinands, 
    determinands.control
  )
  
  data <- wk$data
  determinands <- wk$determinands
  determinands.control <- wk$control

  
  # clean QA data - creates unique qaID to link data and QA files 

  if (data_format == "old") {
  
    wk <- ctsm_clean_QA(info$purpose, QA, data, info$compartment)
    
    data <- wk$data
    QA <- wk$QA

  } else {
    
    data$qaID = rep(0, nrow(data))
    
  }


  cat("\nFurther cleaning of data\n")
  
  id <- c("station", "date", "filtered", "species", "sampleID", "matrix", "determinand")
  id <- intersect(id, names(data))
  ord <- do.call("order", data[id])
  data <- data[ord, ]
  

  # check all stations are in station dictionary - abort if not - means there
  # has been an error in extraction

  ok <- data$station %in% rownames(station.dictionary)
  data <- ctsm.check(
    data, !ok, action = "delete", message = "Stations not in station dictionary", 
    fileName = "unidentified stations", merge.stations = FALSE)
  if (!all(ok)) 
    stop("Error in extraction: stations not in station dictionary", call. = FALSE)


  # drop stations (sediment) or station / species combinations (biota) with no
  # data in the relatively recent period 
  # recent years reduces size of data and removes legacy species not in info files
  
  cat(
    "   Dropping stations with no data between", min(info$recent_years), "and", 
    max(info$recent_years), "\n"
  )

  # recognised_species was replaced by species in line 507 -> species reference table
  
  if ("species" %in% names(data)) {
    data <- within(data, {
      species <- as.character(species)
      species <- ifelse(
        species %in% row.names(info.species), 
        as.character(info.species[species, "reference_species"]), 
        species
      )
      species <- factor(species)
    })
  }
  
  id.names <- intersect(c("station", "species"), names(data))
  id <- do.call("paste", data[id.names])
  ok <- id %in% id[is.recent(data$year)]
  data <- droplevels(data[ok, ])
  

  # drop data corresponding to stations outside the region (mainly OSPAR or HELCOM requirement)
  
  id <- is.na(station.dictionary[as.character(data$station), "region"])
  if (any(id)) {
    cat("   Dropping stations with no associated region in station dictionary\n")
    data <- droplevels(subset(data, !id))
  }


  # add variables that are going to be useful throughout
  # pargroup in ICES extraction, but can also be got from info.determinand
  
  data$group <- ctsm_get_info(
    "determinand", data$determinand, "group", info$compartment, sep = "_"
  )

  if (info$compartment == "biota") {
    data$family <- ctsm_get_info("species", data$species, "species_group")
  }
  
  if (!"pargroup" %in% names(data)) {
    data$pargroup <- ctsm_get_info("determinand", data$determinand, "pargroup")
  }
  

  # drop samples which only have auxiliary data
  
  ok <- with(data, sampleID %in% sampleID[group != "Auxiliary"])
  if (any(!ok)) {
    cat("   Dropping samples with only auxiliary variables\n")
    data <- data[ok, ]
  }  

  # remove species that are not to be assessed and check family and sex appropriate 

  for (varID in c("species", "family", "sex")) data <- ctsm.check0(data, varID, info$compartment)


  # check all determinands have a valid matrix, basis and unit
  
  for (varID in c("basis", "matrix", "unit")) data <- ctsm.check0(data, varID, info$compartment)


  # get method of analysis and method of chemical extraction
  
  if (data_format == "old") {
    id <- c("metoa", "metcx") 
    data[id] <- QA[as.character(data$qaID), id]
  }

  
  # check all data have a valid method of analysis, value and number in pool

  for (varID in c("metoa", "value", "noinp")) data <- ctsm.check0(data, varID, info$compartment)


  # get digestion method for metals in sediments and check that these are valid
  
  if (info$compartment == "sediment") {
    data$digestion <- ctsm_get_digestion(data)
  }


  # check no replicate measurements (some are genuine replicates, mostly from
  # early years) but will just delete these for simplicity (taking the average
  # has been tried, but gets really messy when there are less thans, different
  # detection limits (sediment), uncertainties etc.)
  
  # not sure of the best place for this
  # should come before the merging with auxiliary variables, otherwise get multiple matches
  # should also come before the linking of e.g. CHR with CHRTR, otherwise the duplicates aren't 
  # necessarily the correct ones
  
  data <- ctsm.check(
    data, paste(sampleID, determinand, matrix), action = "delete.dups", 
    message = "Replicate measurements, only first retained", 
    fileName = "replicate measurements")


  # ensure qflag, limit of detection and limit of quantification are consistent
  
  # define function for testing equality - also need a relative component to cope with some tiny
  # value submitted with units g/g
  
  my_near <- function(x, y) {
    abs <- near(x, y) 
    rel <- case_when(
      x > 0 ~ abs((x - y) / y) < 1e-5,
      TRUE ~ TRUE
    )
    abs & rel
  }
    

  # check detection limits associated with data are positive 

  if (print_code_warnings)
    warning('Need to make checking of detection limits determinand specific', call. = FALSE)
  
  data <- ctsm.check(
    data, limit_detection <= 0, action = "make.NA", 
    message = "Non-positive detection limits", 
    fileName = "non positive det limits", missingID = "limit_detection")
  
  data <- ctsm.check(
    data, limit_quantification <= 0, action = "make.NA", 
    message = "Non-positive quantification limits", 
    fileName = "non positive quant limits", missingID = "limit_quantification")
  
  
  # limit_quantification must be greater than limit_detection - otherwise set both to missing

  data <- ctsm.check(
    data,  limit_quantification <= limit_detection, action = "make.NA", 
    message = "Limit of quantification less than limit of detection", 
    fileName = "limits inconsistent", 
    missingID = c("limit_detection", "limit_quantification"))
  
  
  # check for valid values of qflag
  # qflag cannot contain a > (although possible for some biological effects - need to revisit)
  # also check for other unrecognised characters

  if (print_code_warnings)
    warning(
      "Need to make qflag determinand specific when biological effects are introduced", 
      call. = FALSE)

  data <- within(data, {
    levels(qflag) <- c(levels(qflag), "")
    qflag[is.na(qflag)] <- ""
    qflag <- recode(qflag, "<~D" = "D", "D~<" = "D", "<~Q" = "Q", "Q~<" = "Q")
  })
  
  data <- ctsm.check(
    data, ! qflag %in% c("", "D", "Q", "<"), action = "delete", 
    message = "Unrecognised qflag values", fileName = "qflags unrecognised")
  
  
  # if qflag = D, then value must equal detection_limit
  # if qflag = Q, then value must equal limit_quantification
  
  id <- with(data, {
    out1 <- qflag == "D" & (is.na(limit_detection) | !my_near(value, limit_detection))
    out2 <- qflag == "Q" & (is.na(limit_quantification) | !my_near(value, limit_quantification))
    out1 | out2
  })
  
  ctsm.check(
    data, id, action = "warning",
    message = "Qflag D and Q inconsistent with respective limits",
    fileName = "qflags and limits inconsistent")
  
  
  # resolve these inconsistencies
  
  data <- within(data, {
    qflag <- case_when(
      qflag %in% "" ~ "",
      !is.na(limit_detection) & my_near(value, limit_detection) ~ "D",
      !is.na(limit_quantification) & my_near(value, limit_quantification) ~ "Q",
      TRUE ~ "<"
    )
    qflag <- factor(qflag, levels = c("", "D", "Q", "<"))  
  })    

  
  # check limit_detection is less than (or equal to) concentration
  # NB would need to be revised if limits are submitted in different units to concentrations
  
  data <- ctsm.check(
    data, 
    id = qflag %in% c("", "<") & limit_detection > value, 
    action = "make.NA", 
    message = "Detection limit higher than data", 
    fileName = "detection limit high", 
    missingID = c("limit_detection", "limit_quantification")
  )
  

  # convert uncertainty into standard deviations, and remove any associated variables
  
  data <- ctsm.check(
    data, !is.na(uncertainty) & uncertainty <= 0, action = "make.NA", 
    message = "Non-positive uncertainties", 
    fileName = "non positive uncertainties", missingID = "uncertainty")
  
  data <- mutate(
    data, 
    uncertainty_sd = case_when(
      methodUncertainty %in% "U2" ~ uncertainty / 2, 
      methodUncertainty %in% "%" ~ value * uncertainty / 100, 
      TRUE ~ uncertainty
    ),
    uncertainty_rel = 100 * (uncertainty_sd / value)
  )                 

  wk_id <- match("methodUncertainty", names(data))
  wk_n <- ncol(data)
  data <- data[c(
    names(data)[1:wk_id], 
    "uncertainty_sd", "uncertainty_rel", 
    names(data)[(wk_id+1):(wk_n-2)])]
  
  ctsm.check(
    data, !is.na(uncertainty) & uncertainty_rel >= 100, action = "warning", 
    message = "Large uncertainties", fileName = "large uncertainties")
  
  # keep these uncertainties if output = uncertainties, so we can see everything
  
  if (output == "time_series") {
    data <- mutate(
      data, 
      uncertainty_sd = if_else(.data$uncertainty_rel < 100, .data$uncertainty_sd, NA_real_)
    )
  }
    
  data <- within(data, {
    uncertainty <- uncertainty_sd
    rm(methodUncertainty, uncertainty_sd, uncertainty_rel)
  })
  

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
      list(data = data, keep = i, drop = wk$det)
    )
  }  

  # drop any remaining unwanted determinands (from sum and perhaps bespoke functions);
  # could make this more elegant!
  
  id <- c(determinands, ctsm_get_auxiliary(determinands, info$compartment))
  
  data <- filter(data, determinand %in% id)

  # set rownames to NULL(ie. back to auto numeric)
  
  rownames(data) <- NULL

  cat("\nCreating time series data\n")  

  data <- data[setdiff(names(data), c("qalink", "alabo"))]

  
  # create new.unit and concentration columns comprising the details from the
  # determinand file in the information folder, required to get correct unit details
  
  data <- within(data, {
    new.unit <- ctsm_get_info(
      "determinand", determinand, "unit", info$compartment, sep = "_"
    )
    concentration <- value
  })
  

  # convert data to conventional units - see "information determinands.csv"
  # NB value is in original units, concentration is in converted units
  # this applies to all determinands (including auxiliary variables)
  # can't change bases here because need to merge with auxiliary variables

  id <- c("concentration", "uncertainty", "limit_detection", "limit_quantification")

  data[id] <- lapply(data[id], convert_units, from = data$unit, to = data$new.unit) 


  # merge auxiliary data with determinand data
  # weights and sediment normalisers are merged by sampleID and matrix
  # others just by sampleID (for now)

  if (print_code_warnings)
    warning(
      "merging of auxiliary variables in ctsm.create.timeseries; mergeID hard-wired",
      call. = FALSE)
  
  auxData <- droplevels(data[data$group == "Auxiliary", ])
  

  # ad-hoc fix to ensure 'key' auxiliary variables are present
  
  auxData <- within(auxData, {
    levels(determinand) <- switch(
      info$compartment,
      biota = unique(c(levels(determinand), "DRYWT%", "LIPIDWT%", "LNMEA")), 
      sediment = unique(c(levels(determinand), "DRYWT%")),
      levels(determinand)
    )
    
    if (info$purpose %in% "AMAP")
      levels(determinand) <- unique(c(levels(determinand), "AGMEA", "C13D", "N15D"))
  })
  
  auxData <- split(auxData, auxData$determinand)
  
  data <- data[data$group != "Auxiliary", ]

  
  # catch for LNMEA measured in WO and ES for birds - probably shouldn't happen because the 
  # sampleID will differ?  need to check
  
  if ("LNMEA" %in% names(auxData)) {
    if (print_code_warnings)
      warning("need to resolve merging of LNMEA with contaminant data", call. = FALSE)
    
    if (anyDuplicated(auxData[["LNMEA"]][["sampleID"]]))
      stop("sampleID with LNMEA measurements in more than one matrix", call. = FALSE)
  }
  
  # and similarly check that C13D and N15D are only measured once in each sub.sample
  
  if (any(c("C13D", "N15D") %in% names(auxData))) {
    if (print_code_warnings)
      warning("need to resolve merging of C13D and N15D with contaminant data", call. = FALSE)
    
    if (anyDuplicated(auxData[["C13D"]][["sampleID"]]))
      stop("sampleID with C13D measurements in more than one matrix", call. = FALSE)
  
    if (anyDuplicated(auxData[["N15D"]][["sampleID"]]))
      stop("sampleID with N15D measurements in more than one matrix", call. = FALSE)
  }
  
  
  for (i in names(auxData)) {

    if (i %in% c("DRYWT%", "LIPIDWT%", "CORG", "LOIGN")) {
      mergeID <- c("sampleID", "matrix")
      newID <- c(
        "concentration", "qflag", "basis", "limit_detection", "limit_quantification", 
        "uncertainty", "qaID"
      )
      newNames <- c(mergeID, i, paste(i, newID[-1], sep = "."))
    } else if (i %in% c("AL", "LI")) {
      mergeID <- c("sampleID", "matrix")
      newID <- c(
        "concentration", "qflag", "basis", "limit_detection", "limit_quantification", 
        "uncertainty", "qaID", "digestion"
      )
      newNames <- c(mergeID, i, paste(i, newID[-1], sep = "."))
    } else if (i %in% c("C13D", "N15D")) {
      mergeID <- "sampleID"
      newID <- c("concentration", "matrix", "basis")
      newNames <- c(mergeID, i, paste(i, newID[-1], sep = "."))
    } else {
      mergeID <- "sampleID"
      newID <- "concentration"
      newNames <- c(mergeID, i)
    }
    
    wk <- auxData[[i]][c(mergeID, newID)]
    names(wk) <- newNames
    data <- merge(data, wk, all.x = TRUE)
  }
  
  data <- droplevels(data)


  # impute %femalepop when missing and sex = 1 - write out remaining
  # missing values for correction
  
  if (info$compartment == "biota")
    data <- ctsm.imposex.check.femalepop(data) 
  

  # convert data to appropriate basis

  # biota - first convert LIPIDWT% to wet weight (where necessary)
  
  if (info$compartment == "biota") {

    id <- paste0("LIPIDWT%", c("", ".uncertainty", ".limit_detection", ".limit_quantification"))
    
    data[id] <- lapply(
      data[id], 
      convert.basis,  
      from = data[["LIPIDWT%.basis"]], 
      to = rep("W", nrow(data)),
      drywt = data[["DRYWT%"]], 
      drywt.qflag = data[["DRYWT%.qflag"]], 
      lipidwt = NA, 
      lipidwt.qflag = NA
    )
  }


  # output = uncertainties
  # convert data to a wet weight basis (biota and water) or dry weight basis (sediment)
  #   for estimating constant and proportional errors and return
  # assume water data are all W - no DRYWT data for water

  if (output == "uncertainties") {
  
    data$new.basis <- switch(
      info$compartment, 
      sediment = factor(rep("D", nrow(data))),
      factor(rep("W", nrow(data)))
    )
    
    id <- c("concentration", "uncertainty", "limit_detection", "limit_quantification")
    
    if (info$compartment %in% c("biota", "sediment")) {
    
      data[id] <- lapply(
        data[id], 
        convert.basis,  
        from = data$basis, 
        to = data$new.basis,
        drywt = data[["DRYWT%"]], 
        drywt.qflag = data[["DRYWT%.qflag"]], 
        lipidwt = switch(info$compartment, biota = data[["LIPIDWT%"]], NA), 
        lipidwt.qflag = switch(info$compartment, biota = data[["LIPIDWT%.qflag"]], NA), 
        exclude = data$group %in% c("Imposex", "Metabolites", "Effects")
      )
      
      # THIS CAN BE STREAMLINED!!!!!
      
      if (info$compartment == "sediment") {
        
        wk_suffix <- c(
          "", ".uncertainty", ".limit_detection", ".limit_quantification"
        )
        
        id <- paste0("AL", wk_suffix)
        data[id] <- lapply(
          data[id], 
          convert.basis, 
          from = data$AL.basis,
          to = data$new.basis,
          drywt = data[["DRYWT%"]], 
          drywt.qflag = data[["DRYWT%.qflag"]]
        )
        
        id <- paste0("LI", wk_suffix)
        data[id] <- lapply(
          data[id], 
          convert.basis, 
          from = data$AL.basis,
          to = data$new.basis,
          drywt = data[["DRYWT%"]], 
          drywt.qflag = data[["DRYWT%.qflag"]]
        )
        
        id <- paste0("CORG", wk_suffix)
        data[id] <- lapply(
          data[id], 
          convert.basis, 
          from = data$CORG.basis,
          to = data$new.basis,
          drywt = data[["DRYWT%"]], 
          drywt.qflag = data[["DRYWT%.qflag"]]
        )
        
        if ("LOIGN" %in% names(data)) {
          
          id <- paste0("LOIGN", wk_suffix)
          data[id] <- lapply(
            data[id], 
            convert.basis, 
            from = data$LOIGN.basis,
            to = data$new.basis,
            drywt = data[["DRYWT%"]], 
            drywt.qflag = data[["DRYWT%.qflag"]]
          )
        }        
      }
    
    }  
      
    # drop unwanted basis variables (apart from C13D and N15D which haven't been converted)
    # need to retain new.basis (AMAP) when the basis is chosen depending on what has been 
    # submitted
    
    id <- grep(".basis", names(data), value = TRUE)
    id <- setdiff(id, c("new.basis", "C13D.basis", "N15D.basis"))
    id <- setdiff(names(data), id)
    data <- droplevels(data[id])
    
    out  = c(
      out, 
      ctsm_import_value(data, station.dictionary, info, print_code_warnings) 
    )
    
    if (data_format == "old") out$QA <- QA
    
    return(out)

  }
    

  # output = time_series
  # convert data to basis of assessment
  

  if (info$purpose %in% c("HELCOM", "OSPAR")) {
    
    cat("   Converting data to appropriate basis for statistical analysis", fill = TRUE)

    data <- mutate(
      data, 
      new.basis = get_basis(
        info$purpose, info$compartment, .data$group, .data$matrix, .data$determinand,
        species = switch(info$compartment, biota = .data$species, NA)
      ),
      new.basis = factor(new.basis)
    )
  }
  

  if (info$purpose %in% "AMAP") {
    
    cat(
      "   Converting data to common basis (within station, species, matrix and determinand group)", 
      fill = TRUE)
    
    data <- get_basis(info$purpose, data, info$compartment)
  }
    

  # don't convert water because already assumed to be on a wet weight basis
  # no DRYWT% data and convert.basis falls over
  
  if (info$compartment %in% c("biota", "sediment")) {
  
    id <- c("concentration", "uncertainty", "limit_detection", "limit_quantification")
    data[id] <- lapply(
      data[id], 
      convert.basis,  
      from = data$basis, 
      to = data$new.basis,
      drywt = data[["DRYWT%"]], 
      drywt.qflag = data[["DRYWT%.qflag"]], 
      lipidwt = switch(info$compartment, biota = data[["LIPIDWT%"]], NA), 
      lipidwt.qflag = switch(info$compartment, biota = data[["LIPIDWT%.qflag"]], NA), 
      exclude = data$group %in% c("Imposex", "Metabolites", "Effects")
    )
    
  }

    
  if (info$compartment == "biota" & any(c("C13D", "N15D") %in% names(data)))
    message("   Warning: Isotope ratios have not been converted to target basis")


  if (info$compartment == "sediment") {

    wk_suffix <- c(
      "", ".uncertainty", ".limit_detection", ".limit_quantification"
    )
    
    id <- paste0("AL", wk_suffix)
    data[id] <- lapply(
      data[id], 
      convert.basis, 
      from = data$AL.basis,
      to = data$new.basis,
      drywt = data[["DRYWT%"]], 
      drywt.qflag = data[["DRYWT%.qflag"]]
    )
    
    id <- paste0("LI", wk_suffix)
    data[id] <- lapply(
      data[id], 
      convert.basis, 
      from = data$AL.basis,
      to = data$new.basis,
      drywt = data[["DRYWT%"]], 
      drywt.qflag = data[["DRYWT%.qflag"]]
    )
    
    id <- paste0("CORG", wk_suffix)
    data[id] <- lapply(
      data[id], 
      convert.basis, 
      from = data$CORG.basis,
      to = data$new.basis,
      drywt = data[["DRYWT%"]], 
      drywt.qflag = data[["DRYWT%.qflag"]]
    )
    
    if ("LOIGN" %in% names(data)) {
      
      id <- paste0("LOIGN", wk_suffix)
      data[id] <- lapply(
        data[id], 
        convert.basis, 
        from = data$LOIGN.basis,
        to = data$new.basis,
        drywt = data[["DRYWT%"]], 
        drywt.qflag = data[["DRYWT%.qflag"]]
      )
    }        
  }

  
  # drop unwanted basis variables (apart from C13D and N15D which haven't been converted)
  # retain new.basis for time series structure

  id <- grep(".basis", names(data), value = TRUE)
  id <- setdiff(id, c("new.basis", "C13D.basis", "N15D.basis"))
  id <- setdiff(names(data), id)
  data <- droplevels(data[id])

  if (return_early) {
    out  = c(
      out, 
      ctsm.import.value(
        data, 
        station.dictionary, 
        info$compartment, 
        info$purpose, 
        print_code_warnings)
    )
    
    if (data_format == "old") out$QA <- QA
    
    return(out)
  }


  # estimate missing uncertainties

  data$uncertainty <- ctsm.estimate.uncertainty(data, "concentration", info$compartment)

  if (info$compartment == "sediment") {
    data$AL.uncertainty <- ctsm.estimate.uncertainty(data, "AL", info$compartment)
    data$LI.uncertainty <- ctsm.estimate.uncertainty(data, "LI", info$compartment)
    data$CORG.uncertainty <- ctsm.estimate.uncertainty(data, "CORG", info$compartment)
    
    if ("LOIGN" %in% names(data)) {
      data$LOIGN.uncertainty <- ctsm.estimate.uncertainty(data, "LOIGN", info$compartment)
    }    
  }

  if (info$compartment == "biota") {

    # restrict shellfish data to 'sensible' months - must make this region specific

    nok <- with(data, {
      month <- substring(months(as.Date(date)), 1, 3)
      contaminants <- ! group %in% c("Effects", "Imposex", "Metabolites")
      contaminants & 
        family %in% c("Bivalvia", "Gastropoda") & 
        month %in% switch(
          info$purpose, 
          AMAP = c("Apr", "May", "Jun", "Jul"),
          OSPAR = c("Apr", "May", "Jun", "Jul"),
          HELCOM = c("Apr", "May", "Jun", "Jul"))
    })

    if (any(nok)) {
      txt <- switch(
        info$purpose, 
        AMAP = "April through July",
        OSPAR = "April through July", 
        HELCOM = "April through July")
      cat("   Dropping shellfish data from", txt, "(spawning season)\n")
      data <- data[!nok, ]
    }
  }


  # normalise data
  # e.g. sediment data to 5% aluminium and 2.5% organic carbon (OSPAR)
  # biota data to 5% lipid (HELCOM)

  if (is.logical(normalise) && normalise) {
    data <- switch(
      info$compartment,
      biota = 
        ctsm_normalise_biota(data, QA, station.dictionary, normalise.control),
      sediment = 
        ctsm_normalise_sediment(data, station.dictionary, normalise.control),
      water = stop("there is no default normalisation function for water")
    )
  } else if (is.function(normalise)) {
    data <- normalise(data, station.dictionary, normalise.control)
  }
    

  # remove concentrations where cv of uncertainty > 100%
  # tidy up missing data for uncertainties and qflag
  
  data <- within(data, {
    concentration[uncertainty > concentration] <- NA
    uncertainty[is.na(concentration)] <- NA
    qflag[is.na(concentration)] <- NA
  })
    

  # drop groups of data at stations with no data in recent years

  cat("   Dropping groups of compounds / stations with no data between", 
      min(info$recent_years), "and", max(info$recent_years), "\n")
  id_names <- intersect(c("station", "filtered", "species", "group"), names(data))
  id <- do.call("paste", data[id_names])
  ok <- id %in% id[is.recent(data$year) & !is.na(data$concentration)]
  data <- data[ok, ]
  

  # check no replicates (typically across matrices)
  # frustratingly have made ctsm.check too clever to generalise this

  if (info$compartment == "biota") {

    data <- ctsm.check(
      data, paste(sampleID, determinand, matrix), action = "delete.dups", 
      message = "Measurements still replicated, only first retained", 
      fileName = "replicate measurements extra")

  } else {
    
    data <- ctsm.check(
      data, paste(sampleID, determinand), action = "delete.dups", 
      message = "Measurements replicated across matrices, only first retained", 
      fileName = "replicate measurements extra")
    
  }
    
  out  = c(
    out, 
    ctsm_import_value(data, station.dictionary, info, print_code_warnings)
  )
  
  if (data_format == "old") out$QA <- QA
 
  out
}


# default routine to check (and clean) data attributes when creating timeSeries

ctsm.check <- function(
  data, id, action = c("delete", "warning", "make.NA", "delete.dups"), 
  message, fileName, merge.stations = TRUE, missingID) {

  # action options: 
  # delete - delete oddities (in id) 
  # warning - print out oddities, but don't change data 
  # make.NA - replace values of missingID corresponding to oddities with NA 
  # delete.dups - delete any duplicated values identified in id, but print out 
  #   all instances of duplication to see what is going on
  
  action <- match.arg(action)

  # evaluate logical identifying odditie:
  # if action = remove duplicates, then find all the duplicated values in id
  # else evaluate id to obtain logical of oddities

  if (action == "delete.dups") {
    duplicateVar <- eval(substitute(id), data, parent.frame())
    odditiesID <- duplicated(duplicateVar) | duplicated(duplicateVar, fromLast = TRUE)
  } 
  else {
    odditiesID <- eval(substitute(id), data, parent.frame())
    if (!is.logical(odditiesID)) stop("'id' must be logical")
    odditiesID <- odditiesID & !is.na(odditiesID)
  }

  
  # do nothing if no oddities
  
  if (!any(odditiesID)) 
    if (action == "warning") return() else return(data)
 

  # get oddities
  
  oddities <- data[odditiesID, ]


  # merge with station dictionary if required - gets dictionary from top environment

  if (merge.stations) {
    station.dictionary <- evalq(station.dictionary, sys.frame(1))
    oddities.stations <- station.dictionary[as.character(oddities$station), ]
    names(oddities.stations) <- paste("SD", names(oddities.stations), sep = ".")
    oddities <- oddities[setdiff(names(oddities), "station")]        
    oddities <- cbind(oddities.stations, oddities)
  }

  
  # sort out columns that are all na and that are logicals for more robust writing
  
  wk <- vapply(oddities, function (x) all(is.na(x)), NA)
  oddities[wk] <- lapply(oddities[wk], function(x) rep_len("", length(x)))

  wk <- vapply(oddities, is.logical, NA)
  oddities[wk] <- lapply(oddities[wk], as.character)


  # write out oddities
  
  oddity.dir <- evalq(oddity.dir, sys.frame(1))
  oddity.file <- paste0(fileName, ".csv")
  message("   Warning: ", message, 
      switch(action, delete = ": deleted data in '", ": see '"), oddity.file)
  
  write.csv(oddities, file.path(oddity.dir, oddity.file), row.names = FALSE, na = "")
  
  #write.xlsx(oddities, oddity.file, sheetName, append = TRUE, showNA = FALSE, 
  #           row.names = FALSE)

  if (action == "warning") return()
  
  data <- switch(
    action, 
    delete = data[!odditiesID, ],
    make.NA = {data[odditiesID, missingID] <- NA ; data},
    delete.dups = data[!duplicated(duplicateVar), ]
  )
  
  return(droplevels(data))        
}
  

ctsm.check2 <- function(data, action, message, fileName, merge.stations = TRUE) {

  # check that all records have a valid action - if not, need to go back to
  # check function
  
  if (any(is.na(action))) {
    oddity.dir <- evalq(oddity.dir, sys.frame(1))
    oddity.file <- paste0(fileName, ".csv")
    
    write.csv(data[is.na(action), ], file.path(oddity.dir, oddity.file), 
              row.names = FALSE, na = "")
    stop("Not all cases considered when checking ", message, ": see '", oddity.file, 
         "'\n", sep = "")
  }
  
  
  # do nothing if no oddities
  
  if (all(action %in% "none")) return(data)
  
  
  # if no errors or warnings, just delete data and return
  
  if (all(action %in% c("none", "delete"))) {
    data <- data[action %in% "none", ]  
    return(droplevels(data))
  }  
  

  # get oddities
  
  notOK <- action %in% c("error", "warning")
  oddities <- data[notOK, ]
  
  
  # merge with station dictionary if required - gets dictionary from top environment
  
  if (merge.stations) {
    station.dictionary <- evalq(station.dictionary, sys.frame(1))
    oddities.stations <- station.dictionary[as.character(oddities$station), ]
    names(oddities.stations) <- paste("SD", names(oddities.stations), sep = ".")
    oddities <- oddities[setdiff(names(oddities), "station")]        
    oddities <- cbind(oddities.stations, oddities)
  }
  
  
  # make action the first column
  
  oddities <- cbind(action = action[notOK], oddities)
  
  
  # sort out columns that are all na and that are logicals for more robust writing
  
  wk <- vapply(oddities, function (x) all(is.na(x)), NA)
  oddities[wk] <- lapply(oddities[wk], function(x) rep_len("", length(x)))
  
  wk <- vapply(oddities, is.logical, NA)
  oddities[wk] <- lapply(oddities[wk], as.character)
  
  
  # write out oddities
  
  oddity.dir <- evalq(oddity.dir, sys.frame(1))
  oddity.file <- paste0(fileName, ".csv")
  message("   Warning: ", message, ": see '", oddity.file)
  
  write.csv(oddities, file.path(oddity.dir, oddity.file), row.names = FALSE, na = "")
  
  
  # delete rows where appropriate
  
  data <- data[action %in% c("none", "warning"), ]
  
  return(droplevels(data))
}


ctsm_import_value <- function(data, station_dictionary, info, print_code_warnings) {
  
  # import_functions.R
  # identifies timeseries in the data and creates timeSeries metadata object
  
  # order data and select variables of interest

  order.names = c(
    "station", "species", "filtered", "year", "sex", "sample", "sub.sample", "group", "determinand")
  order.names <- intersect(order.names, names(data))
  data <- data[do.call(order, data[order.names]), ]
  
  if (print_code_warnings)
    warning("auxiliary variables hard-wired in ctsm.import.value: need to resolve", call. = FALSE)
  
  auxiliary <- c(
    "AL", "LI", "CORG", "LOIGN", "LNMEA", "AGMEA", "DRYWT%", "LIPIDWT%", 
    "%FEMALEPOP", "CMT-QC-NR", "MNC-QC-NR", "C13D", "N15D"
  )

  out.names = c(
    "station", "latitude", "longitude", "filtered", "species", "sex", "depth",
    "year", "date", "time", "sample", "sub.sample", "sampleID", 
    "matrix", "AMAP_group", "group", "determinand", "basis", "unit", "value", "metoa", "noinp", 
    "concOriginal", "qflagOriginal", "uncrtOriginal", 
    "concentration", "new.basis", "new.unit", "qflag",  
    "qaID", "limit_detection", "limit_quantification", "uncertainty",  
    paste(
      rep(auxiliary, each = 6), 
      c("", ".qflag", ".qaID", ".limit_detection", ".limit_quantification", ".uncertainty"), 
      sep = ""
    ),
    "C13D.basis", "C13D.matrix", "N15D.basis", "N15D.matrix")

  out.names <- intersect(out.names, names(data))
  data <- data[out.names]
  
  row.names(data) <- NULL
  data <- droplevels(data)
  
  
  # get seriesID and timeSeries structure 
  
  id <- c("station", "determinand")
  if (info$compartment == "biota")
    id <- switch(
      info$purpose, 
      AMAP = c(id, "species", "matrix", "AMAP_group"), 
      OSPAR = c(id, "species", "matrix", "sex", "metoa", "AMAP_group"),
      c(id, "species", "matrix", "sex", "metoa")
    )
  if (info$compartment == "water" & info$purpose %in% "HELCOM")
    id <- c(id, "filtered")
  timeSeries <- data[id]
  
  
  # simplify relevant elements of matrix, sex and metoa for level6 and level7, 
  # with a view to then 
  # identifying each time series

  if (info$compartment == "biota") {
    timeSeries <- mutate(
      timeSeries,

      sex = as.character(.data$sex),
      sex = if_else(.data$determinand %in% "EROD", .data$sex, NA_character_),
      sex = factor(sex),
      
      metoa = as.character(metoa),
      .group = ctsm_get_info(
        "determinand", .data$determinand, "group", "biota", sep = "_"
      ),
      metoa = if_else(.group %in% "Metabolites", .data$metoa, NA_character_),
      .group = NULL,
      metoa = factor(metoa)
    )
 }
    
 
  # create seriesID column in data, filled with the concatenated values of timeSeries and 
  # reorder variables

  assign("wk.timeSeries", timeSeries, pos = 1)
    
  data$seriesID <- factor(do.call("pasteOmitNA", timeSeries))

  data <- data[c("seriesID", setdiff(names(data), "seriesID"))]
  
  
  # pick up new.basis and new.unit for each timeSeries 
  
  timeSeries$basis <- data$new.basis
  timeSeries$unit <- data$new.unit
  
  data$new.basis <- data$new.unit <- NULL
  

  # timeSeries is now the unique rows of timeSeries
  
  timeSeries <- droplevels(unique(timeSeries))
  
  id <- setdiff(names(timeSeries), c("basis", "unit"))

  rownames(timeSeries) <- do.call("pasteOmitNA", timeSeries[id])

  
  # change timeSeries output columns to fit the levels of the xml requirements
  
  timeSeries <- changeToLevelsForXML(timeSeries, info$compartment, info$purpose)


  # remove any series that don't have data in recent years
  
  is_recent_data <- (data$year %in% info$recent_years) & !is.na(data$concentration)
  
  id <- data$seriesID[is_recent_data]
  id <- unique(id)
  
  data <- data[data$seriesID %in% id, ]
  data <- droplevels(data)
  
  ok <- row.names(timeSeries) %in% id
  timeSeries <- timeSeries[ok, ]
  timeSeries <- droplevels(timeSeries)
  
  # remove unused stations from the station dictionary

  ok <- row.names(station_dictionary) %in% as.character(timeSeries$station) 
  
  station_dictionary <- droplevels(station_dictionary[ok, ])
  
  
  list(data = data, stations = station_dictionary, timeSeries = timeSeries)
}


changeToLevelsForXML <- function(timeSeries, compartment, purpose) {
  
  # set up columns for populating
    
  timeSeries[c("level6element", "level6name", "level7element", "level7name")] <- NA
    
  if (compartment %in% c("sediment", "water"))
    return(timeSeries)

  if (purpose %in% "AMAP") {
    timeSeries <- within(timeSeries, {
      level6element <- "matrix"
      level6name <- as.character(matrix)
      level7element <- "AMAP_group"
      level7name <- as.character(AMAP_group)
    })
    return(timeSeries)
  }
  
  
  group <- ctsm_get_info(
    "determinand", timeSeries$determinand, "group", compartment, sep = "_"
  )
  
  id <- timeSeries$determinand %in% "EROD"
  if (any(id))
    timeSeries[id, ] <- within(timeSeries[id, ], {
      level6element <- "sex"
      level6name <- as.character(sex)
      level7element <- "matrix"
      level7name <- as.character(matrix)
    })

  id <- timeSeries$determinand %in% "ACHE"
  if (any(id))
    timeSeries[id, ] <- within(timeSeries[id, ], {
      level6element <- "matrix"
      level6name <- as.character(matrix)
    })
  
  id <- group %in% "Metabolites"
  if (any(id))
    timeSeries[id, ] <- within(timeSeries[id, ], {
      level6element <- "METOA"
      level6name <- as.character(metoa)
    })
  
  id <- !group %in% c("Effects", "Imposex", "Metabolites")
  if (any(id))
    timeSeries[id, ] <- within(timeSeries[id, ], {
      level6element <- "matrix"
      level6name <- as.character(matrix)
    })

  id <- !group %in% c("Effects", "Imposex", "Metabolites")
  if (any(id) & purpose %in% "OSPAR") 
    timeSeries[id, ] <- within(timeSeries[id, ], {
      level7element <- "AMAP_group"
      level7name <- as.character(AMAP_group)
    })

  timeSeries
}





ctsm.clean.stations <- function(purpose, ...) {
  cat("\nCleaning station dictionary\n")
  do.call(paste("ctsm.clean.stations", purpose, sep = "."), list(...))
}


ctsm.clean.stations.OSPAR <- function(stations, compartment) {

  # restrict to OSPAR stations used for temporal monitoring 
  # also includes grouped stations: Assessment Grouping for OSPAR MIME

  stations <- filter(
    stations, 
    grepl("OSPAR", .data$programGovernance) & grepl("T", .data$PURPM)
  )

  # and stations that are used for contaminants or biological effects
  
  stations <- filter(
    stations, 
    switch(
      compartment, 
      biota = grepl("CF|EF", .data$dataType), 
      sediment = grepl("CS", .data$dataType), 
      water = grepl("CW", .data$dataType)
    )
  )
  

  # ensure country is consistent
  
  stations <- mutate(stations, country = str_to_title(.data$country))
  
    
  # ensure no backward slashes in station and replace them in name
  
  not_ok <- grepl("\\", stations$station, fixed = TRUE)
  if (any(not_ok))
    stop("backward slashes present in station variable")
  
  stations <- mutate(stations, name = gsub("\\", "/", .data$name, fixed = TRUE))

  
  # create station ID and remove stations that have been replaced
  
  stations <- stations %>% 
    unite(stationID, .data$country, .data$station, remove = FALSE) %>% 
    filter(is.na(.data$replacedBy))
  
  
  # check whether any remaining duplicated stations 
  # if present, select most recent of these
  # crude method to sort by station, startYear and endYear (since NAs are sorted last)
  # and then reverse the ordering of the whole data frame so it works with ctsm.check
  
  stations <- stations %>% 
    arrange(.data$stationID, .data$startYear, .data$endYear) %>% 
    arrange(desc(row_number()))
  
  stations <- ctsm.check(
    stations, stationID, action = "delete.dups",
    message = "Duplicate stations - first one selected",
    fileName = "duplicate stations", merge.stations = FALSE)
  
  
  # identify stations outside the OSPAR region; only a warning because need to
  # retain this information until we have deleted the data that are associated with these stations; 
  # should probably do this in a single function call later on, but low priority

  ctsm.check(
    stations, is.na(OSPAR_region) | is.na(OSPAR_subregion), action = "warning", 
    message = "Stations outside OSPAR region", 
    fileName = "stations outside area", merge.stations = FALSE)
    

  # check all regions are within the appropriate OSPARregion - can sometimes go wrong due to 
  # local shape file errors
  
  ok <- local({
    id1 <- stations$OSPAR_region
    id2 <- info.regions[["OSPAR"]][stations$OSPAR_subregion, "OSPAR_region"]
    (is.na(id1) & is.na(id2)) | id1 == id2
  })
  
  ctsm.check(
    stations, !ok, action = "warning", 
    message = "OSPAR subregion in wrong OSPAR region", 
    fileName = "Region information incorrect", merge.stations = FALSE)
  

  # tidy up output
  
  stations <- arrange(
    stations, 
    .data$OSPAR_region, 
    .data$OSPAR_subregion, 
    .data$country, 
    .data$station
  )  
    
  stations <- column_to_rownames(stations, "stationID")

  col_id <- c(
    "OSPAR_region", "OSPAR_subregion", "country", "station", "code", "name", 
    "latitude", "longitude", "offshore", "MSTAT", "WLTYP"
  )
  stations <- stations[col_id]

  stations
}

ctsm.clean.stations.AMAP <- function(stations, compartment) {
  
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
      compartment, 
      biota = grepl("CF", .data$dataType) | grepl("EF", .data$dataType), 
      sediment = grepl("CS", .data$dataType), 
      water = grepl("CW", .data$dataType)
    )
  )  
  
  # assume missing OSPAR region information corresponds to AMAP stations in 'region 0'
  
  # add in regional information for East AMAP area
  
  stations <- dplyr::mutate_at(stations, c("OSPARregion", "region"), as.character)
  
  stations <- within(stations, {
    OSPARregion[is.na(OSPARregion)] <- "0"
    region[OSPARregion %in% "0"] <- "East AMAP area"
  })
  
  stations <- dplyr::mutate_at(stations, c("OSPARregion", "region"), as.factor)
  
  
  # create station ID and remove stations that have been replaced
  
  stations <- stations %>% 
    unite(stationID, .data$country, .data$station, remove = FALSE) %>% 
    filter(is.na(.data$replacedBy))
  
  
  # check whether any remaining duplicated stations 
  # if present, select most recent of these
  # crude method to sort by station, startYear and endYear (since NAs are sorted last)
  # and then reverse the ordering of the whole data frame so it works with ctsm.check
  
  stations <- stations %>% 
    arrange(.data$stationID, .data$startYear, .data$endYear) %>% 
    arrange(desc(row_number()))

  stations <- ctsm.check(
    stations, stationID, action = "delete.dups", 
    message = "Duplicate stations - first one selected", 
    fileName = "duplicate stations", merge.stations = FALSE)
  
  
  # identify stations outside the AMAP region; only a warning because need to
  # retain this information until we have deleted the data that are associated with these stations; 
  # should probably do this in a single function call later on, but low priority
  
  ctsm.check(
    stations, !(OSPARregion %in% c("0", "1")) | is.na(region), action = "warning", 
    message = "Stations outside AMAP region", 
    fileName = "stations outside area", merge.stations = FALSE)
  
  
  # check all regions (within the OSPAR area) are within the appropriate OSPARregion
  # (can sometimes go wrong due to local shape file errors)
  
  ok <- local({
    id1 <- as.character(stations$OSPARregion)
    region <- as.character(stations$region)
    id2 <- info.regions[["OSPAR"]][region, "OSPARregion"]
    id1 == "0" | (is.na(id1) & is.na(id2)) | id1 == id2 
  })
  
  ctsm.check(
    stations, !ok, action = "warning", 
    message = "Region in wrong OSPAR region", 
    fileName = "Region information incorrect", merge.stations = FALSE)
  
  
  # tidy up output
  
  stations <- stations %>% 
    droplevels() %>% 
    arrange(.data$OSPARregion, .data$region, .data$country, .data$station) %>% 
    column_to_rownames("stationID")
  
  col_id <- c(
    "OSPARregion", "region", "country", "station", "code", "name", "latitude", "longitude", 
    "offshore", "MSTAT", "WLTYP", "ICES_ecoregion"
  )
  stations <- stations[col_id]
  
  stations
}

ctsm.clean.stations.HELCOM <- function(stations, compartment) {
  
  # restrict to stations that have a HELCOM region
  
  stations <- drop_na(stations, .data$HELCOM_subbasin)
  
  
  # ensure country is consistent
  
  stations <- mutate(stations, country = str_to_title(.data$country))
  

  # ensure no backward slashes in station and replace them in name
  
  not_ok <- grepl("\\", stations$station, fixed = TRUE)
  if (any(not_ok))
    stop("backward slashes present in station variable")
  
  
  # create station ID and remove stations that have been replaced
  
  stations <- unite(
    stations, 
    stationID, 
    .data$country, 
    .data$station, 
    remove = FALSE
  ) 
  
  stations <- filter(stations, is.na(.data$replacedBy))
  

  # check whether any remaining duplicated stations 
  # if present, select most recent of these
  # crude method to sort by station, startYear and endYear (since NAs are sorted last)
  # and then reverse the ordering of the whole data frame so it works with ctsm.check

  stations <- stations %>% 
    arrange(.data$stationID, .data$startYear, .data$endYear) %>% 
    arrange(desc(row_number()))
  
  stations <- ctsm.check(
    stations, stationID, action = "delete.dups",
    message = "Duplicate stations - first one selected",
    fileName = "duplicate stations", merge.stations = FALSE)


  # tidy up output

  stations <- arrange(
    stations,
    .data$HELCOM_subbasin, 
    .data$HELCOM_L3, 
    .data$HELCOM_L4, 
    .data$country, 
    .data$station
  ) 
  
  stations <- column_to_rownames(stations, "stationID")
  
  col_id <- c(
    "HELCOM_subbasin", "HELCOM_L3", "HELCOM_L4", 
    "country", "station", "code", "name", "latitude", "longitude",
    "MSTAT", "WLTYP"
  )
  stations <- stations[col_id]
  
  stations  
}


ctsm.clean.contaminants <- function(purpose, ...) {
  cat("\nCleaning contaminant and biological effects data\n")
  data <- do.call(paste("ctsm.clean.contaminants", purpose, sep = "."), list(...))
  
  id <- c(
    "station", "year", "date", "time", "latitude", "longitude", "depth", "species", 
    "filtered", "sampleID", "sample", "sub.sample", 
    "replicate", "upload", "noinp", "sex", 
    "determinand", "pargroup", "matrix", "basis", "metoa", "metcx", "unit", 
    "value", "qflag", "limit_detection", "limit_quantification", "uncertainty", 
    "methodUncertainty", "alabo", "qalink", "AMAP_group"
  )
  
  data[intersect(id, names(data))]
}

ctsm.clean.contaminants.OSPAR <- function(data, compartment, ...) {

  # check whether submitted station has matched to station correctly
  # for Denmark, France, Ireland, Norway, Spain (2005 onwards), Sweden, UK
   
  # have to use sd_name because there is then the grouped station amalgamation

  odd <- with(data, {
    country %in% c("Denmark", "France", "Ireland", "Norway", "Sweden", "United Kingdom") | 
      (country %in% "Spain" & year > 2005)
  })
  
  odd <- odd & with(data, !is.na(submitted.station) & is.na(sd_name))
  
  ctsm.check(
    data, odd, action = "warning", 
    message = "Submitted.station not recognised by dictionary", 
    fileName = "submitted station not recognised", merge.stations = FALSE)
  

  # drop data with no stations  

  cat("   Dropping data with no stations\n")  

  ok <- with(data, !is.na(station_name))
  data <- droplevels(data[ok, ])
  
 
  # create new station that is a combination of country and station - corresponds 
  # to the row.names of the station dictionary - and unique sampleID 

  data <- data %>% 
    unite(station, .data$country, .data$station_name, remove = FALSE) %>%
    mutate(
      station = factor(station),
      sampleID = switch(
        compartment, 
        sediment = factor(sample, labels = ""), 
        biota = factor(sub.sample, labels = ""), 
        water = factor(sample, labels = "")
      )
    )
  
  # drop sample and sub.sample as replicate gives emough information 
  
  data <- data[setdiff(names(data), c("sample", "sub.sample"))]
  
  data
}

ctsm.clean.contaminants.AMAP <- ctsm.clean.contaminants.OSPAR

ctsm.clean.contaminants.HELCOM <- function(data, compartment, ...) {
  
  # check whether submitted station has matched to station correctly for Denmark and Sweden
      
  # (differs from OSPAR because extraction doesn't yet have grouped station amalgamation)
  
  odd <- data$country %in% c("Denmark", "Sweden") 

  odd <- odd & with(data, !is.na(submitted.station) & is.na(station_name))
  
  ctsm.check(
    data, odd, action = "warning", 
    message = "Submitted.station not recognised by dictionary", 
    fileName = "submitted station not in dictionary", merge.stations = FALSE)
  
  
  # drop data with no stations  
  
  cat("   Dropping data with no stations\n")  
  
  ok <- with(data, !is.na(station_name))
  data <- droplevels(data[ok, ])
  
  
  # create new station that is a combination of country and station - corresponds 
  # to the row.names of the station dictionary - and unique sampleID 
  
  data <- data %>% 
    unite(station, .data$country, .data$station_name, remove = FALSE) %>%
    mutate(
      station = factor(station),
      sampleID = switch(
        compartment, 
        sediment = factor(sample, labels = ""), 
        biota = factor(sub.sample, labels = ""), 
        water = factor(sample, labels = "")
      )
    )
  
  # drop sample and sub.sample as replicate gives emough information 
  
  data <- data[setdiff(names(data), c("sample", "sub.sample"))]
  
  data
}






ctsm_clean_QA <- function(purpose, ...) {
  cat("\nCleaning QA data\n")
  do.call(paste("ctsm_clean_QA", purpose, sep = "_"), list(...))
}

ctsm_clean_QA_OSPAR <- function(QA, data, compartment) {
  
  # filter QA by data_type 
  
  id <- switch(compartment, biota = "CF", sediment = "CS", water = "CW")
  
  QA <- dplyr::filter(QA, .data$data_type %in% id)

  
  # retain useful variables  

  QA <- QA[c("alabo", "year", "determinand", "metcx", "metoa", "qalink")]
  
  
  # relabel organotins (which have been relabelled in adjustment file)
  
  QA <- ctsm_TBT_convert(QA, action = "relabel")
  

  # drop crm data with no useful information
  
  QA <- dplyr::filter(QA, !(is.na(.data$metcx) & is.na(.data$metoa)))
  

  # retain unique information

  QA <- unique(QA)


  # qalink and determinand should be unique combinations 
  
  QA <- dplyr::arrange(QA, dplyr::across(c("qalink", "determinand", "metcx", "metoa")))
  
  QA <- ctsm.check(
    QA, 
    paste(qalink, determinand), 
    action = "delete.dups",
    message = "Conflicting QA information - first one selected",
    fileName = "conflicting QA", 
    merge.stations = FALSE
  )
  
  ctsm_link_QA_OSPAR(QA, data, compartment)
}

ctsm_clean_QA_HELCOM <- ctsm_clean_QA_OSPAR

ctsm_clean_QA_AMAP <- ctsm_clean_QA_OSPAR



ctsm_link_QA <- function(purpose, ...)
  do.call(paste("ctsm_link_QA", purpose, sep = "."), list(...))

ctsm_link_QA_OSPAR <- function(QA, data, compartment) {
  
  # create unique qaID to link data and QA files 
  # should be able to use qalink and determinand, but alabo and year give
  # useful extra information

  var_id <- c("determinand", "qalink", "alabo", "year")
  
  data <- tidyr::unite(data, "qaID", dplyr::all_of(var_id), remove = FALSE)
  
  QA <- tidyr::unite(QA, "qaID", dplyr::all_of(var_id), remove = FALSE)
  
  
  # restrict QA to those values found in data
  
  QA <- dplyr::filter(QA, .data$qaID %in% data$qaID)
  

  # get digestion method based on metcx - only needed for sediment

  # if (compartment == "sediment") {
  #   QA <- ctsm.link.QA.digestion(QA, compartment)
  # }


  # tidy up
  
  QA <- tibble::column_to_rownames(QA, "qaID")

  list(QA = QA, data = data)
}

ctsm_link_QA_AMAP <- ctsm_link_QA_OSPAR

ctsm_link_QA_HELCOM <- ctsm_link_QA_OSPAR



ctsm_get_digestion <- function(data) {

  # import_functions.R
  # get digestion method for metals based on metcx
  
  data$digestion <- rep(NA_character_)


  # only apply ctsm_get_info on metals data to avoid all the potential issues
  # associated with organics submissions
  
  is_metal <- data$pargroup %in% "I-MET"  
  
  data[is_metal, "digestion"] <- ctsm_get_info(
    "methodExtraction",
    data[is_metal, "metcx"], 
    "digestion", 
    na_action = "input_ok"
  )
  

  # trap mercury digestion wth metcx of NON  
  
  id <- data$metoa %in% "AAS-CV" & data$metcx %in% "NON"
  if (any(id)) {
    message("   Warning: Ad-hoc way of dealing with mercury digestion that has metcx of NON")
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

  ctsm.check(
    data, not_ok, action = "warning", 
    message = "Missing digestion methods", 
    fileName = "digestion_errors", merge.stations = FALSE)

  data$digestion  
}




ctsm.initialise.oddities <- function(path) {

  # create oddity directories and backup if necessary
  
  if (!(path %in% list.dirs(full.names = FALSE, recursive = FALSE))) 
    dir.create(path)

  compartment <- evalq(info$compartment, sys.frame(1))  
  output <- file.path(path, compartment)

  if (!(output %in% list.dirs(full.names = FALSE))) {
    dir.create(output)
    cat("\nInfo: oddities will be written to '", output, "'\n", sep="")
    return(output)
  } 

  backup <- file.path(path, paste(compartment, "backup", sep = "_"))
  
  if (!(backup %in% list.dirs(full.names = FALSE))) dir.create(backup) 
      
  cat("\nInfo: oddities will be written to '", output, "' with previous oddities\n", 
      "      backed up to '", backup, "'\n", sep="")
  
  old.files <- dir(output, full.names = TRUE)
  file.copy(from = old.files, to = backup, overwrite = TRUE)
  file.remove(old.files)
    
  output  
}


ctsm.imposex.check.femalepop <- function(data) {

  # makes missing values == 100% if sex == F
  
  if (!any(data$group %in% "Imposex"))
    return(data)

  impID <- data$group %in% "Imposex"
  
  replaceID <- impID & is.na(data[["%FEMALEPOP"]]) & data$sex == "F"
  data[replaceID, "%FEMALEPOP"] <- 100

  missingID <- impID & is.na(data[["%FEMALEPOP"]])
  ctsm.check(
    data[impID, ], missingID, action = "warning", 
    message = "Missing FEMALEPOP values", fileName = "imposex missing FEMALEPOP")
  
  return (data)
}


ctsm_check_determinands <- function(compartment, data, determinands, control = NULL) {

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
  

  # check all determinands are recognised by info.determinands
  # note all auxiliaries are guaranteed to be info.determinands
  
  id <- c(determinands, get_control_dets(control))
  
  id <- unique(id)
  
  ok <- id %in% row.names(info.determinand)
  if (!all(ok)) 
    stop('Not found in determinand information file: ', paste(id[!ok], collapse = ", "))
 

  # simplify determinands and determinand.control so they only contain values that are
  # required
  
  id <- unique(data$determinand)
  
  if (!is.null(control))
    id <- c(id, names(control))
  
  determinands <- intersect(determinands, id)

  
  auxiliary <- ctsm_get_auxiliary(determinands, compartment)
  
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
  # matrix and, if so, delete drop - note that ctsm.check doesn't do the
  # deleting because it isn't necessarily the first instance that is retained
  
  ID <- with(data, paste(sampleID, matrix))
  
  dropID <- data$determinand %in% drop 
  keepID <- data$determinand %in% keep
  
  dups <- ID %in% intersect(ID[dropID], ID[keepID])

  if (printDuplicates) {
    
    dropTxt <- paste(drop, collapse = ", ")
    
    dupsID <- dups & (dropID | keepID)
    
    ctsm.check(
      data, dupsID, action = "warning",  
      message = paste(keep, "and", dropTxt, "submitted in same sample - deleting", dropTxt, 
                      "data"), 
      fileName = paste("determinand link", keep), ...)
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
  
  data <- within(data, {
    determinand <- as.character(determinand)
    determinand[determinand %in% drop] <- keep
    determinand <- factor(determinand)
  })
  
  data
}  


determinand.link.imposex <- function(data, keep, drop, ...) {
  
  stopifnot(length(keep) == 1, length(drop) == 1)

  detID <- c(keep, drop)
  
  # for imposex, indices and stages aren't linked by sampleID, but by visit
  # will assume, for simplicity, that only one visit per year
  
  visitID <- with(data, paste(station, year))
  
  
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
  
  ctsm.check(
    data, dups, action = "warning",  
    message = paste("inconsistent", keep, "and", drop, "submitted in same year"), 
    fileName = paste("determinand link", keep), ...)
  
  # delete indices submitted in same visit as individual data (whether consistent or not)

  dups <- tapply(data$determinand, visitID, function(x) all(detID %in% x))
  dups <- visitID %in% names(dups)[dups]
  
  data <- data[!(dups & data$determinand %in% drop), ]
    
  
  # relabel any remaining indices as stage data
    
  if (any(data$determinand %in% drop)) 
    message("   Data submitted as ", drop, "relabelled as ", keep)
  else return(data)
  
  data <- within(data, {
    determinand <- as.character(determinand)
    determinand[determinand %in% drop] <- keep
    determinand <- factor(determinand)
  })
  
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

  data <- within(data, {
    determinand <- as.character(determinand)
    determinand[determinand %in% drop] <- keep
    determinand <- factor(determinand)
  })
  
  data
})  


determinand.link.sum <- function(data, keep, drop, ...) {
  
  stopifnot(length(keep) == 1, length(drop) > 1)
  
  if (!any(data$determinand %in% drop)) 
    return(data)


  # identify samples with drop and not keep, which are the ones that will be summed
  # if keep already exists, then don't need to do anything
  # don't delete drop data because might want to assess them individually
  
  ID <- with(data, paste(sampleID, matrix))
  
  dropID <- data$determinand %in% drop 
  keepID <- data$determinand %in% keep
  
  sum_ID <- ID %in% setdiff(ID[dropID], ID[keepID])

  if (length(sum_ID) == 0)
    return(data)
  
  
  dropTxt <- paste(drop, collapse = ", ")
  cat("   Data submitted as", dropTxt, "summed to give", keep, fill = TRUE)


  # get relevant sample matrix combinations
  
  data <- split(data, with(data, determinand %in% drop & sum_ID))
  
  ID <- with(data[["TRUE"]], paste(sampleID, matrix))

  summed_data <- by(data[["TRUE"]], ID, function(x) {
    
    # check all bases are the same 

    stopifnot(n_distinct(x$basis) == 1)
    
    if (!all(drop %in% x$determinand)) return(NULL)

    
    # adjust values if units vary
    # ideally use unit in info.determinand, but makes it more awkward because
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
    # ensures a sensible qaID, metoa, etc.
    
    out <- x[which.max(x$value), ]
    
    out$determinand <- keep
    
    # sum value and limit_detection, make it a less-than if all are less-thans, and take 
    # proportional uncertainty from maximum value (for which uncertainty is reported)
    
    out$value <- sum(x$value)
    out$limit_detection <- sum(x$limit_detection)
    out$limit_quantification <- sum(x$limit_quantification)
    
    if ("" %in% x$qflag)
      out$qflag <- ""
    else if (n_distinct(x$qflag) == 1) 
      out$qflag <- unique(x$qflag) 
    else 
      out$qflag <- "<"

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
  
  ID <- with(data, paste(sampleID, matrix))
  
  dropID <- data$determinand %in% drop 
  keepID <- data$determinand %in% keep
  
  sum_ID <- ID %in% setdiff(ID[dropID], ID[keepID])
  
  if (length(sum_ID) == 0)
    return(data)
  
  dropTxt <- paste(drop, collapse = ", ")
  cat("   Data submitted as", dropTxt, "summed to give", keep, fill = TRUE)
  
  
  # get relevant sample matrix combinations
  
  data <- split(data, with(data, determinand %in% drop))
  
  ID <- with(data[["TRUE"]], paste(sampleID, matrix))
  
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
    # ensures a sensible qaID, metoa, etc.
    
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
    
    if ("" %in% x$qflag)
      out$qflag <- ""
    else if (n_distinct(x$qflag) == 1) 
      out$qflag <- unique(x$qflag) 
    else 
      out$qflag <- "<"
    out$qflag <- if(all(x$qflag %in% "<")) "<" else ""
    
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

  
ctsm_normalise_sediment <- function(data, station_dictionary, control) {
  
  # normalises sediment concentrations
  
  # method supplied by control
  
  ctsm_normalise_default <- list(
    metals = list(method = "pivot", normaliser = "AL", extra = NULL), 
    organics = list(method = "simple", normaliser = "CORG", value = 2.5), 
    exclude = NULL
  )
  
  control <- modifyList(ctsm_normalise_default, control)
  
  
  # save non-normalised concentrations and qflags for plotting purposes later on
  # also save uncertainties just in case
  
  data <- mutate(
    data, 
    concOriginal = .data$concentration,     
    qflagOriginal = .data$qflag,
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
    exclude_id <- row.names(station_dictionary)[exclude_id]
    exclude_id <- as.character(data$station) %in% exclude_id
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
            "determinand", normaliser, "unit", "sediment", sep = "_"
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
        
        pivot <- info.pivotValues
        pivot <- pivot[pivot$determinand %in% c(as.character(data$determinand), normaliser), ]
        rownames(pivot) <- with(pivot, paste(determinand, digestion))
        pivot <- droplevels(pivot)
        
        
        # get digestion for both contaminant and normaliser (nn indicates not required)
        # if missing digestion for normaliser, assume the strongest method for the 
        #   same sample matrix combination
        # might both be missing if QA does not contain the metcx code 
        
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
          
          linkID <- with(data, paste(sampleID, matrix))
          
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
          station_id <- row.names(station_dictionary)[station_id]
          id <- as.character(data$station) %in% station_id
          
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
      
      notOK <- getNdata("qflag") %in% c("<", "D", "Q")
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


ctsm_normalise_sediment_HELCOM <- function(data, station_dictionary, control) {
  
  # normalises sediment concentrations
  
  # method supplied by control
  
  ctsm_normalise_default <- list(
    metals = list(method = "pivot", normaliser = "AL", extra = NULL), 
    copper = list(method = "hybrid", normaliser = "CORG", value = 5),
    organics = list(method = "simple", normaliser = "CORG", value = 5), 
    exclude = NULL
  )
  
  control <- modifyList(ctsm_normalise_default, control)
  
  
  # save non-normalised concentrations and qflags for plotting purposes later on
  # also save uncertainties just in case
  
  data <- mutate(
    data, 
    concOriginal = .data$concentration,     
    qflagOriginal = .data$qflag,
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
    exclude_id <- row.names(station_dictionary)[exclude_id]
    exclude_id <- as.character(data$station) %in% exclude_id
  } else {
    exclude_id <- FALSE
  }
  
  if (any(exclude_id)) {
    excluded_data <- data[exclude_id, ]
    data <- data[!exclude_id, ]
  }
  
  
  # make ad-hoc change to deal with LOIGN
  # must undo at the end of the code

  data <- mutate(
    data, 
    .tmp = CORG,
    .tmp.qflag = CORG.qflag,
    .tmp.uncertainty = CORG.uncertainty,
    CORG = if_else(is.na(.tmp), 0.35 * LOIGN, CORG),
    CORG.qflag = if_else(
      is.na(.tmp), 
      as.character(LOIGN.qflag), 
      as.character(CORG.qflag)
    ),
    CORG.qflag = factor(CORG.qflag),
    CORG.uncertainty = if_else(is.na(.tmp), 0.35 * LOIGN.uncertainty, CORG.uncertainty)
  )

  
  # split into metals (CD, PB), copper and organics and then normalise each 
  # with AL, CORG (LOIGN) and CORG (LOIGN) respectively#
  
  groupID <- case_when(
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
            "determinand", normaliser, "unit", "sediment", sep = "_"
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
        
        pivot <- info.pivotValues
        pivot <- pivot[pivot$determinand %in% c(as.character(data$determinand), normaliser), ]
        rownames(pivot) <- with(pivot, paste(determinand, digestion))
        pivot <- droplevels(pivot)
        
        
        # get digestion for both contaminant and normaliser (nn indicates not required)
        # if missing digestion for normaliser, assume the strongest method for the 
        #   same sample matrix combination
        # might both be missing if QA does not contain the metcx code 
        
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
          
          linkID <- with(data, paste(sampleID, matrix))
          
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
        
        pivot <- info.pivotValues
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
      
      notOK <- getNdata("qflag") %in% c("<", "D", "Q")
      if (any(notOK)) {
        message('   Removing sediment data where normaliser is a less than')
        concentration[notOK] <- NA
      }
      
      data$concentration <- concentration
      data$uncertainty <- uncertainty
      
      data
    })
  
  data <- unsplit(data, groupID)


  data <- mutate(
    data, 
    CORG = .tmp,
    CORG.qflag = .tmp.qflag,
    CORG.uncertainty = .tmp.uncertainty,
    .tmp = NULL,
    .tmp.qflag = NULL,
    .tmp.uncertainty = NULL
  )

  if (any(exclude_id)) {
    data <- dplyr::bind_rows(data, excluded_data)
  }
  
  data
}


ctsm_normalise_biota_HELCOM <- function(data, station_dictionary, control) {
  
  # normalises biota concentrations on a lipid weight basis to 5% lipid
  
  # method supplied by control
  
  ctsm_normalise_default <- list(
    lipid = list(method = "simple", value = 5), 
    other = list(method = "none"), 
    exclude = NULL
  )
  
  control <- modifyList(ctsm_normalise_default, control)
  
  
  # save non-normalised concentrations and qflags for plotting purposes later on
  # also save uncertainties just in case
  
  data <- mutate(
    data, 
    concOriginal = .data$concentration,     
    qflagOriginal = .data$qflag,
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
    exclude_id <- row.names(station_dictionary)[exclude_id]
    exclude_id <- as.character(data$station) %in% exclude_id
  } else {
    exclude_id <- FALSE
  }
  
  if (any(exclude_id)) {
    excluded_data <- data[exclude_id, ]
    data <- data[!exclude_id, ]
  }
  
  
  # split into contaminants on a lipid weight basis and others
  
  groupID <- if_else(data$new.basis %in% "L", "lipid", "other")
  
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

      data$concentration <- data$concentration * control$value / 100
      data$uncertainty <- data$uncertainty * control$value / 100
      
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
  
  Css <- if_else(ok, (Cm - Cx) * (Nss - Nx) / (Nm - Nx) + Cx, NA_real_)
  
  var_Css <- if_else(ok, ((Cm - Cx) / (Nm - Nx))^2, NA_real_)
  var_Css <- (Nss - Nx)^2 * (var_Cm + var_Nm * var_Css) + (Nss - Nm)^2 * (var_Cx + var_Nx * var_Css)
  var_Css <- if_else(ok, var_Css / ((Nm - Nx)^2), NA_real_)
  
  data.frame(Cm, Nm, Css, var_Css)
}



ctsm.estimate.uncertainty <- function(data, response_id, compartment) {

  # estimating missing uncertainties
  
  # only return uncertainty so can modify data object at will
  
  # check response_id is valid (could be "concentration" or an auxiliary variable)
  # check compartment is valid and not a variable in data
  
  stopifnot(
    is.character(response_id),
    length(response_id) == 1,
    response_id %in% names(data),
    
    is.character(compartment),
    length(compartment) == 1
  )
  

  # standardise variable names when response_id is an auxiliary variable

  if (response_id != "concentration") {
    data <- mutate(data, determinand = rep(response_id, nrow(data)))

    var1 <- c("uncertainty", "qflag", "limit_detection", "limit_quantification", "qaID")  
    var2  <- paste(response_id, var1, sep = ".")
    
    data[c("concentration", var1)] <- data[c(response_id, var2)]
  }
    

  # is there anything to do?

  missing_id <- is.na(data$uncertainty)
  
  if (!any(missing_id)) 
    return(uncertainty)

  
  # get two components of variance
  
  data$sd_constant <- ctsm_get_info(
    "uncertainty", data$determinand, "sd_constant", compartment, na_action = "output_ok")
  
  data$sd_variable = ctsm_get_info(
    "uncertainty", data$determinand, "sd_variable", compartment, na_action = "output_ok")
  
  
  # adjust sd_constant to correct basis (
  # for biota sd_constant is on a wet weight basis but the sample might be on 
  # a dry or lipid weight

  if (compartment == "biota") {
    data$sd_constant <- convert.basis(
      data$sd_constant, 
      rep("W", length(data$sd_constant)), 
      data$new.basis, 
      data[["DRYWT%"]], 
      data[["DRYWT%.qflag"]], 
      data[["LIPIDWT%"]], 
      data[["LIPIDWT%.qflag"]], 
      exclude = data$group %in% c("Imposex", "Effects", "Metabolites")
    )
  }

  
  # adjust sd_constant to use limit and qflag information where possible 
  # however qflag < is difficult to work with as it is inconsistent with both lod and loq so
  #   take lower of concentration / 3 and sd_constant
  # qflag = "" is also difficult to work with because we can't trust the lod and loq 
  #   option - again take lower of many options
  
  data$sd_constant <- case_when(
    data$qflag %in% "D" ~ data$limit_detection / 3,
    data$qflag %in% "Q" ~ data$limit_quantification / 10,
    data$qflag %in% "<" ~ pmin(data$sd_constant, 
                               data$concentration / 3, 
                               na.rm = TRUE),
    data$qflag %in% ""  ~ pmin(data$sd_constant, 
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
    uncertainty = if_else(
      is.na(.data$uncertainty), .data$estimated_uncrt, .data$uncertainty
    )
  )

  data$uncertainty
}



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
      ! "methodUncertainty" %in% names(data))
    stop("methodUncertainty not provided to allow conversion of uncertainty")
  
  
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
      with(data[id, ], stopifnot(is.na(methodUncertainty) | methodUncertainty %in% "U2"))
    
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


ctsm_TBT_remove_early_data <- function(ctsm_obj, recent_year = min(ctsm_obj$info$recent_years)) {
  
  tin <- c("TBTIN", "TPTIN", "DBTIN", "DPTIN", "MBTIN", "MPTIN")
  cation <- c("TBSN+", "TPSN+", "DBSN+", "DPSN+", "MBSN+", "MPSN+")
  
  id <- with(ctsm_obj$data, determinand %in% c(tin, cation))
  id <- factor(id, levels = c("TRUE", "FALSE"))
  
  data <- split(ctsm_obj$data, id)
  
  tin_data <- data[["TRUE"]]
  
  # drop data with no stations
  
  ok <- with(tin_data, !is.na(station_name))
  tin_data <- tin_data[ok, ]
  
  # drop time series with no data in recent years - avoids sorting out old data that won't get used
  
  wk <- with(tin_data, tapply(year, station_name[, drop = TRUE], max))
  ok <- with(tin_data, station_name %in% names(wk)[wk >= recent_year])
  tin_data <- tin_data[ok, ]
  
  ctsm_obj$data <- rbind(data[["FALSE"]], tin_data)
  
  ctsm_obj
} 


   
