# main import and data processing functions ----

library(readxl)

# read data  ---- 

#' Read HARSAT data
#'
#' Reads in contaminant and effects data, the station dictionary and various
#' reference tables. For data from the ICES webservice, it matches data to
#' stations in the station dictionary.  It also allows the user to set control
#' parameters that dictate the assessment process.
#'
#' @param compartment A string: `"biota"`, `"sediment"` or `"water"`
#' @param purpose A string specifying whether to use the default set up for
#'   `"OSPAR"`, `"HELCOM"`, or `"AMAP"` or to use a customised setup `"custom"`
#' @param contaminants A file reference for the contaminant data
#' @param stations A file reference for the station data
#' @param data_dir The directory where the data files can be found (sometimes
#'   supplied using 'file.path'). Defaults to "."; i.e. the working directory.
#' @param data_format A string specifying whether the data were extracted from
#'   the ICES webservice (`"ICES"` - the default) or are in the simplified
#'   format designed for other data sources (`"external"`).
#' @param info_files A list of files specifying reference tables which override
#'   the defaults. See examples.
#' @param info_dir The directory where the reference tables can be found
#'   (sometimes supplied using 'file.path'). Defaults to "."; i.e. the working
#'   directory
#' @param extraction A date saying when the extraction was made. Optional. This
#'   should be provided according to ISO 8601; for example, 29 February 2024
#'   should be supplied as "2024-02-29". If the contaminant data were extracted
#'   from the ICES webservice and the download file name has not been changed,
#'   the extraction data will be taken from the contaminant file name.
#' @param max_year An integer giving the last monitoring year that should be
#'   included in the assessment. Data from monitoring years after `max_year`
#'   will be deleted. If not specified `max_year` is taken to be the last
#'   monitoring year in the contaminant data file.
#' @param oddity_dir The directory where the 'oddities' will be written
#'   (sometimes supplied using 'file.path'). This directory (and subdirectories)
#'   will be created if it does not already exist.
#' @param control A list of control parameters that override the default values
#'   used to run the assessment. These include the reporting window; the way in
#'   which data are matched to stations following an ICES extraction;
#'   information about reporting regions, and so on. See Details.
#'
#' @returns A list with the following components:
#' * `call` The function call.
#' * `info` A list containing the reference tables and the control parameters.
#' * `data` A data frame containing the contaminant (and effects) data. For
#'   `external` data, this is identical to the input data file apart from some
#'   extra empty columns which have been added. For `ICES` data, some existing
#'   columns have been renamed (otherwise they are untouched) and some
#'   additional columns have been constructed. The key ones of these are:
#'   - `station_code` the code of the station in the station dictionary that
#'   best matches the data
#'   - `station_name` the name of the station
#'   - `species` (biota) the species based on `worms_accepted_name` where
#'   available and `speci_name` otherwise
#'   - `filtration` (water) whether the sample was `filtered` or `unfiltered`
#'   based on `method_pretreatment`
#'   - `retain` a logical indicating whether each record would have been
#'   retained under the previous ICES extraction protocol. For example, `retain`
#'   will be `FALSE` if the vflag entry is `"S"` or suspect. Records for which
#'   `retain == FALSE` are deleted later in `tidy_data`
#' * `stations`
#'
#' @details
#'
#'   ## Control parameters
#'
#'   Many aspects of the assessment process can be controlled through the
#'   parameters stored in `info$control`. This is a list populated with default
#'   values which can then be overwritten, if required, using the `control`
#'   argument.
#'
#'   ## External data
#'
#'   If `data_format = "external"`, a simplified data and station file can
#'   be supplied. These should be .csv files, preferably saved with a UTF-8 
#'   encoding. All missing value should be supplied as empty cells, not as `NA`
#'   or some other code.    
#'   
#'   The data file has one row for each measurement. The mandatory variables 
#'   are: 
#'
#'   * `country`: identifies the source of the data. For international 
#'   assessments, this is typically the country of origin, but for national 
#'   assessments it could be a local monitoring authority. It is read in as a 
#'   character string. 
#'   * `station_code`: the code of the station where the data were collected. 
#'   This can be numeric, but is read in as a character string. 
#'   * `station_name`: the name of the station where the data were collected. 
#'   This is usually more meaningful to a user than `station_code`. 
#'   * `year`: an integer giving the monitoring year.
#'   * `sample`: an identifier that is used to match measurements from the same 
#'   individual (biota), sediment sample or water sample. This can be numeric, 
#'   but will be read in as a character string. Take care over this variable as
#'   if it is incorrectly specified you could lose a lot of data.
#'   * `determinand`: the code identifying the thing measured. This should 
#'   match a record in the determinand reference table. 
#'   * `matrix`: the code identifying the tissue (biota) or size fraction 
#'   (sediment) of the sample. Use ICES codes at present. For water it should 
#'   always be set to `"WT"`. 
#'   * `unit`: the unit of the measurement. Use ICES codes at present.
#'   * `value`: the numeric value of the measurement. 
#'    
#'   The values of `country`, `station_code` and `station_name` (in each row)
#'   must match an entry in the station file. No missing values are allowed in 
#'   the mandatory data variables.
#'  
#'   The station file has one row for each station. The mandatory variables are: 
#'   
#'   * `country`: see data file variables 
#'   * `station_code`: see data file variables
#'   * `station_name`: see data file variables
#'   * `station_latitude`: the nominal latitude of the station (nominal because
#'   samples are rarely taken at exactly the same location each year). This 
#'   should be provided as decimal degrees (or a numeric value based on some 
#'   other projection).
#'   * `station_longitude`: the nominal longitude of the station. This should 
#'   be provided as decimal degress (or a numeric value based on some other 
#'   projection).
#'
#'   The optional variables (including the semi-optional regional information)
#'   will be documented at a later date.
#'   
#' @export
read_data <- function(
  compartment = c("biota", "sediment", "water"), 
  purpose = c("OSPAR", "HELCOM", "AMAP", "custom"), 
  contaminants, 
  stations, 
  data_dir = ".", 
  data_format = c("ICES", "external"),
  info_files = list(),
  info_dir = ".",
  extraction = NULL, 
  max_year = NULL,
  oddity_dir = "oddities",
  control = list()) {

  # import_functions.R

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
        "'extraction' must be a single character string of the form ", 
        "\"yyyy-mm-dd\": e.g. \"",  lubridate::today(), "\"", 
        call. = FALSE
      )
    }
  }  
  
  
  # turn extraction into date object

  if (!is.null(extraction)) {
    extraction <- suppressWarnings(lubridate::ymd(extraction))
    if (is.na(extraction)) {
      stop(
        "extraction not recognised: it must be a valid date of the form ", 
        "\"yyyy-mm-dd\": e.g. \"", lubridate::today(), "\"",
        call. = FALSE
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
  
  cntrl <- control_default(purpose, compartment)
  
  cntrl <- control_modify(cntrl, control)

  if (any(names(cntrl) %in% names(info))) {
    id <- names(cntrl)
    id <- id[id %in% names(info)]
    warning(
      "\n conflict between function arguments and control elements ", 
      "- results may be unexpected:\n ",
      paste(id, collapse = ", "), "\n",
      call. = FALSE, immediate. = TRUE)
  }
  
  info <- append(info, cntrl)
  
  
  # read in reference tables
  
  info <- read_info(info, info_dir, info_files)
  

  # read in station dictionary, contaminant and biological effects data and QA data
  
  stations <- read_stations(stations, data_dir, info)

  data <- read_contaminants(contaminants, data_dir, info)

  if (data_format == "ICES") {

    data <- add_stations(data, stations, info)
    
    stations <- finalise_stations(stations, info)
    
    data <- finalise_data(data, info)

  }
    

  # populate max_year if not set in function call
  
  if (is.null(info$max_year)) {
    info$max_year <- max(data$year)
    cat(
      "\nArgument max_year taken to be the maximum year in the data:", 
      info$max_year,
      "\n"
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
  
  out
}


#' Reads an input file, given a particular encoding, and possibly additional hints
#' 
#' @param file the file name
#' @param header logical (default `TRUE`), whether or not a header line is expected
#' @param sep character (defailt `,`), field separator character
#' @param quote character (default `"`), field quote character
#' @param dec character (default `.`), decimal character
#' @param fill logical (default `TRUE`), if rows have unequal lengtg, additional blank fields are added
#' @param comment.char character (default is empty), if not empty, a comment character
#' @param strip.white local (default `TRUE`), strips leading and trailing white space from fields
#' @param ... any additional parameters to [`utils::read.table`]
#' @return a data frame
safe_read_file <- function(file, header=TRUE, sep=",", quote="\"", dec=".", fill=TRUE, comment.char="", strip.white=TRUE, ...) {

  ## Handle Excel
  extension <- toupper(tools::file_ext(file))
  if (extension == 'XLS' || extension == 'XLSX') {
    result <- readxl::read_excel(file, col_names = TRUE)
    return(result)
  }

  ## Check the encoding, and warn if incompatible -- we allow both UTF-8 and ASCII, 
  ## since UTF-8 is a superset of ASCII, and guess_encoding will only guess one of
  ## them. This only warns, it does not stop a workflow.  
  enc <- readr::guess_encoding(file)
  encodings <- toupper(enc$encoding)
  if (length(intersect(c("ASCII", "UTF-8"), encodings)) == 0) {
    warning(
      "file encoding isn't ASCII or UTF-8: ", 
      file,
      ", please use iconv or save this file as UTF-8"
    )
  }

  table <- utils::read.table(
    file=file,
    header=header,
    sep=sep,
    quote=quote,
    dec=dec,
    fill=fill,
    comment.char=comment.char,
    strip.white=strip.white,
    fileEncoding='utf-8',
    ...
  )
  table
}


control_default <- function(purpose, compartment) {
  
  # import functions
  # sets up default values that control the assessment
  
  # reporting_window is set to 6 to match the MSFD reporting cycle
  # series with no data in the most recent_window are excluded
  
  # region$id are the names or the relevant regions in the data
  # region$names are the names that will be used for reporting 
  
  # region$all is a logical that determines whether all data (and stations)
  # must have an associated region

  # bivalve_spawning_season is a character vector of months when contaminant
  # data for bivalves and gastropds will be deleted because they are in the 
  # spawning season
  
  # use_stage is a logical which determines whether, for biota, stage is used
  # to populate subseries
  
  # relative_uncertainties is a 2-vector giving the range of acceptable 
  # relative uncertainties for log-normally distributed data; the default is
  # to accept relative uncertainties greater than (but not equal) to 1% and 
  # less than (but not equal to) 100%
  
  region <- list()
  
  region$id <- switch(
    purpose, 
    OSPAR = c("ospar_region", "ospar_subregion"),
    HELCOM = c("helcom_subbasin", "helcom_l3", "helcom_l4"),
    NULL
  )

  region$names <- switch(
    purpose, 
    OSPAR = c("OSPAR_region", "OSPAR_subregion"),
    HELCOM = c("HELCOM_subbasin", "HELCOM_L3", "HELCOM_L4"),
    NULL
  )
  
  region$all <- !is.null(region$id)
    
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
  
  use_stage <- FALSE
  
  relative_uncertainty <- c(1, 100)

  add_stations <- switch(
    purpose, 
    OSPAR = list(
      method = "both",
      area = "OSPAR",
      datatype = TRUE,
      temporal = TRUE,
      governance_type = "both",
      governance_id = c("OSPAR", "AMAP"),
      group = TRUE,
      check_coordinates = FALSE
    ),
    HELCOM = list(
      method = "both",
      area = "HELCOM",
      datatype = FALSE,
      temporal = FALSE,
      governance_type = "none",
      governance_id = NULL,
      group = FALSE,
      check_coordinates = FALSE
    ),
    AMAP = list(
      method = "both",
      area = NULL,
      datatype = FALSE,
      temporal = FALSE,
      governance_type = "none",
      governance_id = NULL,
      group = FALSE,
      check_coordinates = FALSE
    ),
    custom = list(
      method = "name",
      area = NULL,
      datatype = FALSE,
      temporal = FALSE,
      governance_type = "none",
      governance_id = NULL,
      group = FALSE,
      check_coordinates = FALSE
    )
  )
    

  list(
    reporting_window = 6L, 
    region = region,
    add_stations = add_stations,
    bivalve_spawning_season = bivalve_spawning_season,
    use_stage = use_stage,
    relative_uncertainty = relative_uncertainty
  )
}


control_modify <- function(control_default, control) {

  # import functions
  # updates default control structure with user specification and does basic 
  # error checking
  
  # region: initialise names if id has been modified and names is not specified
  
  if ("region" %in% names(control)) {
    wk <- control$region
    if ("id" %in% names(wk) && !("names" %in% names(wk))) {
      control$region <- append(wk, list(names = wk$id))
    }
  }

  
  control <- modifyList(control_default, control, keep.null = TRUE)

  if (length(control$region$id) != length(control$region$names)) {
    stop(
      "error in control argument: length of region$id and region$names must be ",
      "identical", 
      call. = FALSE
    )
  }

  
  if (!is.null(control$bivalve_spawning_season)) {
    months <- c(
      "January", "February", "March", "April", "May", "June", "July", "August",
      "September", "October", "November", "December"
    )
    if (!all(control$bivalve_spawning_season %in% months)) {
      stop(
        "error in control argument: invalid months in bivalve_spawning_season", 
        call. = FALSE
      )
    }
  }

  
  if (length(control$relative_uncertainty) != 2L || 
      control$relative_uncertainty[1] < 0 || 
      control$relative_uncertainty[1] > control$relative_uncertainty[2]) {
    stop(
      "error in control argument: invalid range of acceptable relative ", 
      "uncertainties",
      call. = FALSE
    )
  }
  
  control
}


#' Find the path for an information file
#' 
#' Locates a requested information file, searching the
#' information file path. If the requested file cannot be found and
#' `required` is not false, stops with an error.
#' 
#' @param name A string: the name of the file, e.g., `thresholds_biota.csv`
#' @param path A vector of strings, directories to search. The information directory
#'   for the package is automatically searched if we haven't found the 
#'   file anywhere else
#' @returns A string, the absolute path for the file, or `NULL` if the file 
#'   cannot be found anywhere.
locate_information_file <- function(name, path) {

  # No point doing this intelligently. The goal is to find the file with
  # least steps.
  for(directory in path) {
    search <- normalizePath(file.path(directory, name), mustWork = FALSE)
    if(file.exists(search)) {
      cat(paste("Found in path", name, search, "\n"))
      return(search)
    }
  }

  # If we fail to find it, fall back to system.file, which isn't clearly
  # documented but suggests it will return the full path.
  search_file <- file.path("information", name)
  search <- system.file(search_file, package = "harsat", mustWork = FALSE)
  if(file.exists(search)) {
    cat(paste("Found in package", name, search, "\n"))
    return(search)
  }
    warning(paste("Missing file in path:", name, path))
  return(NULL)
}

read_info <- function(info, path, info_files) {
  
  # location: import_functions.R
  # purpose: reads in reference tables, using default files unless overruled 
  #  by info_control
  
  # default reference tables

  ## If the path is a string, make it a vector path
  if(is.character(path)) {
    path <- c(path)
  }

  files <- list(
    determinand = locate_information_file("determinand.csv", path),
    species = locate_information_file("species.csv", path),
    thresholds = locate_information_file(paste0("thresholds_", info$compartment, ".csv"), path)
  )
  
  files$method_extraction <- locate_information_file("method_extraction.csv", path)
  files$pivot_values <- locate_information_file("pivot_values.csv", path)
  files$matrix <- locate_information_file("matrix.csv", path)
  files$imposex <- locate_information_file("imposex.csv", path)
  
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
    files$determinand, info$compartment
  )
  
  info$matrix <- ctsm_read_matrix(files$matrix)
  
  if (info$compartment == "biota") {
    info$species <- ctsm_read_species(files$species)
    info$imposex <- ctsm_read_imposex(files$imposex)
  }

  if (info$compartment == "sediment") {
    info$method_extraction <- ctsm_read_method_extraction(
      files$method_extraction
    ) 
    info$pivot_values <- ctsm_read_pivot_values(files$pivot_values)
  }  
      
  if (!is.null(files$thresholds)) {
    info$thresholds <- ctsm_read_thresholds(
      files$thresholds, info$compartment
    )
  }
  
  info
}


#' Generates and logs a file digest
#' 
#' @description
#' Writes a line on the console containing the file name and its
#' MD5 digest. This can be used to track changes in input files, for
#' reproducibility reasons. MD5 is good enough for this as we aren't 
#' after cryptograohic levels of security, and its cheaper to compute. 
#' 
#' @param infile the input file name
report_file_digest <- function(infile) {
  md5 <- digest::digest(infile, algo='md5')
  cat("MD5 digest for: '", infile, "': '", md5, "'\n", sep = "")
}


#' Reads station dictionary
#'
#' @param file A file reference for the station dictionary.
#' @param data_dir A path to the directory holding the station dictionary.
#'   Defaults to the working directory.
#' @param info A list containing at least the following two elements:
#' * compartment: `"biota"`, `"sediment"` or `"water"`
#' * data_format: `"ICES"` or `"external"`
#' @returns A data frame containing the station dictionary.
read_stations <- function(file, data_dir = ".", info) {

  # import functions

  infile <- file.path(data_dir, file)
  cat("Reading station dictionary from:\n '", infile, "'\n", sep = "")

  
  if (info$data_format == "ICES") {

    var_id <- c(
      station_code = "character",
      station_country = "character", 
      helcom_subbasin = "character", 
      helcom_l3 = "character", 
      helcom_l4 = "character", 
      ices_ecoregion = "character", 
      ospar_region = "character", 
      ospar_subregion = "character", 
      is_amap_area = "logical", 
      is_helcom_area = "logical", 
      is_ospar_area = "logical", 
      station_name = "character", 
      station_longname = "character", 
      station_replacedby = "character",
      station_asmtmimeparent = "character",
      station_activefromdate = "character", 
      station_activeuntildate = "character", 
      station_programgovernance = "character", 
      station_governance = "character", 
      station_purpm = "character",             
      station_latitude = "numeric",
      station_latituderange = "numeric",
      station_longitude = "numeric", 
      station_longituderange = "numeric",
      station_geometry = "character", 
      station_datatype = "character",         
      station_wltyp = "character", 
      station_mstat = "character",             
      station_notes = "character", 
      station_deprecated = "logical"
    )
    
    report_file_digest(infile)
    stations <- safe_read_file(
      infile, 
      sep = "\t",
      quote = "",
      na.strings = c("", "NULL"), 
      strip.white = TRUE,
      colClasses = var_id
    )
    
  }    
  

  if (info$data_format == "external") {  
    
    var_id <- c(
      country = "character",
      station_name = "character",
      station_code = "character",
      station_longname = "character",
      station_latitude = "numeric", 
      station_longitude = "numeric",
      station_type = "character",
      waterbody_type = "character"
    )
    
    if (!is.null(info$region$id)) {
      extra <- rep("character", length(info$region$id))
      names(extra) <- info$region$id
      var_id <- c(var_id, extra)
    }
    
    required <- c(
      "country", "station_name", "station_code", "station_latitude", 
      "station_longitude", info$region$id
    )
    
    
    # check required variables are present in data
    
    report_file_digest(infile)
	stations <- safe_read_file(infile, strip.white = TRUE, nrows = 1)
    
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
    
    stations <- safe_read_file(
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


#' Reads contaminant data 
#'
#' A quick way of reading the contaminant data (which also avoids having to do 
#' any station matching if `info$data_format == "ICES`).
#'
#' @param file A file reference for the contaminant data.
#' @param data_dir A path to the directory holding the contaminant data.
#'   Defaults to the working directory.
#' @param info A list containing at least the following two elements:
#' * compartment: `"biota"`, `"sediment"` or `"water"`
#' * data_format: `"ICES"` or `"external"`
#' @returns A data frame containing the contaminant data.
#' 
#' @export
read_contaminants <- function(file, data_dir = ".", info) {

  # silence non-standard evaluation warnings
  .data <- NULL
  
  # import functions
  # read in contaminant (and biological effects) data
  
  infile <- file.path(data_dir, file)
  cat("\nReading contaminant and effects data from:\n '", 
      infile, "'\n", sep = "")


  if (info$data_format == "ICES") {

    var_id <- c(
      country = "character",
      mprog = "character", 
      helcom_subbasin = "character",     
      helcom_l3 = "character",
      helcom_l4 = "character",
      ices_ecoregion = "character",
      ospar_region = "character", 
      ospar_subregion = "character", 
      is_amap_monitoring = "logical",
      is_helcom_monitoring = "logical", 
      is_medpol_monitoring = "logical", 
      is_ospar_monitoring = "logical",  
      is_amap_area = "logical", 
      is_helcom_area = "logical",
      is_ospar_area = "logical",
      rlabo = "character",
      slabo = "character",
      alabo = "character",               
      statn = "character",
      myear = "integer",
      date = "Date",                
      latitude = "numeric",
      longitude = "numeric",
      dephu = "numeric",               
      dephl = "numeric",
      purpm = "character",
      finfl = "character",
      param = "character",
      pargroup = "character",
      matrx = "character",              
      basis = "character",
      value = "numeric",
      munit = "character",              
      detli = "numeric",
      lmqnt = "numeric",
      uncrt = "numeric",              
      metcu = "character",
      qflag = "character",
      vflag = "character",              
      metoa = "character",
      metcx = "character",
      metpt = "character",
      metst = "character",
      metps = "character",
      metfp = "character",               
      smtyp = "character",
      smpno = "character",
      subno = "character",              
      dcflgs = "character",
      tblanalysisid = "integer",
      tblparamid = "integer",          
      tblsampleid = "character",
      tblspotid = "integer",
      tbluploadid = "integer"
    )

    if (info$compartment == "biota") {

      extra <- c(
        rlist = "character",
        speci = "character",
        speci_name = "character",
        aphiaid = "integer",
        worms_name = "character",
        aphiaid_accepted = "integer",
        worms_accepted_name = "character",
        sexco = "character",
        stage = "character",
        noinp = "integer",
        bulkid = "character",
        tblbioid = "character",
        accessionid = "integer"
      )

      var_id <- c(var_id, extra)

    }

    report_file_digest(infile)
    data <- safe_read_file(       
      infile, 
      sep = "\t",
      # quote = "",             ## Quotes shouldn't be allowed, but disabling them breaks all reads
      na.strings = c("", "NULL"), 
      strip.white = TRUE,
      colClasses = var_id
    )
    
    # add in extra required variable
    
    if (info$compartment == "biota" && info$use_stage) {
      data$subseries <- data$stage
    } else {
      data$subseries <- NA_character_
    }

    
    # convert param (determinand) and munit (unit) to upper / lower case
    
    data <- dplyr::mutate(
      data, 
      param = toupper(.data$param),
      munit = tolower(.data$munit)
    )
    
        
    return(data)
  }  
  

  if (info$data_format == "external") {

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
  

  if (info$data_format == "external") {  

    # check required variables are present in data

    report_file_digest(infile)    
    data <- safe_read_file(infile, strip.white = TRUE, nrows = 1)
  
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

    data <- safe_read_file(
      infile, 
      na.strings = c("", "NULL"),
      strip.white = TRUE,
      colClasses = var_id[ok]
    )
    
  }
  

  # create missing (non-required) variables 
  
  if (info$data_format == "external") {
  
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
      
  }  
  

  
  # additional validations (for external data)
  
  if (info$data_format == "external") { 
  
    # every valid `uncertainty` must have a valid `unit_uncertainty`
    uncertainty_present <- which(complete.cases(data$uncertainty))
    uncertainty_present_valid_units <- 
      data$unit_uncertainty[uncertainty_present] %in% c("%", "U2", "SD")
    if (! all(uncertainty_present_valid_units)) {
      stop(
        "Missing or invalid uncertainty units for specified uncertainty values. ",
        "Please check that all uncertainty values have a valid unit: %, U2, or SD"
      )
    }
  }    
  
  
  # check regional identifiers are in the extraction 

  # if (info$data_format == "ICES_old") {
  #   
  #   pos <- names(data) %in% info$region$id
  #   if (sum(pos) != length(info$region$id)) {
  #     stop("not all regional identifiers are in the data extraction")
  #   }
  # 
  #   names(data)[!pos] <- tolower(names(data)[!pos])     # just in case!
  #   
  # } else 
  
  if (info$data_format == "external") {
    names(data) <- tolower(names(data))  
  }
  
  
  # create more useful names
  # for biota, tblsampleid is the species, tblbioid gives the subsample (individual) 
  # which with matrix gives the unique sample id
  # for sediment, tblsampleid and matrix give the unique sample id
  # to streamline, could make sample = tblbioid (biota) or tblsampleid (sediment)

  
  # variables common across purpose and compartments
  

  # ensure further consistency
  
  data <- dplyr::mutate(
    data, 
    determinand = toupper(.data$determinand),
    unit = tolower(.data$unit)
  )
  
  
  data
}


#' Add stations to contaminant data from an ICES extraction
#'
#' Adds the station name and station code to the contaminant data from an ICES
#' extraction. This is done by either matching the station names submitted with
#' the data to the station dictionary, or by matching the sample coordinates to
#' the station dictionary, or a combination of both.
#'
#' @param data A data frame with the contaminant data from an ICES extraction
#' @param stations A data frame with the ICES station dictionary
#' @param info A HARSAT information list which must contain the elements
#'   `purpose`, `compartment`, and `add_stations`. The latter is a list of
#'   control parameters supplied through `control_default` or `control_modify`
#'   which control how the station matching is achieved. See details.
#'
#' @details `info$add_stations` is a list of control parameters that modify the
#'   station matching process:
#' * method: a string specifying whether the stations are matched by `"name"`,
#'   `"coordinates"`, or `"both"`. If `info$purpose` is `"custom"`, `method` is
#'   restricted to either `"name"` (the default) or `"coordinates"`. If
#'   `info$purpose` is `"OSPAR"`, `"HELCOM"` or `"AMAP"`, then method is set to
#'   `"both"` by default and stations are matched by name or coordinates
#'   according to rules specified by OSPAR, HELCOM or AMAP data providers.
#'   Specifically, stations are matched by name for Denmark, France (biota and
#'   water - all years; sediment 2009 onwards), Ireland, Norway, Portugal, Spain
#'   (2005 onwards), Sweden, The Netherlands (2007 onwards), United Kingdom. All
#'   other stations are matched by coordinates.
#' * area: a vector of strings containing one or more of `"OSPAR"`, `"HELCOM"`
#'   and `"AMAP"`; this restricts the stations to those in the corresponding
#'   convention area(s); NULL matches to all stations in the station dictionary
#' * datatype: a logical specifying whether the stations should be restricted
#'   to those with an appropriate datatype. If `TRUE`, a contaminant measurement
#'   in biota (for example) will only be matched to stations with
#'   `station_datatype` containing the string `"CF"`. Similarly, a biological
#'   effect measurement in biota will only be matched to stations with
#'   `station_datatype` containing the string `"EF"`
#' * temporal: a logical with `TRUE` indicating that stations should be
#'   restricted to those with `station_purpm` containing the string `"T"`
#' * governance_type: a string: `"none"`, `"data"`, `"stations"` or `"both"`.
#'   `"none"` means data and station governance are both ignored. `"data"` means
#'   that matching will be restricted by data governance but not station
#'   governance; for example if `governance_id == c("OSPAR", "AMAP")`, then data
#'   will only be matched to a station if one of `is_ospar_monitoring` and
#'   `is_amap_monitoring` is `TRUE`, with all stations considered regardless of
#'   station governance. `"stations"` mean that matching will be restricted by
#'   station governance but not by data governance; for example if
#'   `governance_id == c("OSPAR", "AMAP")`, then the stations will be restricted
#'   to those where `station_programgovernance` contains either `"OSPAR"` or
#'   `"AMAP"`, with all data considered regardless of data governance. `both`
#'   uses both data and station governance. If `governance_id` contains a single
#'   value, then the matching is strict. However, if `governance_id` contains
#'   multiple values, then the matching is more complicated. For example, if
#'   `governance_id == c("OSPAR", "AMAP")`, then measurements with
#'   `is_ospar_monitoring == TRUE` and `"is_amap_monitoring == FALSE"` are
#'   matched to stations where `station_programgovernance` contains `"OSPAR";
#'   measurements with `is_ospar_monitoring == FALSE` and `is_amap_monitoring ==
#'   TRUE` are matched with stations where `station_programgovernance` contains
#'   `"AMAP"`; but measurements where `is_ospar_monitoring == TRUE` and
#'   `is_amap_monitoring == TRUE` are matched to stations where
#'   `station_programgovernance` contains either `"OSPAR"` or `"AMAP"`.
#' * governance_id: a vector of strings containing one or more of `"OSPAR"`,
#'   `"HELCOM"` and `"AMAP"`.
#' * grouping: a logical with `TRUE` indicating that stations will be grouped
#'   into meta-stations as specified by `station_asmtmimeparent` in the station
#'   dictionary. Defaults to `FALSE` apart from when `info$purpose == "OSPAR"`.
#' * check_coordinates: a logical with `TRUE` indicating that, when
#'   stations are matched by name, the sample coordinates must also be within
#'   the station geometry. No implemented yet, so defaults ot `FALSE`.
#'
#' @returns A data frame containing the contaminant data augmented by variables
#'   containing the station code and the station name
#' 
add_stations <- function(data, stations, info){

  # silence non-standard evaluation warnings
  .data <- .id <- .year <- .replaced <- .grouped <- .order <- NULL
  country <- station_geometry <- datatype <- NULL

  cat("\nMatching data with station dictionary\n", sep = "")

  
  control <- info$add_stations
  
  
  # get ordering variable, so output is in the original order
  
  data <- dplyr::mutate(data, .order = 1:nrow(data))
  
  
  # transform .year to the year in date; not the same as monitoring year 
  
  data$.year <- substr(data$date, 0, 4)
  data$.year <- as.numeric(data$.year)
  x <- dplyr::rename(data, station_name = "statn")
  stations <- dplyr::filter(stations, .data$station_deprecated == "FALSE")

  
  # restrict stations to convention area
  
  if (!is.null(control$area)) {
    cat(
      " - restricting to stations in these convention areas:", 
      paste(control$area, collapse = ", "), 
      "\n"
    )
    id <- match(control$area, c("OSPAR", "HELCOM", "AMAP"))
    id <- c("is_ospar_area", "is_helcom_area", "is_amap_area")[id]
    ok <- apply(stations[id], 1, any)
    stations <- stations[ok, ]
  } else {
    cat(" - no restriction of stations to a convention area\n")
  }
  

  # restrict stations by datatype
  # need further restrictions below to differentiate contaminants and effects
  
  if (control$datatype) {
    id <- switch(
      info$compartment, 
      biota = c("CF", "EF"), 
      sediment = c("CS", "ES"),
      water = c("CW", "EW")
    )
    cat(" - restricting to station with data types", id[1], "or", id[2], "\n")
    id <- paste(id, collapse = "|")
    stations <- dplyr::filter(stations, grepl(id, .data$station_datatype))
  } else {
    cat(" - no restriction of stations by data type\n")
  }
  

  # restrict stations by purpose of monitoring
  
  if (control$temporal) {
    cat(" - restricting to stations marked for temporal monitoring\n")
    stations <- dplyr::filter(stations, grepl("T", .data$station_purpm))
  } else {
    cat(" - no restriction of stations by purpose of monitoring\n")
  }
  
  
  # messages about other restrictions 

  if (control$governance_type %in% c("data", "both")) {
    cat(
      " - restricting to data with program governance:", 
      paste(control$governance_id, collapse = ", "), 
      "\n"
    )
    stations <- dplyr::filter(stations, grepl("T", .data$station_purpm))
  } else {
    cat(" - no restriction of data by program governance\n")
  }
  
  if (control$governance_type %in% c("stations", "both")) {
    cat(
      " - restricting to stations with program governance:", 
      paste(control$governance_id, collapse = ", "), 
      "\n"
    )
    stations <- dplyr::filter(stations, grepl("T", .data$station_purpm))
  } else {
    cat(" - no restriction of stations by program governance\n")
  }
  

  # ensure in the stations file only year is in the validity dates variables
  # make this numeric
  
  stations <- dplyr::mutate(
    stations, 
    station_activefromdate = substr(.data$station_activefromdate, 0, 4),
    station_activefromdate = as.numeric(.data$station_activefromdate),
    station_activeuntildate = substr(.data$station_activeuntildate, 0, 4),
    station_activeuntildate = as.numeric(.data$station_activeuntildate)
  )

  
  # deal with degenerate case where station_longituderange or 
  # station_latituderange equals zero (so polygon collapses to a line)

  ok <- stations$station_longituderange > 0 & stations$station_latituderange > 0
  if (!all(ok)) {

    stations <- split(stations, ok)

    stations[["FALSE"]] <- dplyr::mutate(
      stations[["FALSE"]], 
      station_longituderange = pmax(.data$station_longituderange, 0.000001),
      station_latituderange = pmax(.data$station_latituderange, 0.000001)
    )
     
    stations[["FALSE"]]$station_geometry <- mapply(
      FUN = function(longitude, latitude, longituderange, latituderange) {
        point <- c(longitude, latitude)
        
        range <- c(longituderange, latituderange)
        mult <- c(1, 1, -1, 1, -1, -1, 1, -1, 1, 1)
        
        polygon <- rep(point, 5) + rep(range, 5) * mult
        polygon <- matrix(polygon, ncol = 2, byrow = TRUE)
      
        geom <- sf::st_geometrycollection(
          c(sf::st_point(point), sf::st_polygon(list(outer = polygon)))
        )  
      
        sf::st_as_text(geom)
      },
      longitude = stations[["FALSE"]]$station_longitude,
      latitude = stations[["FALSE"]]$station_latitude,
      longituderange = stations[["FALSE"]]$station_longituderange,
      latituderange = stations[["FALSE"]]$station_latituderange,
      USE.NAMES = FALSE
    )
      
    stations <- unsplit(stations, ok)      

  }
  
  
  # create the datatype variable for matching (if required) 

  wk <- switch(info$compartment, biota = "F", sediment = "S", water = "W")

  x <- dplyr::mutate(
    x, 
    datatype = dplyr::case_when(
      startsWith(.data$pargroup, "O-")                 ~ paste0("C", wk), 
      startsWith(.data$pargroup, "OC-")                ~ paste0("C", wk), 
      startsWith(.data$pargroup, "OP-")                ~ paste0("C", wk), 
      .data$pargroup %in% c("I-MET", "I-RNC")          ~ paste0("C", wk),
      .data$pargroup %in% c("B-MBA", "B-TOX", "B-END") ~ paste0("E", wk), 
      .data$pargroup %in% c("B-GRS", "B-HST")          ~ paste0("D", wk), 
      TRUE                                             ~ "AUX"
    )
  )
  

  # split data into those records that are to be matched by name (x_name) and
  # those that are to be matched by coordinates
  
  if (control$method == "name") {
    
    x_name <- x
    x_coor <- x[FALSE, ]
  
  } else if (control$method == "coordinates") {
    
    x_name <- x[FALSE, ]
    x_coor <- x
  
  } else {
    
    # based on specification by OSPAR and HELCOM data assessors
    # .id is a logical giving records to match by name
    
    # these countries always match by name

    x <- dplyr::mutate(
      x, 
      .id = country %in% c(
        "Denmark", "Ireland", "Norway", "Portugal", "Sweden", "United Kingdom"
      ) 
    )
      
    # these countries match by name for some compartments / years
        
    x <- switch(
      info$compartment,
      biota = dplyr::mutate(
        x, 
        .id = .id | 
          (.data$country == "France") |
          (.data$country == "Spain" & .data$.year > 2004) |
          (.data$country == "The Netherlands" & .data$.year > 2006) |
          (.data$country == "Germany" & .data$.year > 2022)
      ),
      sediment = dplyr::mutate(
        x, 
        .id = .id |   
          (.data$country == "France" & .data$.year > 2008) |
          (.data$country == "Spain" & .data$.year > 2004) |
          (.data$country == "The Netherlands" & .data$.year > 2006) |
          (.data$country == "Germany" & .data$.year > 2022)
      ),
      water = dplyr::mutate(
        x, 
        .id = .id |  
          (.data$country == "France") |
          (.data$country == "Spain" & .data$.year > 2004) |
          (.data$country == "The Netherlands" & .data$.year > 2006) |
          (.data$country == "Germany" & .data$.year > 2022)
      ),
    )
    
    # and Germany always matches by name for HELCOM biota
    
    if (info$compartment == "biota" && info$purpose == "HELCOM") {
      x <- dplyr::mutate(x, .id = .id | (.data$country %in% "Germany"))
    }
    
    
    # split into two groups
      
    x_name <- dplyr::filter(x, .id)
    x_coor <- dplyr::filter(x, !.id)
    
    x_name <- dplyr::select(x_name, -.id)
    x_coor <- dplyr::select(x_coor, -.id)

  }
    
     

  ## Start of loop for coordinates matching
  
  cat(" - matching", nrow(x_coor), "records by coordinates\n")
  
  if (nrow(x_coor) > 0L) {
  
    id <- c(
      "country", "latitude",'longitude', '.year', 'datatype', 'is_ospar_monitoring', 
      'is_amap_monitoring', 'is_helcom_monitoring'
    ) 
    x_coor_unique <- unique(x_coor[id])
    
    res <- data.frame()
    # check <- data.frame()
    
    for(i in 1:nrow(x_coor_unique)) {
      
      file <- x_coor_unique[i,]
      
      # restrict stations to appropriate country and years 
      
      stations_subset <- dplyr::filter(
        stations, 
        .data$station_country == file$country, 
        .data$station_activefromdate <= file$.year, 
        is.na(.data$station_activeuntildate) | 
          .data$station_activeuntildate >= file$.year
      )
      
      # ensure match to contaminants or effects datatype 
      
      if (control$datatype && file$datatype != "AUX") {
        stations_subset <- dplyr::filter(
          stations_subset, 
          grepl(file$datatype, .data$station_datatype)
        )
      }
      
      # match by program governance
      # have already (partially) restricted by station governance
      
      if (control$governance_type %in% c("data", "both")) {
        
        id <- match(control$governance_id, c("OSPAR", "HELCOM", "AMAP"))
        id <- c("is_ospar_monitoring", "is_helcom_monitoring", "is_amap_monitoring")[id]
        
        # identify which conventions the measurement is marked for
        
        file_gov <- sapply(file[id], all)
        
        # if type == "data", only need to update stations subset if the measurement
        # does not satisfy the data governance requirement (might be easier to return)
        
        if (!any(file_gov)) {
          stations_subset <- stations_subset[FALSE, ]
        } else if (control$governance_type == "both") {
          station_gov <- control$governance_id[file_gov]
          station_gov <- paste(station_gov, collapse = "|")
          stations_subset <- dplyr::filter(
            stations_subset,
            grepl(station_gov, .data$station_programgovernance)
          )
        }
        
      } 
      
      sd <- sf::st_as_sf(
        stations_subset,
        wkt = "station_geometry",
        crs = sf::st_crs(4326)
        # 4326 is the EPSG code for the datum/projection used (https://epsg.io/4326). 
        # It needs to be specified so that the spatial functions know what 
        # coordinate system should be used when calculating distance/area etc. 
        # The 4326 is the WGS84 system used by most GPS systems        
      )
      
      dpoint <- sf::st_point(c(file$longitude,file$latitude))
      dpoint_sfc <- sf::st_sfc(dpoint)
      dpoint_sfc_4326 <- sf::st_set_crs(dpoint_sfc, 4326)
      
      sd1 <- sd[dpoint_sfc_4326, op = sf::st_intersects]
      
      if(nrow(sd1)> 1){
        
        part <- data.frame()
        for(i in 1:nrow(sd1)) {
          sd2 <- sd1[i,]
          sdpoint <- sf::st_point(c(sd2$station_longitude,sd2$station_latitude))
          # make a simple feature collection sfc of sdpoint
          sdpoint_sfc <- sf::st_sfc(sdpoint)
          # make a version of sdpoint_sfc using crs= 4326
          sdpoint_sfc_4326 <- sf::st_set_crs(sdpoint_sfc,4326)
          # calculate distance and turn to scalar
          sd_dist <- sf::st_distance(dpoint_sfc_4326, sdpoint_sfc_4326)
          attributes(sd_dist) <- NULL
          sd2$dist <- sd_dist
          part <- rbind(part, sd2)
        }
        part <- as.data.frame(part)
        
        part <- dplyr::select(part, - station_geometry)
        
        sd1 <- dplyr::filter(part, .data$dist == min(.data$dist))
        # if both distances are the same, keep the station with the highest code, 
        # should be the newest one
        sd1 <- dplyr::filter(sd1, .data$station_code == max(.data$station_code)) 
        result <- dplyr::mutate(file, station_code = sd1$station_code)
        res <- rbind(res, result)
        
      } else if(nrow(sd1)==1) {
        
        result <- dplyr::mutate(file, station_code = sd1$station_code)
        res <- rbind(res, result)
        
      } else {
        
        result <- dplyr::mutate(file,  station_code = NA)
        res <- rbind(res, result)
        
      }
      
    }
    
    
    # all unique 
    # res_unique<- unique(res)
    
    x_coor_matched <- dplyr::left_join(
      x_coor, 
      res, 
      by = c(
        "country", "latitude", "longitude", ".year", "datatype", 
        "is_amap_monitoring", "is_helcom_monitoring", "is_ospar_monitoring"
      )
    )
    # Remove created variables
    x_coor_matched <- dplyr::select(x_coor_matched, -datatype, -.year)
    
  }
        
  ## Start of loop for name-matching
  
  cat(" - matching", nrow(x_name), "records by station name\n")

  if (nrow(x_name) > 0L) {
  
    id <- c(
      "country", '.year', 'station_name', 'datatype', 'is_ospar_monitoring', 
      'is_amap_monitoring', 'is_helcom_monitoring'
    )
    x_name_unique <- unique(x_name[id])
    
    res <- data.frame()
    
    for(i in 1:nrow(x_name_unique)) {
      file <- x_name_unique[i,]
      
      # restrict stations to appropriate country and years 
      
      stations_subset <- dplyr::filter(
        stations, 
        .data$station_country == file$country, 
        .data$station_activefromdate <= file$.year, 
        is.na(.data$station_activeuntildate) | 
          .data$station_activeuntildate >= file$.year
      )
      
      # restrict stations to station_name match
      
      stations_subset <- dplyr::filter(
        stations_subset, 
        .data$station_name == file$station_name
      )
      
      # ensure match to contaminants or effects datatype 
      
      if (control$datatype && file$datatype != "AUX") {
        stations_subset <- dplyr::filter(
          stations_subset, 
          grepl(file$datatype, .data$station_datatype)
        )
      }
      
      # match by program governance
      # have already (partially) restricted by station governance
      
      if (control$governance_type %in% c("data", "both")) {
        
        id <- match(control$governance_id, c("OSPAR", "HELCOM", "AMAP"))
        id <- c("is_ospar_monitoring", "is_helcom_monitoring", "is_amap_monitoring")[id]
        
        # identify which conventions the measurement is marked for
        
        file_gov <- sapply(file[id], all)
        
        # if type == "data", only need to update stations subset if the measurement
        # does not satisfy the data governance requirement (might be easier to return)
        
        if (!any(file_gov)) {
          stations_subset <- stations_subset[FALSE, ]
        } else if (control$governance_type == "both") {
          station_gov <- control$governance_id[file_gov]
          station_gov <- paste(station_gov, collapse = "|")
          stations_subset <- dplyr::filter(
            stations_subset,
            grepl(station_gov, .data$station_programgovernance)
          )
        }
      }
      
      if(nrow(stations_subset) == 0L){
        result <- dplyr::mutate(file, station_code = NA)
        res <- rbind(res, result)
      }
      
      if(nrow(stations_subset) == 1L){
        result <- dplyr::mutate(file, station_code = stations_subset$station_code)
        res <- rbind(res, result)
      }
      
      if (nrow(stations_subset) > 1L) {
        stop(
          "Multiple station matches when joining by name - contact HARSAT development team", 
          call. = FALSE
        )
      }
    }
    
    # res_unique<- unique(res)
    
    x_name_matched <- dplyr::left_join(
      x_name, 
      res, 
      by = c(
        "country", "station_name", ".year", "datatype", "is_amap_monitoring", 
        "is_helcom_monitoring", "is_ospar_monitoring"
      )
    )
    
    x_name_matched <- dplyr::select(x_name_matched, -datatype, -.year)
    
  }
    
  
  # full_match equals data file plus a station_code variable
  
  if (nrow(x_coor) == 0L) {
    full_match <- x_name_matched
  } else if (nrow(x_name) == 0L) {
    full_match <- x_coor_matched
  } else {
    full_match <- rbind(x_coor_matched, x_name_matched)
  }
  
  full_match <- dplyr::rename(full_match, statn = "station_name")
  
  
  # get station_name associated with station_code
  # rename as sd_code_match and sd_name_match
  
  full_match <- dplyr::left_join(
    full_match, 
    stations[c("station_code", "station_name")],
    by = "station_code"
  )
  
  full_match <- dplyr::rename(
    full_match,
    sd_code_match = "station_code",
    sd_name_match = "station_name"
  )
  
  
  # get 'current' station choice 
  
  full_match$sd_code_current <- full_match$sd_code_match
  full_match$sd_name_current <- full_match$sd_name_match
  
  
  # update current stations with any replacements
  
  full_match <- dplyr::left_join(
    full_match,
    stations[c("station_code", "station_replacedby")],
    by = c("sd_code_current" = "station_code")
  )  
  
  # get station_name associated with station_replacedby
  # important: this will be missing if station_replacedby is no longer in
  # the station dictionary because e.g. it has the wrong data type
  
  full_match <- dplyr::left_join(
    full_match,
    stations[c("station_code", "station_name")],
    by = c("station_replacedby" = "station_code")
  )  
  
  # rename as station_code_replaced and station_name_replaced
  
  full_match <- dplyr::rename(
    full_match,
    sd_code_replaced = "station_replacedby",
    sd_name_replaced = "station_name"
  )
  
  # remove replaced codes that are no longer in the station dictionary
  # update 'current' station code
  
  full_match <- dplyr::mutate(
    full_match,
    .replaced = !is.na(.data$sd_name_replaced),
    sd_code_replaced = ifelse(
      .replaced, 
      .data$sd_code_replaced,
      NA_character_
    ),
    sd_code_current = ifelse(
      .replaced, 
      .data$sd_code_replaced, 
      .data$sd_code_current
    ),
    sd_name_current = ifelse(
      .replaced, 
      .data$sd_name_replaced, 
      .data$sd_name_current
    ),
    .replaced = NULL
  )
  
  # check no double replacements required
  
  check <- dplyr::left_join(
    full_match[c("sd_code_replaced")], 
    stations[c("station_code", "station_replacedby")],
    by = c("sd_code_replaced" = "station_code")
  )
  
  if (!all(is.na(check$station_replacedby))) {
    stop(
      "two levels of station replacement required - contact HARSAT development team", 
      call. = FALSE
    )
  }

  
  if (control$group) {

    cat(" - grouping stations using station_asmtmimegovernance\n")
    
    full_match <- dplyr::left_join(
      full_match,
      stations[c("station_code", "station_asmtmimeparent")],
      by = c("sd_code_current" = "station_code")
    )  
    
    # get station_name associated with station_asmtmimeparent
    # important: this will be missing if station_asmtmimeparent is no longer in
    # the station dictionary because e.g. it has the wrong data type

    full_match <- dplyr::left_join(
      full_match,
      stations[c("station_code", "station_name")],
      by = c("station_asmtmimeparent" = "station_code")
    )  
    
    # rename as station_code_grouped and station_name_grouped

    full_match <- dplyr::rename(
      full_match,
      sd_code_grouped = "station_asmtmimeparent",
      sd_name_grouped = "station_name"
    )
    
    # remove grouped codes that are no longer in the station dictionary
    # update 'current' station code

    full_match <- dplyr::mutate(
      full_match,
      .grouped = !is.na(.data$sd_name_grouped),
      sd_code_grouped = ifelse(
        .grouped, 
        .data$sd_code_grouped,
        NA_character_
      ),
      sd_code_current = ifelse(
        .grouped, 
        .data$sd_code_grouped, 
        .data$sd_code_current
      ),
      sd_name_current = ifelse(
        .grouped, 
        .data$sd_name_grouped, 
        .data$sd_name_current
      ),
      .grouped = NULL
    )
    
    # check no double groupings required
    
    check <- dplyr::left_join(
      full_match[c("sd_code_grouped")], 
      stations[c("station_code", "station_asmtmimeparent")],
      by = c("sd_code_grouped" = "station_code")
    )
    
    if (!all(is.na(check$station_asmtmimeparent))) {
      stop(
        "two levels of station grouping required - contact HARSAT development team", 
        call. = FALSE
      )
    }

  }
  

  # rename current stations as station_code and station_name
  
  full_match <- dplyr::relocate(
    full_match,
    station_code = "sd_code_current",
    station_name = "sd_name_current",
    .after = dplyr::last_col()
  )


  cat(
    " -", sum(!is.na(full_match$station_code)), "of", nrow(full_match), 
    "records have been matched to a station\n"
  )
  
  
  # reorder
  
  full_match <- dplyr::arrange(full_match, .order)
  full_match <- dplyr::select(full_match, -.order)

  full_match
  
}



finalise_data <- function(data, info) {
  
  # silence non-standard evaluation warnings
  .data <- NULL
  retain <- species <- NULL

  # common to all compartments

  data <- dplyr::rename(
    data, 
    station_submitted = "statn", 
    year = "myear",
    sample_latitude = "latitude",
    sample_longitude = "longitude",
    depth_upper = "dephu", 
    depth_lower = "dephl",
    determinand = "param", 
    matrix = "matrx", 
    unit = "munit",
    limit_detection = "detli",
    limit_quantification = "lmqnt",
    uncertainty = "uncrt",
    unit_uncertainty = "metcu", 
    censoring = "qflag",
    method_analysis = "metoa",
    method_extraction = "metcx",
    method_pretreatment = "metpt",
    replicate = "tblparamid", 
    upload = "tbluploadid"
  )

  
  # compartment specific
  
  if (info$compartment %in% c("sediment", "water")) {
    
    data <- dplyr::rename(data, sample = "tblsampleid")
    
  } else if (info$compartment == "biota") {
    
    data <- dplyr::rename(
      data,
      sample_cluster = "tblsampleid", 
      sample = "tblbioid", 
      sex = "sexco",
      n_individual = "noinp"
    )    
  
  }
  
  
  # apply ICES filters
  
  # exclude suspect values in all compartments 
  
  data <- dplyr::mutate(
    data, 
    retain = !grepl("S", .data$vflag)
  )

  
  # biota:
  # create species column; use worms_accepted_name if it exists, and speci_name
  #   otherwise
  # exclude records with missing species
  
  if (info$compartment == "biota") {
    
    data = dplyr::mutate(
      data,
      species = .data$worms_accepted_name,
      species = ifelse(is.na(.data$species), .data$speci_name, .data$species),
      retain = retain & !is.na(species)
    )

    data <- dplyr::relocate(data, "species", .after = "worms_accepted_name")

  }
    

  # sediment:
  # exclude samples where the upper depth > 0
  # exclude samples where the matrix is not one of:
  #   SED20, SED62, SED63, SED90, SED100, SED500, SED1000, SED2000, SEDTOT
  
  if (info$compartment == "sediment") {
    
    ok_matrix <- paste0(
      "SED", 
      c("20", "62", "63", "90", "100", "500", "1000", "2000", "TOT")
    )

    data = dplyr::mutate(
      data, 
      retain = .data$retain & !is.na(.data$depth_upper) & 
        dplyr::near(.data$depth_upper, 0),
      retain = .data$retain & .data$matrix %in% ok_matrix, 
    )

  }

    
  # water: 
  # exclude samples where the upper depth > 5.5m
  # exclude samples where matrix == "SPM"
  # exclude samples where there is no filtration information
  # create filtration column
  
  if (info$compartment == "water") {

    data = dplyr::mutate(
      data, 
      retain = .data$retain & !is.na(.data$depth_upper) & 
        .data$depth_upper <= 5.5,
      retain = .data$retain & .data$matrix == "WT", 
      retain = .data$retain & !is.na(.data$method_pretreatment)
    )
    
    wk <- strsplit(data$method_pretreatment, "~|-")
    data$filtration <- sapply(
      wk, 
      function(x) {
        if (length(x) == 1L && is.na(x)) {
          NA_character_
        } else if (any(c("NF", "NONE") %in% x)) {
          "unfiltered" 
        } else {
          "filtered"
        }
      }
    )

    data <- dplyr::relocate(data, "filtration", .after = "matrix")
        
  }

  data
}


finalise_stations <- function(stations, info) {

  # rename columns that will be passed down into later funcions
  
  stations <- dplyr::rename(
    stations,
    country = "station_country",
    station_type = "station_mstat",
    waterbody_type = "station_wltyp"
  )

  stations
}
  



ctsm_read_QA <- function(file, path, purpose) {
  
  # silence non-standard evaluation warnings
  crmco <- crmmb <- crmmv <- NULL
  myear <- munit <- metoa <- metcx <- NULL
  dtype <- param <- tblanalysisid <- .data <- NULL

  # import functions
  # read in method data (and additional crm information)

  # read in data
  
  infile <- file.path(path, file)
  cat("\nReading QA data from '", infile, "'\n", sep = "")
  
  report_file_digest(infile)
  crm <- safe_read_file(
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


# tidy data ----

#' Tidy the data
#' 
#' Reduces the size of the extraction by removing redundant variables.
#' Any ad-hoc changes will usually be made between [`read_data`] and [`tidy_data`].
#' The output is in the correct format for [`create_timeseries`].
#'
#' @param ctsm_obj the time series data, typically returned from [`read_data`].
#' @export
tidy_data <- function(ctsm_obj) {
  
  # import_functions.R
    
  info <- ctsm_obj$info
  data <- ctsm_obj$data
  stations <- ctsm_obj$stations

  # set up oddity directory and back up any previous oddity files
  
  initialise_oddities(info$oddity_dir, info$compartment)
  
  
  # tidy station dictionary and contaminant data

  
  # drop data that are suspect or do not meet requirements 

  if (info$data_format %in% "ICES") {
    
    ok <- data$retain
    if (!all(ok)) {
      cat(
        "\nDropping", 
        sum(!ok), 
        "records from data flagged for deletion. Possible reasons are:\n", 
        "- vflag = suspect\n", 
        switch(
          info$compartment, 
          biota = "- species missing\n",
          sediment = paste(
            "- upper depth > 0m\n", 
            "- unusual matrix\n"
          ),
          water = paste(
            "- upper depth >5.5m\n", 
            "- filtration missing\n",
            "- matrix = 'SPM'\n"
          )
        )
      )
      data <- data[ok, ]
    }

  }    
    

  # drop data with no station_code and stations that are not in the data
  
  if (info$data_format %in% c("ICES", "external")) {
    
    ok <- !is.na(data$station_code)
    if (!all(ok)) {
      cat(
        "\nDropping", 
        sum(!ok), 
        "records from data that have no valid station code\n" 
      )
      data <- data[ok, ]
    }
    
    ok <- stations$station_code %in% data$station_code
    if (!all(ok)) {
      cat(
        "\nDropping", 
        sum(!ok), 
        "stations that are not associated with any data\n" 
      )
      stations <- stations[ok, ]
    }
    
  }
  
  
  stations <- tidy_stations(stations, info)

  data <- tidy_contaminants(data, info)


  ctsm_obj$stations <- stations
  ctsm_obj$data <- data
  
  ctsm_obj
}


tidy_stations <- function(stations, info) {

  # silence non-standard evaluation warnings
  .data <- NULL
  
  cat("\nCleaning station dictionary\n")
  
  # ensure country is consistent
  
  stations <- dplyr::mutate(
    stations, 
    country = stringr::str_to_title(.data$country)
  )
  
  
 # replace backward slash with _ in station (long) name

  stations <- dplyr::mutate(
    stations, 
    station_longname = gsub("\\", "_", .data$station_longname, fixed = TRUE)
  )
    


  # select useful columns
  
  
  col_id <- c(
    info$region$id, "country", "station_name", "station_code", "station_longname", 
    "station_latitude", "station_longitude", "station_type", "waterbody_type"
  )
  

  stations <- stations[col_id]
  
  
  # order data
  
  col_id <- c(info$region$id, "country", "station_name")
  
  stations <- dplyr::arrange(stations, dplyr::across(dplyr::all_of(col_id))) 
  
  stations
}






tidy_contaminants <- function(data, info) {
  
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
      "option in read_data.",
      call. = FALSE
    )
    
  }
  


  # add required variables 'replicate' and 'pargroup' (not present in external data)
  
  # have removed 'pargroup' for now as it can result in a lot of unrecognised 
  # determinands - 'pargroup' is added later in create_timeseries
  
  if (info$data_format == "external") {
    data$replicate <- seq(from = 1, to = nrow(data), by = 1)
    # data$pargroup <- ctsm_get_info(info$determinand, data$determinand, "pargroup")
  }
  
      
  data$sample <- as.character(data$sample)

  
  # retain useful variables
  
  var_id <- c(
    "station_code", "sample_latitude", "sample_longitude", 
    "year", "date", "time", "depth", 
    "species", "sex", "n_individual", "subseries", "sample", "replicate", 
    "determinand", "pargroup", "matrix", "basis", "filtration", 
    "method_analysis", "method_extraction", "method_pretreatment",
    "unit", "value", "censoring", "limit_detection", "limit_quantification", 
    "uncertainty", "unit_uncertainty", "alabo", "qalink"
  )
  
  data <- data[intersect(var_id, names(data))]
  
  
  data
}





#' Create a time series
#' 
#' Cleans the data and turns it into time series structures ready for assessment
#' 
#' @param ctsm.obj the CTSM object, as returned from `tidy_data`
#' @param determinands the determinands to use, by default derived by
#'   calling `ctsm_get_determinands`, which takes values from 
#'   the determinand reference table
#' @param determinands.control determinands control values, when needed
#' @param oddity_path a directory to write the oddities to
#' @param return_early a boolean (default `FALSE`), if `TRUE`, returns early 
#' @param print_code_warnings a boolean (default `FALSE`)
#' @param get_basis a basis function, by default [`get_basis_default`]
#' @param normalise boolean or function, if TRUE, uses a default normalization
#' @param normalise.control a list of control data for the normalization function
#' 
#' @return a completed timeseries object, which can be used for assessments
#' 
#' @export
create_timeseries <- function(
  ctsm.obj, 
  determinands = ctsm_get_determinands(ctsm.obj$info), 
  determinands.control = NULL, 
  oddity_path = "oddities", 
  return_early = FALSE, 
  print_code_warnings = FALSE, 
  get_basis = get_basis_default,
  normalise = FALSE, 
  normalise.control = list()) {

  # silence non-standard evaluation warnings
  .data <- .month <- .not_ok <- group <- value <- NULL
  determinand <- uncertainty <- uncertainty_sd <- uncertainty_rel <- species_group <- NULL

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

  oddity_path <- initialise_oddities(info$oddity_dir, info$compartment)


  # checks station dictionary:
  # - no backward slashes in station_name or station_longname
  # - no duplicated or missing station_code 
  
  ctsm_check_stations(station_dictionary)
  
  
  # determinand validation
  
  # retains determinands of interest, including auxiliary determinands and those
  # required by determinands.control$variables 
  # checks all determinands of interest are recognised by info$determinand

  wk <- ctsm_check_determinands(info, data, determinands, determinands.control)
  
  data <- wk$data
  determinands <- wk$determinands
  determinands.control <- wk$control

  
  cat("\nCleaning data\n")
  
  id <- c(
    "station_code", "date", "filtration", "species", "sample", "matrix", 
    "determinand"
  )
  id <- intersect(id, names(data))
  ord <- do.call("order", data[id])
  data <- data[ord, ]

  # merge with country and station_name from station dictionary to give more 
  # meaningful output in oddity files

  var_id <- c("country", "station_code", "station_name")
  var_id <- intersect(var_id, names(station_dictionary))

  data <- dplyr::left_join(data, station_dictionary[var_id], by = "station_code")
  data <- dplyr::relocate(data, dplyr::all_of(var_id))


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
  

  # species validation
  
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
  
  if (info$region$all) {
    
    data <- dplyr::left_join(
      data, 
      station_dictionary[c("station_code", info$region$id)],
      by = "station_code"
    )
    
    if (any(is.na(data[info$region.id]))) {
      cat("   Dropping stations with missing region information in station dictionary\n")
      data <- tidyr::drop_na(data, dplyr::all_of(info$region$id))
    }
    
    data <- data[setdiff(names(data), info$region$id)]
  }


  # add variables that are going to be useful throughout
  # pargroup is already in an ICES extraction, but might need to supply it
  # for external data (from info$determinand)
  
  data$group <- ctsm_get_info(
    info$determinand, data$determinand, "group", info$compartment, sep = "_"
  )

  if (info$compartment == "biota") {
    data$species_group <- ctsm_get_info(info$species, data$species, "species_group")
  }
  
  if (!"pargroup" %in% names(data)) {
    data$pargroup <- ctsm_get_info(info$determinand, data$determinand, "pargroup")
  }
  
  # NB distribution will be missing for auxiliary data

  data$distribution <- ctsm_get_info(
    info$determinand, data$determinand, "distribution", na_action = "output_ok"
  )
  
  
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

  data <- check_censoring(data, info, print_code_warnings)
  

  # ensure uncertainties are plausible

  data <- check_uncertainty(data, info, type = "reported")

  
  # convert all uncertainties to unit SD

  data <- dplyr::mutate(
    data, 
    uncertainty = dplyr::case_when(
      .data$unit_uncertainty %in% "U2" ~ .data$uncertainty / 2, 
      .data$unit_uncertainty %in% "%" ~ .data$value * .data$uncertainty / 100, 
      .default = .data$uncertainty
    ), 
    unit_uncertainty = "SD" 
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
    
    args = list(data = data, info = info, keep = i, drop = wk$det)
    if ("weights" %in% names(wk)) {
      args = c(args, list(weights = wk$weights))
    }
    
    data <- do.call(linkFunction, args)
  }  

  # drop any remaining unwanted determinands (from sum and perhaps bespoke functions);
  # could make this more elegant!
  
  id <- c(determinands, ctsm_get_auxiliary(determinands, info))
  
  data <- dplyr::filter(data, .data$determinand %in% id)


  # remove 'extra' data from time series which have a subseries classification 
  # for some but not all records

  data <- check_subseries(data, info)
    

  # set rownames to NULL(ie. back to auto numeric)
  
  rownames(data) <- NULL

  cat("\nCreating time series data\n")  


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
    data <- ctsm.imposex.check.femalepop(data, info)
  }
  

  # convert data to basis of assessment
  
  data <- ctsm_convert_to_target_basis(data, info, get_basis)
  

  if (return_early) {
    out <- c(
      out, 
      output_timeseries(data, station_dictionary, info, extra = "alabo")
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

  # check that all normal and lognormal data have uncertainties
  
  data <- ctsm_check(
    data, 
    distribution %in% c("normal", "lognormal") & !is.na(concentration) & 
      is.na(uncertainty), 
    action = "delete", 
    message = "Missing uncertainties which cannot be imputed", 
    file_name = "missing_uncertainties", 
    info = info
  )
  

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
      biota = stop("there is no default normalisation function for biota"),
      ## TODO: implement normalise_biota
      # biota = 
      #   normalise_biota(data, station_dictionary, info, normalise.control),
      sediment = stop("there is no default normalisation function for sediment"),
      ## TODO: implement normalise_sediment
      # sediment = 
      #   normalise_sediment(data, station_dictionary, info, normalise.control),
      ## TODO: implement normalise_water
      water = stop("there is no default normalisation function for water")
    )
  } else if (is.function(normalise)) {
    data <- normalise(data, station_dictionary, info, normalise.control)
  }
    
  # check whether implausible uncertainties have been calculated during the 
  #   data processing (e.g. during normalisation)
  # if so - make concentration, uncertainty and censoring missing
  
  data <- check_uncertainty(data, info, type = "calculated")
  
  
  # final check to ensure all normal and lognormal data have an uncertainty
  
  notok <- data$distribution %in% c("normal", "lognormal") & 
    !is.na(data$concentration) & is.na(data$uncertainty)
  
  if (any(notok)) {
    stop(
      "uncertainties missing where they should be present: \n", 
      "contact HARSAT development team")
  }
  
  
  # drop groups of data at stations with no data in recent years

  cat("   Dropping groups of compounds / stations with no data between", 
      min(info$recent_years), "and", max(info$recent_years), "\n")
  id_names <- intersect(c("station_code", "filtration", "species", "group"), names(data))
  id <- do.call("paste", data[id_names])
  ok <- id %in% id[is.recent(data$year) & !is.na(data$concentration)]
  data <- data[ok, ]
  

  # create seriesID and timeSeries object
  
  out <- c(
    out, 
    output_timeseries(data, station_dictionary, info)
  )
  
  out
}



#' Extracts timeSeries
#'
#' Gets the timeSeries component of a `harsat` object, optionally having added
#' extra information about each station
#'
#' @param harsat_obj A `harsat` object following a call to [`create_timeseries`]`.
#' @param add logical (default `TRUE`), if `TRUE`, adds extra information 
#'   about each station; if `FALSE`, simply returns the existing timeseries,
#'
#' @return A data.frame containing the timeSeries component with (optionally)
#'   extra information about each station.
#' @export
get_timeseries <- function(harsat_obj, add = TRUE) {
  
  if (!all(c("timeSeries", "stations") %in% names(harsat_obj))) {
    stop("the time series have not yet been created")
  }
  
  if (!add) {
    return(harsat_obj$timeSeries)
  }
  
  out <- dplyr::left_join(
    harsat_obj$timeSeries, 
    harsat_obj$stations[c("station_code", "station_name", "country")], 
    by = "station_code") 
  
  out <- dplyr::relocate(
    out, 
    c("station_name", "country"), 
    .after = "station_code"
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
  

output_timeseries <- function(data, station_dictionary, info, extra = NULL) {
  
  # silence non-standard evaluation warnings
  .data <- .group <- seriesID <- NULL

  # import_functions.R
  
  # identifies timeseries in the data and creates timeSeries metadata object


  # order data 

  id = c(
    "station_code", "species", "filtration", "year", "sex", "sample", "group", 
    "determinand"
  )
  
  data <- dplyr::arrange(data, dplyr::across(dplyr::any_of(id)))


  # select variables of interest

  id <- c(
    "station_code", "sample_latitude", "sample_longitude",  
    "species", "sex", "depth",
    "year", "date", "time", "sample",   
    "matrix", "filtration", "subseries", "group", "determinand", "basis", 
    "unit", "value", 
    "method_analysis", "n_individual", 
    "concOriginal", "censoringOriginal", "uncrtOriginal", 
    "concentration", "new.basis", "new.unit", "censoring",  
    "limit_detection", "limit_quantification", "uncertainty"
  )
  
  if (!is.null(extra)) {
    id <- c(id, extra)
  }
  
  auxiliary <- ctsm_get_auxiliary(data$determinand, info)
  auxiliary_id <- paste0(
    rep(auxiliary, each = 5), 
    c("", ".censoring", ".limit_detection", ".limit_quantification", ".uncertainty") 
  )
    
  id <- c(id, auxiliary_id)

  data <- dplyr::select(data, dplyr::any_of(id))

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
    id <- c(id, "filtration", "subseries") 
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
    dplyr::all_of(names(timeSeries)), 
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
  
  timeSeries <- tibble::column_to_rownames(timeSeries, "seriesID")

  
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


initialise_oddities <- function(path, compartment) {

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



determinand.link.check <- function(data, info, keep, drop, printDuplicates = TRUE, ...) {

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
      info = info,
      ...
    )
  }
  
  data[!(dups & dropID), ]
}  
  

determinand.link.replace <- function(data, info, keep, drop, ...) {

  # core function for relabelling determinand 'drop' as determinand 'keep'
  # most of the work is checking that there aren't data submitted as both for the same
  # sample
  
  stopifnot(length(keep) == 1, length(drop) == 1)
  
  if (any(data$determinand %in% drop)) 
    cat("   Data submitted as", drop, "relabelled as", keep, "\n")
  else return(data)
  
  
  # check for samples with both drop and keep and, if they exist, delete drop

  data <- determinand.link.check(data, info, keep, drop, ...)
  
  
  # relabel the levels so that drop becomes keep
  
  id <- data$determinand %in% drop
  data$determinand[id] <- keep
  
  data
}  


determinand.link.imposex <- function(data, info, keep, drop, ...) {
  
  stopifnot(length(keep) == 1, length(drop) == 1)

  detID <- c(keep, drop)
  
  # for imposex, indices and stages aren't linked by sample, but by visit
  # will assume, for simplicity, that only one visit per year
  
  visitID <- with(data, paste(station_code, year))
  
  
  # find visits when both individuals and stages reported and check consistent
  
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
    info = info,
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

determinand.link.BBKF <- function(data, info, keep, drop, ...) {
  
  stopifnot(
    identical(keep, "BBKF"), 
    identical(sort(drop), c("BBF", "BBJF", "BBJKF", "BKF"))
  )
  
  # first sum samples with both BBF and BKF
  
  data <- determinand.link.sum(data, info, "BBKF", c("BBF", "BKF"))
  
  # now sum samples with both BBJF and BKF to give BBJKF
  
  data <- determinand.link.sum(data, info, "BBJKF", c("BBJF", "BKF"))
  
  # now replace BBJKF with BBKF
  
  data <- determinand.link.replace(data, info, "BBKF", "BBJKF")
  
  data
}



assign("determinand.link.LIPIDWT%", function(data, info, keep, drop, ...) {

  stopifnot(identical(keep, "LIPIDWT%"), identical(sort(drop), c("EXLIP%", "FATWT%")))

  # if multiple values present, choose FATWT%, then LIPIDWT%, then EXLIP% (from Foppe)
  
  data <- determinand.link.check(
    data, info, keep = "LIPIDWT%", drop = "EXLIP%", printDuplicates = FALSE, ...
  )
  data <- determinand.link.check(
    data, info, keep = "FATWT%", drop = "EXLIP%", printDuplicates = FALSE, ...
  )
  data <- determinand.link.check(
    data, info, keep = "FATWT%", drop = "LIPIDWT%", printDuplicates = FALSE, ...
  )

  if (!any(data$determinand %in% drop)) return(data)

  # relabel the levels so that drop becomes keep
  
  cat("   Data submitted as EXLIP% or FATWT% relabelled as LIPIDWT%", "\n")

  id <- data$determinand %in% drop
  data$determinand[id] <- keep
  
  data
})  


determinand.link.sum <- function(data, info, keep, drop, ...) {
  
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

  if (sum(sum_ID) == 0L)
    return(data)
  
  
  dropTxt <- paste(drop, collapse = ", ")
  cat("   Data submitted as", dropTxt, "summed to give", keep, fill = TRUE)


  # get relevant sample matrix combinations
  
  data <- split(data, with(data, determinand %in% drop & sum_ID))
  
  ID <- with(data[["TRUE"]], paste(sample, matrix))

  summed_data <- by(data[["TRUE"]], ID, function(x) {
    
    # check all bases are the same 

    stopifnot(dplyr::n_distinct(x$basis) == 1)
    
    if (!all(drop %in% x$determinand)) return(NULL)

    
    # adjust values if units vary
    # ideally use unit in info$determinand, but makes it more awkward because
    # have to pass in compartment
    
    if (dplyr::n_distinct(x$unit) > 1) {

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
    else if (dplyr::n_distinct(x$censoring) == 1) 
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




determinand.link.TEQDFP <- function(data, info, keep, drop, weights) {

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
  
  if (sum(sum_ID) == 0)
    return(data)
  
  dropTxt <- paste(drop, collapse = ", ")
  cat("   Data submitted as", dropTxt, "summed to give", keep, fill = TRUE)
  
  
  # get relevant sample matrix combinations
  
  data <- split(data, with(data, determinand %in% drop))
  
  ID <- with(data[["TRUE"]], paste(sample, matrix))
  
  summed_data <- by(data[["TRUE"]], ID, function(x) {
    
    # check all bases are the same 
    
    if (!all(drop %in% x$determinand)) return(NULL)

    stopifnot(dplyr::n_distinct(x$basis) == 1)
    
    
    # convert to ug/kg and then to TEQ
    
    id <- c("value", "uncertainty", "limit_detection", "limit_quantification")
    
    x[id] <- lapply(x[id], convert_units, from = x$unit, to = "ug/kg")
    
    TEQ <- weights[as.character(x$determinand)]
    
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
    else if (dplyr::n_distinct(x$censoring) == 1) 
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
  nSummed <- if (is.null(summed_data)) 0 else nrow(summed_data)
  nLost <- nTotal - nSummed
  
  if (nLost > 0) 
    message("     ", nLost, " of ", nTotal, " samples lost due to incomplete submissions")

  
  # combine data for both drop and keep and then add back into main data set
  
  data[["TRUE"]] <- rbind(data[["TRUE"]], summed_data)

  data <- do.call(rbind, data)
  
  data
}  


check_censoring <- function(data, info, print_code_warnings) {
  
  # silence non-standard evaluation warnings
  value <- limit_detection <- limit_quantification <- NULL

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


check_uncertainty <- function(data, info, type = c("reported", "calculated")) { 

  # import_functions.r
  
  # uncertainties must be non-negative for all data
  # uncertainties must be strictly positive for normal or lognormal data
  # relative uncertainties must be within specified range (1, 100) default for
  #   lognormal data
  
  # type = reported is used for submitted data
  # type = calculated is used to check whether implausible uncertainties have
  #   been created in e.g. the normalisation process
  
  type <- match.arg(type)
  
  
  # calculate relative uncertainties for lognormal data
  # use value for reported data and concentration for calculated data
  
  id <- switch(type, reported = "value", calculated = "concentration")

  data <- dplyr::mutate(
    data, 
    .ok = .data$distribution %in% "lognormal", 
    relative_uncertainty = dplyr::case_when(
      .ok & .data$unit_uncertainty %in% "SD" ~ 
        100 * .data$uncertainty / .data[[id]], 
      .ok & .data$unit_uncertainty %in% "U2" ~ 
        100 * .data$uncertainty / (2 * .data[[id]]),
      .ok & .data$unit_uncertainty %in% "%" ~ .data$uncertainty,
      .default = NA_real_
    ), 
    .ok = NULL
  )
    
  data <- dplyr::mutate(
    data,
    reason = dplyr::case_when(
      .data$uncertainty < 0                                        ~ "negative", 
      .data$distribution %in% c("normal", "lognormal") & 
        .data$uncertainty == 0                                     ~ "zero",
      .data$distribution %in% "lognormal" & 
        .data$relative_uncertainty <= info$relative_uncertainty[1] ~ "small",
      .data$distribution %in% "lognormal" & 
        .data$relative_uncertainty >= info$relative_uncertainty[2] ~ "large", 
      .default = "ok"
    )
  )    

  data <- dplyr::relocate(
    data, 
    "relative_uncertainty", 
    .after = "unit_uncertainty"
  )
  
  data <- dplyr::relocate(data, "reason")

  if (type == "reported") {
    message <- "Implausible uncertainties reported with data"
    file_name <- "implausible_uncertainties_reported"
    missing_id <- "uncertainty"
  } 
  
  if (type == "calculated") {
    message <- "Implausible uncertainties calculated in data processing"
    file_name <- "implausible_uncertainties_calculated"
    missing_id <- c("concentration", "uncertainty", "censoring")
  }
      
  data <- ctsm_check(
    data, 
    reason != "ok", 
    action = "make.NA", 
    message = message, 
    file_name = file_name, 
    missing_id = missing_id,
    info = info
  )
  
  data$reason <- NULL
  data$relative_uncertainty <- NULL

  data
}
  
  
  
check_subseries <- function(data, info) {

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
    dplyr::any_of(c("station_code", "species", "determinand", "matrix")), 
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


#' Normalises sediment concentrations, OSPAR vwersion
#' 
#' @param data the data object
#' @param station_dictionary the station dictionary
#' @param info the information object
#' @param control control values
#' @export
normalise_sediment_OSPAR <- function(data, station_dictionary, info, control) {
  
  # silence non-standard evaluation warnings
  .data <- NULL

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
        
        # if ("PNL" %in% Cdigestion) {
        #   warning(
        #     "ad-hoc fix to deal with new NL method - must resolve for next assessment", 
        #     call. = FALSE
        #   )
        # }
        
        
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

#' Normalises sediment concentrations, HELCOM version
#' 
#' @param data the data object
#' @param station_dictionary the station dictionary
#' @param info the information object
#' @param control control values
#' @export
normalise_sediment_HELCOM <- function(data, station_dictionary, info, control) {

  # silence non-standard evaluations 
  .data <- NULL
  CORG <- CORG.censoring <- CORG.uncertainty <- NULL
  LOIGN <- LOIGN.censoring <- LOIGN.uncertainty <- NULL
  .tmp <- .tmp.censoring <- .tmp.uncertainty <- NULL
  
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

#' Normalises biota concentrations, HELCOM vwersion
#' 
#' @param data the data object
#' @param station_dictionary the station dictionary
#' @param info the information object
#' @param control control values
#' @export
normalise_biota_HELCOM <- function(data, station_dictionary, info, control) {
  
  # silence non-standard evaluation warnings
  .data <- NULL

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

  # silence non-standard evaluation warnings
  .data <- NULL

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
#' 
#' @param data the data
#' @param subset optional, a subset expression
#' @param action one of `"relabel"`, `"convert"`, `"change_unit"` -- 
#'   note that `change_unit` moves units from tin units to conventional units
#' @param from either `"tin"` or `"cation"`
#' @param convert_var one of `"value"`, `"limit_detection"`, 
#'   `"limit_quantification"`, `"uncertainty"`
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


   
