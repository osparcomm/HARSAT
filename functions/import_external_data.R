# function to read in external data
# based on AMAP mercury data

add_non_ICES_data <- function(
  ICES_data, AMAP_data, AMAP_stations, 
  keep = c("all", "AMAP"), path = ".") {

  library("dplyr")
  library("readxl")
  library("tidyr")
  
  keep <- match.arg(keep)
  

  # utily function
  
  my_paste <- function(data) do.call(paste, data)
  

  # read in AMAP data
  
  file_id <- file.path(path, AMAP_data)
  
  id <- read_excel(file_id, n_max = 0)
  id <- names(id)


  # check all required columns are present
  
  required <- c(
    "country", "station_code", "station_name", "species", "sex", "determinand", "matrix", "basis", 
    "unit", "uncertainty", "methodUncertainty", "qflag", "sub.sample", "year", "value", 
    "limit_detection", "limit_quantification", "AMAP_group") 
 
  ok <- required %in% id
  
  if (!all(ok)) 
    stop("following columns missing from AMAP data file: ", paste(required[!ok], collapse = ", "))
    

  col_types <- rep("skip", length(id))
  col_types[id %in% c(
    "country", "station_code", "station_name", "species", "sex", "determinand", "matrix", "basis", 
    "unit", "methodUncertainty", "qflag", "sub.sample", "AMAP_group")] <- "text"
  col_types[id %in% c(
    "year", "value", "uncertainty", "limit_detection", "limit_quantification")] <- "numeric"

  AMAP_data <- read_excel(file_id, col_types = col_types, na = "")
  AMAP_data <- as.data.frame(AMAP_data)


  # populate qalink with 0 (generic link if no info) and alabo = "amap" (distinct because
  # all genuine alabo values are upper case)
  
  AMAP_data <- mutate(
    AMAP_data, 
    qalink = 0,
    alabo = "amap"
  )
  

  # populate pargroup - need to get this from info file in future
  
  if (!all(AMAP_data$determinand %in% c("HG", "DRYWT%", "LNMEA", "LIPIDWT%")))
    stop("need to code pargroup flexibly")

  AMAP_data <- mutate(
    AMAP_data, 
    pargroup = recode(
      .data$determinand, 
      "HG" = "I-MET", 
      "DRYWT%" = "B-BIO", 
      "LNMEA" = "B-BIO", 
      "LIPIDWT%" = "B-BIO")
  )
  
  
  # data checks
  # ensure no missing values in key variables
  
  var_id <- c("country", "station_code", "station_name", "year", "species", "determinand", "matrix", 
    "unit", "value", "sub.sample")

  ok <- sapply(AMAP_data[var_id], function(x) !any(is.na(x)))
    
  if (!all(ok))
    stop("following variables have missing values: ", paste(var_id[!ok], collapse = ", "))
  
  stopifnot(
    with(AMAP_data, ifelse(!is.na(uncertainty), !is.na(methodUncertainty), TRUE))
  )
  

  # ensure sub.sample is correctly populated
  
  wk <- with(AMAP_data, table(sub.sample, determinand, matrix))
  if (!all(wk %in% c(0, 1)))
    stop("replicate measurements in the same sub.sample, matrix combination are not allowed")
    
  ss_check <- unique(AMAP_data[c("sub.sample", "station_code", "station_name", "year", "species")])
  wk <- table(ss_check$sub.sample)
  if (!all(wk == 1)) {
    dups <- filter(ss_check, sub.sample %in% names(wk)[wk > 1])
    dups <- arrange(dups, sub.sample)
    cat("Error: the following sub.sample identifiers are duplicated:\n")
    print(dups)
    stop(call. = FALSE)
  }
  
  wk <- substring(AMAP_data$sub.sample, 1, 1)
  if (!all(wk %in% toupper(letters)))
    stop("sub.sample must start with a character")

  
  # ensure consistency of AMAP_group within sub.samples
  
  wk <- with(AMAP_data, tapply(AMAP_group, sub.sample, function(x) length(unique(x))))
  if (!all(wk %in% 1L))
    stop("multiple AMAP_groups in the same sub.sample are not allowed")
  
  
  # some timeseries have AMAP_group classifications for some but not all records
  # classify the latter as "undefined"
  
  AMAP_data <- group_by(AMAP_data, station_code, species, determinand, matrix)
  
  AMAP_data <- mutate(
    AMAP_data, 
    AMAP_group = if (all(is.na(AMAP_group))) AMAP_group else replace_na(AMAP_group, "undefined")
  )
  
  AMAP_data <- ungroup(AMAP_data)
  
  AMAP_data <- as.data.frame(AMAP_data)
  
  
  # read in AMAP stations

  file_id <- file.path(path, AMAP_stations)
  
  id <- read_excel(file_id, n_max = 0)
  id <- names(id)
  
  # check all required columns are present
  
  required <- c(
    "code", "OSPAR_region", "OSPAR_subregion", "country", "programGovernance", "station", "latitude", 
    "longitude") 
  
  ok <- required %in% id
  
  if (!all(ok)) 
    stop("following columns missing from AMAP station file: ", paste(required[!ok], collapse = ", "))
  

  col_types <- rep("skip", length(id))
  col_types[
    id %in% c("code", "OSPAR_region", "OSPAR_subregion", "country", "programGovernance", "station")] <- "text"
  col_types[id %in% c("latitude", "longitude")] <- "numeric"

  AMAP_stations <- readxl::read_excel(file_id, col_types = col_types, na = "")
  AMAP_stations <- as.data.frame(AMAP_stations)


  # data checks
  
  if (any(is.na(AMAP_stations)))
    stop("missing values are not allowed in the station file")
  
  wk <- substring(AMAP_stations$code, 1, 1)
  if (!all(wk %in% toupper(letters)))
    stop("the station code must start with a character")

  if (!all(AMAP_stations$OSPARregion %in% c(0, 1)))
    stop("OSPARregion must be 0 or 1")
  
  purpose <- ICES_data$info$purpose
  if (!all(grepl(purpose, AMAP_stations$programGovernance)))
    stop("stations must be marked for ", purpose)

  region_id <- c(
    "Barents Sea", "Greenland-Scotland ridge", "Norwegian Sea", "East of Iceland", "West of Iceland")
  wk <- with(AMAP_stations, ifelse(OSPAR_region %in% 1, OSPAR_subregion %in% region_id, TRUE))
  if (!all(wk))
    stop("region in OSPARregion 1 must be one of", paste(region_id, collapse = ", "))

  
  # check for no conflicts in station dictionary
  
  id <- c("country", "station")
  if (any(my_paste(AMAP_stations[id]) %in% my_paste(ICES_data$stations[id])))
    stop("AMAP station replicates a station in the ICES station dictionary", call. = FALSE)
  
  
  # check for consistency between data and stations

  id <- sort(unique(AMAP_data$station_code))
  ok <- id %in% c(AMAP_stations$code, as.character(ICES_data$stations$code))
  if (!all(ok))
    stop(
      "The following stations are in the data file but not the station file or the ICES station ", 
      "dictionary: \n       ", 
      paste(id[!ok], collapse = ", "), 
      call. = FALSE
    )
  
  stopifnot(
    AMAP_data$station_name %in% c(AMAP_stations$station, as.character(ICES_data$stations$station))
  )
    
  wk <- unique(AMAP_data[c("country", "station_code", "station_name")])
  ok <- my_paste(wk) %in% 
    c(my_paste(AMAP_stations[c("country", "code", "station")]), 
      my_paste(ICES_data$stations[c("country", "code", "station")]))

  if (!all(ok)) {
    cat("Error: the following combinations in the data are not in ICES or AMAP station files:\n")
    print(wk[!ok, ])
    stop(call. = FALSE)
  }
  

  # fill in PURPM and dataType
  
  AMAP_stations <- mutate(
    AMAP_stations,
    PURPM = "T",
    dataType = "CF"
  )
    

  # construct QA file
    
  AMAP_QA <- data.frame(qalink = 0)
  
  
  # bind with ICES data
  
  ICES_data$data <- bind_rows(ICES_data$data, AMAP_data)

  ICES_data$stations <- bind_rows(ICES_data$stations, AMAP_stations)

  ICES_data$QA <- bind_rows(ICES_data$QA, AMAP_QA)

  if (keep == "all") 
    return(ICES_data)
  
  # select only AMAP data
  
  id <- with(ICES_data$data, substring(sub.sample, 1, 1) %in% toupper(letters))
  ICES_data$data <- ICES_data$data[id, ]
  
  ICES_data
}