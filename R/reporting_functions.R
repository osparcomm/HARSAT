# mapping functions ----

ctsm.Mercator <- function(x) atanh(sin(x * pi / 180))

#' Calculates a polar projection
#' 
#' Calculates easting and northing in Lambert Azimuthal Equal 
#' Area North Polar aspect projection
#' 
#' @param latitude the latitude
#' @param longitude the longitude
#' @return a list with projected `longitude` (easting) and `latitude` (northing) values
#' @export
ctsm.projection <- function(latitude, longitude) {
  # calculates easting and northing in Lambert Azimuthal Equal Area North Polar aspect projection
  
  # constants and functions required

  radian <- function(theta) pi * theta / 180

  phi0 <- radian(52)
  phi0 <- radian(90)
  lambda0 <- radian(10)
  lambda0 <- radian(0)

  a <- 6378137
  e <- 0.081819191
  FE <- 0
  FN <- 0


  qcalc <- function(phi)
  {
    x <- sin(phi)
    (1 - e^2) * (x / (1 - (e * x)^2) - (1 / (2 * e)) * log((1 - e * x) / (1 + e * x)))
  }


  # calculations

  phi <- radian(latitude)
  lambda <- radian(longitude)
  
  q <- qcalc(phi)
  q0 <- qcalc(phi0)
  qP <- qcalc(pi / 2)
 
  Rq <- a * sqrt(qP / 2)
  beta <- asin(q / qP)
  beta0 <- asin(q0 / qP)
  rho <- a * sqrt(qP - q)
 
  B <- Rq * sqrt(2 / (1 + sin(beta0) * sin(beta) + cos(beta0) * cos(beta) * cos(lambda - lambda0)))
  D <- a * cos(phi0) / (sqrt(1 - e^2 * sin(phi0)^2) * Rq * cos(beta0))

#  E <- FE + B * D * cos(beta) * sin(lambda - lambda0)
#  N <- FN + (B / D) * (cos(beta0) * sin(beta) - sin(beta0) * cos(beta) * cos(lambda - lambda0))
 
  E <- FE + rho * sin(lambda - lambda0)
  N <- FN - rho * cos(lambda - lambda0)
   
  list(longitude = E, latitude = N)
}  


# support functions ----

#' Subsets an assessment object
#' 
#' Selects specific time series and simplifies the data, stations and 
#' assessment components to match
#'  
#' @param assessment_obj An assessment object resulting from a call to
#'   run_assessment.
#' @param subset A vector specifying the timeseries to be retained. An
#'   expression will be evaluated in the timeSeries component of assessment_obj; 
#'   use 'series' to identify individual timeseries.
#'
#' @returns a new assessment object, after applying the subset
#.
#' @export
subset_assessment <- function(assessment_obj, subset) {
  
  # reporting_functions.R
  # subsets an assessment object by filtering on the timeSeries component
  
  timeSeries <- assessment_obj$timeSeries
  
  timeSeries <- tibble::rownames_to_column(timeSeries, "series")
  ok <- eval(substitute(subset), timeSeries, parent.frame())
  timeSeries <- timeSeries[ok, ]
  series_id <- timeSeries$series
  
  row.names(timeSeries) <- NULL
  timeSeries <- tibble::column_to_rownames(timeSeries, "series")
  
  assessment_obj$timeSeries <- timeSeries

  
  # update other components to be consistent
  
  assessment_obj$assessment <- assessment_obj$assessment[series_id]
  
  ok <- assessment_obj$data$seriesID %in% series_id
  assessment_obj$data <- assessment_obj$data[ok, ]
  
  ok <- assessment_obj$stations$station_code %in% timeSeries$station_code 
  assessment_obj$stations <- assessment_obj$stations[ok, ]
  row.names(assessment_obj$stations) <- NULL
  
  assessment_obj
}  




ctsm.web.AC <- function(assessment_ob, classification) {
  
  # identifies which AC are used for each determinand

  assessment <- assessment_ob$assessment
  
  # gets series ID for each timeseries by determinand
  
  assessment_id <- split(
    rownames(assessment_ob$timeSeries), 
    assessment_ob$timeSeries$determinand, 
    drop = TRUE
  )
  
  # identity all AC that are relevant
  
  AC_id <- names(classification[["below"]])
  stopifnot(AC_id %in% assessment_ob$info$AC)
  
  # loop over determinands

  out <- sapply(assessment_id, USE.NAMES = TRUE, simplify = FALSE, FUN = function(id) {

    # AC used by series
  
    AC_series <- lapply(assessment[id], function(i) {
      AC <- i$AC
      AC <- AC[AC_id]
      AC <- !is.na(AC)
    })
    
    AC_series <- dplyr::bind_rows(AC_series)
    
    
    # AC used by determinand
    
    AC_used <- apply(AC_series, 2, any)
    
    
    # is BAC the only AC for some series
    
    if ("BAC" %in% AC_id) {
      BAC_only <- rowSums(AC_series) == 1 & AC_series[["BAC"]]
      BAC_only <- any(BAC_only)
      AC_used <- c(AC_used, "BAC_only" = BAC_only)
    }  
    
    
    # are there some series with no AC
    
    AC_none <- rowSums(AC_series) == 0
    AC_none <- any(AC_none)

    c(AC_used, "none" = AC_none)
  })
  
  out <- dplyr::bind_rows(out, .id = "determinand")
  
  out <- as.data.frame(out)
  
  tibble::column_to_rownames(out, "determinand")
}




# summary table ----

#' Write assessment summary to a csv file
#'
#' @description
#'
#' Creates a data frame summarising the assessment of each time series and
#' writes it to a csv file. The summary includes:
#'
#' * meta-data such as the monitoring location and number of years of data for
#' each time series
#' * the fitted values in the last monitoring year with associated upper
#' one-sided 95% confidence limits
#' * the trend assessments (p-values and trend estimates)
#' * the status assessments (if there any thresholds)
#' * (optionally) a symbology summarising the trend (shape) and status (colour)
#' of each time series. This is experimental. 
#'
#' @param harsat_obj A harsat object following a call to `run_assessment`.
#' @param output_file The name of the output csv file. If using NULL, the file
#'   will be called `biota_summary.csv`, `sediment_summary.csv` or `water_summary.csv`
#'   as appropriate. By default the file will be written to the working
#'   directory. If a file name is provided, a path to the output file can also
#'   be provided (e.g. using `file.path`). The `output_dir`` option can also be
#'   used to specify the output file directory.
#' @param output_dir The output directory for `output_file`. The default is the
#'   working directory. Any file path provided in `output_file`, will be
#'   appended to `output_dir`. The resulting output directory must already
#'   exist.
#' @param export Logical. `TRUE` (the default) writes the summary table to a csv
#'   file. `FALSE` returns the summary table as an R object (and does not write to
#'   a csv file).
#' @param threshold_groups A names list of valid thresholds that allows 
#'   thresholds of the same 'type' to be reported together. See details.
#' @param collapse_AC `r lifecycle::badge("deprecated")` Use `threshold_groups` 
#'   instead.  
#' @param extra_output A character vector specifying extra summary metrics 
#'   to be included in the output. Currently only recognises "power" to give the 
#'   seven power metrics computed for lognormally distributed data. Defaults to 
#'   `NULL`; i.e. no extra output. 
#' @param symbology Experimental. A character string "default" or a user-defined 
#'   function that specifies a symbology typically used to characterise the 
#'   patterns of change in and the status of each time series. Defaults to 
#'   `NULL`; i.e. no symbology. Multiple symbologies can be applied. See 
#'   details.
#' @param symbology_control Experimental. A named list of control options for 
#'   the symbology. See details.
#' @param determinandGroups optional, a list specifying `labels` and `levels`
#'   to rename the existing determinand groups. The life of this argument is 
#'   limited.
#' @param timestamp Logical. `FALSE` (the detault) does nothing. `TRUE` outputs
#'   an extra column with the date and time that the summary file was created. 
#'   Most likely only useful if you are using `append` (see below). 
#' @param append Logical. `FALSE` (the default) overwrites any existing summary
#'   file. `TRUE` appends data to it, creating it if it does not yet exist.
#'
#' @returns a summary object, when `export` is `FALSE`
#' 
#' @section Default symbology:
#' 
#' `symbology = "default"` calls a pre-defined symbology that generates a 
#' 'shape' and a 'colour' to characterise the status of each time series. Its
#' behaviour is controlled using `symbology_control`, a named list with the 
#' following elements:
#'
#' * `shape`: a list with names `none`, `mean`, `flat`, `up`, `down` giving the
#' shape associated with each pattern of change. Here, `none` corresponds to 
#' insufficient data to fit a parametric model; `mean` to sufficient data to fit
#' a parametric model but not to assess for trends; `flat` to no significant 
#' change in level (concentration) over time; `up` to a #' significant increase 
#' in level over time; `down` to a significant decrease in level over time. 
#' Their default values are `"small_open_circle"`, 
#' `"small_filled_circle"`, `"large_filled_circle"`, `"upward_triangle"`,
#' `"downward_triangle"`
#' * `alpha` is the size of the test for change; default = `0.05`
#' * `change` determines whether the change is based on the recent time window 
#' (typically the last twenty years) or the whole time series; options 
#' `"recent"` (default) and `"overall"` 
#' * `colour` (default `NULL`) is a named list that characterises the status of
#' a time series based on specified thresholds. The list names must match 
#' (a subset of) the names of the thresholds used in the assessment or, if the 
#' thresholds have been grouped, the group names. Each threshold must have two 
#' elements: `below` gives the colour if the time series is 
#' significantly below the threshold (p < 0.05); `above` gives the colour 
#' otherwise. If multiple thresholds are used, they must be ordered from best
#' to worst status. See examples
#' * `no_threshold` (default `"black"`) is the colour used when no thresholds are
#' applied to a time series. Another option might be to use `NA_character_`
#' * `adjust_nonparam` (default `TRUE`) is a logical that allows the symbology
#' to be adjusted for short time series (often dominated by less-than values)
#' where a non-parametric test for status can be applied
#' * `names` (default `list(colour = "colour", shape = "shape")`) allows the 
#' names of the symbology columns in the summary table to be adjusted; this can
#' be important if multiple symbologies are applied
#' 
#' @section Custom symbologies:
#' 
#' Users can apply custom symbologies by letting `symbology` be a user-supplied
#' function of the form `fn(summary, info, control)` where:
#' 
#' * `summary` is the summary table before applying the symbology; for 
#' convenience, there is an additional column `method` which doesn't 
#' appear in the final summary table but takes, in particular, values `"none"` 
#' and `"mean"` corresponding respectively to no parametric model and 
#' insufficient data to fit a trend (see `shape` above)
#' * `info` contains the contents of `harsat_obj$info`; i.e. all the reference
#' tables and additional information about the assessment. For convenience,
#' it also contains a temporary element `.threshold_group` which has the names
#' of the threshold groups (which can differ from the thresholds themselves)
#' * `control` contains `symbology_control` and allows the user to pass 
#' additional information to the function
#' 
#' The output of the function must be a data frame with one column called
#' `series` that contains the series identifier of each time series 
#' (not necessarily in the same order as the summary table) and  
#` the remaining columns giving the symbology. There can be any number of 
#' symbology columns, but they must not share any of the existing names in the 
#' summary table.
#' 
#' See the examples for more inspiration.
#' 
#' @section Multiple symbologies:
#' 
#' Multiple symbologies can be applied by specifying a named list whch can be 
#' a mixture of default symbologies and custom symbologies. `symbology_control`
#' must then be a named list (with the same names) giving control information 
#' for each symbology. It is important to ensure that each symbology gives 
#' output columns with different names. See examples.
#' 
#' @examples
#' 
#' # Default symbology with one threshold: the EQS. The colour will be "green"
#' # if the time series is significantly below the EQS in the last monitoring
#' # year and "red" otherwise 
#' \dontrun{
#' write_summary_table(
#'   water_assessment,
#'   symbology = "default",
#'   symbology_control = list(
#'     colour = list(EQS = list(below = "green", above = "red"))
#'   )
#' )
#' }
#'
#' # Now applied using the overall change instead of the recent change.
#' \dontrun{
#' write_summary_table(
#'   water_assessment,
#'   symbology = "default",
#'   symbology_control = list(
#'     colour = list(EQS = list(below = "green", above = "red")),
#'     change = "overall"
#'   )
#' )
#' }
#' 
#' # If we only want to change one shape, then we only need to specify that one
#' \dontrun{
#' write_summary_table(
#'   water_assessment,
#'   symbology = "default",
#'   symbology_control = list(
#'     colour = list(EQS = list(below = "green", above = "red")),
#'     shape = list(flat = "square")
#'   )
#' )
#' }
#'
#' # Assessment thresholds grouped into BAC and EAC equivalents.
#' # Symbology now has two thresholds giving:
#' # "blue" if significantly below the BAC
#' # "orange" if not significantly below the BAC and there is no EAC
#' # "green" if significantly below the EAC but not the BAC
#' # "red" otherwise  
#' \dontrun{
#' write_summary_table(
#'   sediment_assessment,
#'   threshold_groups = list(
#'     BAC = "BAC", 
#'     EAC = c("EAC", "ERL", "EQS", "FEQG")
#'   ),
#'   symbology = "default", 
#'   symbology_control = list(
#'     colour = list(
#'       BAC = list(below = "blue", above = "orange"),
#'       EAC = list(below = "green", above = "red")
#'     ), 
#'   )
#' )
#' }
#' 
#' # Assessment thresholds grouped into BAC and EAC equivalents. Human health
#' # thresholds grouped as HQS
#' # Two symbologies applied, one for environmental thresholds, the other for
#' # health thresholds
#' # Note the named lists and that the output names are specified 
#' \dontrun{
#' write_summary_table(
#'   biota_assessment,
#'   threshold_groups = list(
#'     BAC = c("BAC", "NRC"),
#'     EAC = c("EAC", "FEQG", "LRC", "QSsp"), 
#'     HQS = c("MPC", "QShh")
#'   ),
#'   symbology = list(env = "default", health = "default"), 
#'   symbology_control = list(
#'     env = list(
#'       colour = list(
#'         BAC = list(below = "blue", above = "orange"),
#'         EAC = list(below = "green", above = "red")
#'       ), 
#'       names = list(shape = "shape_env", colour = "colour_env")
#'     ),
#'     health = list(
#'       colour = list(HQS = list(below = "green", above = "red")), 
#'       names = list(shape = "shape_health", colour = "colour_health")
#'     )
#'   )
#' )
#' }
#' 
#' # Custom symbology that only reports time series where there is sufficient
#' # information to assess trends and which colours the time series by whether, 
#' # for each determinand, mean concentrations in the last monitoring year are 
#' # below or above the median concentration observed across time series 
#' \dontrun{
#' symbology_user <- function(summary, info, control) {
#'   summary <- dplyr::mutate(
#'     summary,
#'     shape = dplyr::case_when(
#'       is.na(p_overall_change)   ~ NA_character_,
#'       p_overall_change > 0.05   ~ "circle",
#'       overall_change > 0        ~ "upward_triangle",
#'       overall_change < 0        ~ "downward_triangle"
#'     )
#'   )
#'   summary <- summary |> 
#'     dplyr::group_by(determinand) |>
#'     dplyr::mutate(
#'       .shape = !is.na(shape),
#'       colour = dplyr::case_when(
#'         !.shape                                          ~ NA_character_,
#'         mean_last_year <= median(mean_last_year[.shape]) ~ "blue",
#'        mean_last_year > median(mean_last_year[.shape])  ~ "red"
#'       )
#'     ) |>
#'     dplyr::ungroup()
#'   summary[c("series", "shape", "colour")]
#' }
#'
#' write_summary_table(
#'   biota_assessment,
#'   symbology = symbology_user
#' )
#' }
#' 
#'
#' @export
write_summary_table <- function(
  harsat_obj, 
  output_file = NULL, 
  output_dir = ".", 
  export = TRUE,
  threshold_groups = NULL, 
  collapse_AC = lifecycle::deprecated(), 
  extra_output = NULL, 
  symbology = NULL, 
  symbology_control = list(),
  determinandGroups = NULL, 
  timestamp = FALSE, 
  append = FALSE) {

  if (lifecycle::is_present(collapse_AC)) {
    lifecycle::deprecate_warn(
      "1.0.4", 
      "write_summary_table(collapse_AC)", 
      "write_summary_table(threshold_groups)")
      threshold_groups <- collapse_AC
  }
  
  
  # get summary from assessment, combine with timeseries object, and rename
  # variables into more reader_friendly form
  
  summary <- make_summary_table(harsat_obj, extra_output, determinandGroups)


  # group thresholds to simplify summary table
  # also updates list of thresholds to work with in the symbology stage
  # note that the names of the thresholds might have changed so can no longer
  #   be picked up from info$AC

  thresholds <- harsat_obj$info$AC

  if (!is.null(threshold_groups)) {
    tmp <- group_thresholds(summary, threshold_groups, thresholds)
    summary <- tmp$summary
    thresholds <- tmp$thresholds
  }
  

  # apply symbology - shape and colour of plotting symbols

  if (!is.null(symbology)) {
    summary <- make_symbology(
      summary, 
      harsat_obj$info, 
      symbology, 
      symbology_control, 
      thresholds
    )
  }

  
  # remove 'method', a convenience variable for constructing symbologies, from 
  # summary table
  
  summary$method <- NULL
  
  
  # results
  
  # if export = FALSE return summary data frame

  if (!export) {
    return(summary)
  }
  

  # otherwise write to .csv file
  
  # get default output_file 
  
  if (is.null(output_file)) {
    output_file <- paste0(harsat_obj$info$compartment, "_summary.csv")
  } 
  
  # check output_file has valid extension
  
  if (!endsWith(output_file, ".csv")) {
    stop(
      "\nThe output file '", output_file, "' does not have a .csv extension.\n", 
      "Check the information supplied to argument 'output_file'.",
      call. = FALSE
    )
  }    
    
  # combine output_file and output_dir and check output directory exists
  
  output_file <- file.path(output_dir, output_file)
  
  wk <- dirname(output_file)
  if (!dir.exists(wk)) {
    stop(
      "\nThe output directory '", wk, "' does not exist.\n", 
      "Create it or check the information supplied to argument 'output_dir'",
      " is correct.",
      call. = FALSE
    )
  }

  
  # create timestamp if timestamp = TRUE

  if (timestamp) {  
    summary$timestamp <- Sys.time()
    summary <- dplyr::relocate(summary, timestamp)
  }    

  # headers on a new file aren't created if append = TRUE
  
  if (!file.exists(output_file)) {
    append <- FALSE
  }
  
  # if append = TRUE check that column names are identical and warn if there
  # are series that are going to be repeated
  
  if (append) {
    old_summary <- safe_read_file(output_file)
    if (!identical(names(old_summary), names(summary))) {
      stop(
        "\nCannot append because the names of the new summary table differ ",
        "from those of the\n", 
        "existing summary file.",
        call. = FALSE
      )
    }
    
    if (any(summary$series %in% old_summary$series)) {
      warning(
        "Some time series in the new summary table are already reported in ",
        "the existing\n", 
        "summary file: you should check what is going on.",
        call. = FALSE
      )
    }
  }
  
  readr::write_excel_csv(summary, output_file, na = "", append = append)
  
  return(invisible())
}


#' @export
make_summary_table <- function(harsat_obj, extra_output, determinandGroups) {
  
  # gets summary from each assessment, augments this with anything in 
  # extra_output, and merges with timeseries 

  info <- harsat_obj$info
  
  # merge stations with timeseries

  # get determinand group
  # NB detGroup is currently a character if determinandGroups is null and a 
  # factor otherwise - the use of the factor assists with ordering, but probably
  # not needed here as this is primarily an OHAT requirement - have raised
  # an issue to tidy this up

  # NB Use of determinandGroups has a limited life expectancy.
  
  timeseries <- harsat_obj$timeSeries
  
  timeseries <- tibble::rownames_to_column(timeseries, "series")
  
  timeseries$detGroup <- ctsm_get_info(
    info$determinand, 
    timeseries$determinand, 
    "group", 
    info$compartment, 
    sep = "_"
  )
  
  if (!is.null(determinandGroups)) {
    
    if (!all(timeseries$detGroup %in% determinandGroups$levels)) {
      stop('some determinand groups present in data, but not in groups argument')
    }
    
    timeseries$detGroup <- factor(
      timeseries$detGroup, 
      levels = determinandGroups$levels, 
      labels = determinandGroups$labels, 
      ordered = TRUE
    )
    
    timeseries$detGroup <-   timeseries$detGroup[, drop = TRUE]
  }  
  

  # join with stations

  timeseries <- dplyr::left_join(
    timeseries, 
    harsat_obj$stations,
    by = "station_code"
  )
  

  # order assessment so that it is compatible with timeseries - 

  assessment <- harsat_obj$assessment
  
  assessment <- assessment[timeseries$series]
  
  summary <- sapply(
    assessment, 
    function(x) {
      out <- x$summary 
      if ("power" %in% extra_output && !is.null(x$power)) {
        out <- cbind(out, x$power)
      }
      out$method <- x$method
      out
    }, 
    simplify = FALSE
  )
  
  if (any(is.null(summary)) | (length(summary) != nrow(timeseries))) {
    stop("coding error - contact harsat development team")
  }
  
  summary <- dplyr::bind_rows(summary)
  
  summary <- cbind(timeseries, summary)

  
  # rename variables 
  # much of this can be simplified by keeping these names consistent throughout 
  # the code
  
  summary <- dplyr::rename(
    summary, 
    determinand_group = "detGroup", 
    n_year_all = "nyall",
    n_year_fit = "nyfit",
    n_year_positive = "nypos",
    first_year_all = "firstYearAll",
    first_year_fit = "firstYearFit",
    last_year = "lastyear",
    detectable_trend = "dtrend",
    mean_last_year = "meanLY",
    climit_last_year = "clLY"
  )
  
  if ("power" %in% extra_output) {
    summary <- dplyr::rename(
      summary, 
      power_dt_obs = "dtrend_obs",
      power_dt_seq = "dtrend_seq",
      power_dt_ten = "dtrend_ten",
      power_ny_seq = "nyear_seq",
      power_pw_obs = "power_obs",
      power_pw_seq = "power_seq",
      power_pw_ten = "power_ten",
    )
  }

  if (!is.null(info$AC)) {
    thresholds <- info$AC

    pos <- match(thresholds, names(summary))
    names(summary)[pos] <- paste0(thresholds, "_value")
    
    for (suffix in c("diff", "achieved", "below")) {
      in_id <- paste0(thresholds, suffix)
      out_id <- paste0(thresholds, "_", suffix)
    
      pos <- match(in_id, names(summary))
      names(summary)[pos] <- out_id
    }
  }
    

  # reorder variables and sort 

  var_id <- c(
    "series", 
    info$region$id,  
    "country", 
    "station_code", "station_name", "station_longname", 
    "station_latitude", "station_longitude", "station_type", "waterbody_type", 
    "determinand", "determinand_group", "species", "filtration",
    "matrix", "basis", "unit", "sex", "method_analysis", 
    "normaliser", "normaliser_value", "normaliser_unit",
    "subseries"
  ) 
  
  summary <- dplyr::relocate(summary, dplyr::any_of(var_id))
  
  if ("power" %in% extra_output) {
    var_id <- c(
      "power_dt_obs", "power_dt_seq", "power_dt_ten", 
      "power_ny_seq", 
      "power_pw_obs", "power_pw_seq", "power_pw_ten"
    )
    summary <- dplyr::relocate(
      summary, 
      dplyr::all_of(var_id), .after = "detectable_trend") 
  }
  
  
  var_id <- c(
    info$region$id, "country", "station_name", "species", "determinand_group", 
    "determinand", "matrix", "filtration"
  )
  
  summary <- dplyr::arrange(
    summary, 
    dplyr::pick(dplyr::any_of(var_id))
  )

  
  # rename region variables if required
  
  if (!identical(info$region$id, info$region$names)) {
    pos <- match(info$region$id, names(summary))
    names(summary)[pos] <- info$region$names
  }
  

  summary
}


make_symbology <- function(
    summary, info, symbology, symbology_control, thresholds) {
  

  # add thresholds to info 
  # safer approach than picking up (guessing) them from the summary names 
  # means that they can be accessed by user-defined symbologies
  # thresholds might have been grouped so existing information in info might 
  #   not be sufficient
  
  if (".threshold_group" %in% names(info)) {
    stop(
      "\nwrite_summary_table wants to create .threshold_group in\n", 
      "harsat_obj$info but it already exists; this is unexpected, so check\n",
      "your script and then contact the harsat development team", 
      call. = FALSE
    )
  }
  
  info$.threshold_group <- thresholds
  

  # symbology can either be a single character string, a single function, 
  # or a named list of the above

  if (!is.list(symbology)) {
    
    symbology <- list(only = symbology)
    symbology_control <- list(only = symbology_control)
    
  } else {
    
    if (is.null(names(symbology))) {
      stop(
        "\nsymbology - multiple symbologies must be specified by a named ",
        "list",
        call. = FALSE
      )
    }
    
    if (is.null(names(symbology_control)) || 
        !identical(sort(names(symbology_control)), sort(names(symbology)))) {
      stop(
        "\nsymbology_control - multiple symbologies have been specified, so\n",
        "symbology_control must be a named list with the same names as ",
        "symbology",
        call. = FALSE
      )
    }
    
  }
  
  for (i in names(symbology)) {
    
    if (is.character(symbology[[i]])) {
      
      if (symbology[[i]] != "default") {
        stop(
          "\nsymbology: currently the only pre-defined symbology is \"default\"",
          call. = FALSE
        )
      }
      symbology_fn <- symbology_default
      
    } else if (is.function(symbology[[i]])) {
      
      id <- formalArgs(symbology[[i]])
      if (!identical(id, c("summary", "info", "control"))) {
        stop(
          "\nsymbology: user-defined function arguments should be 'summary'",
          " 'info' and\n'control'",
          call. = FALSE
        )
      }
      
      symbology_fn <- symbology[[i]]
      
    } else {
      
      stop(
        "\nsymbology: symbology is not a recognised character string or a\n", 
        "user-defined function",
        call. = FALSE
      )
      
    }
    
    result <- symbology_fn(summary, info, symbology_control[[i]])
    
    if (!("series" %in% names(result))) {
      stop(
        "\nerror in applying symbology: series not in result",
        call. = FALSE
      )
    }
    
    id <- setdiff(names(result), "series")
    
    if (any(id %in% names(summary))) {
      stop(
        "\nsymbology: the names of the symbology variables already exist;\n",
        "if applying the default symbology, adjust control$names;\n", 
        "if applying a user-defined symbology, adjust the output names;\n",
        "if applying multiple symbologies, check you are not using the same\n",
        "names twice", 
        call. = FALSE
      )
    }
    
    summary <- dplyr::left_join(
      summary, 
      result, 
      by = "series", 
      relationship = "one-to-one"
    )
    
    summary <- dplyr::relocate(
      summary, 
      dplyr::all_of(id), 
      .before = "n_year_all"
    )
    
  }    
  
  summary
}



symbology_default <- function(
    summary, 
    info, 
    control = list()) {
  
  # silence non-standard evaluation warnings

  .data <- NULL
  

  # set up and modify control structures
  
  control = symbology_default_cntrl(control, info)
    
  
  # shape = trend 
  
  # a trend is estimated unless method is "none" or "mean"
  # note, method is inconsistently applied apart from "none" and "mean"; look 
  #   at imposex assessments - this needs to be resolved
  
  # trend can either be 'overall', based on p_overal_change, or 'recent', based
  #   on 'p_recent_change' 
  # note p_recent_change might not exist even if a trend has been fitted to the 
  #   whole time series if there are too few years of data in the recent window

  trend_id <- paste0(control$change, "_change")
  trend_p <- paste0("p_", trend_id)
    
  summary <- dplyr::mutate(
    summary,
    shape = dplyr::case_when(
      .data$method %in% "none"           ~ control$shape$none,
      .data$method %in% "mean"           ~ control$shape$mean,
      is.na(.data[[trend_p]])            ~ control$shape$flat,
      .data[[trend_p]] >= control$alpha  ~ control$shape$flat,
      .data[[trend_id]] > 0              ~ control$shape$up,
      .data[[trend_id]] < 0              ~ control$shape$down
    )
  )
  
  
  # colour = status
  
  # status is based on  
  # - upper confidence limit (if a parametric model has been fitted) 
  # - mean_last_year (otherwise)
  
  if (!is.null(control$colour)) {
    
    thresholds <- names(control$colour)
    
    # need to ensure thresholds are ordered correctly from best to worst status
    # gets complicated if good_status equates to high concentrations for some
    #   determinands; to get around this simply change sign on the value of the 
    #   threshold when good_status = high

    good_status <- ctsm_get_info(
      info$determinand, 
      summary$determinand, 
      "good_status"
    )
    
    t_value <- paste0(thresholds, "_value")    
    
    wk <- summary[t_value]
    wk[] <- lapply(wk, "*", dplyr::if_else(good_status == "low", 1, -1))
    
    ok <- apply(wk, 1, function(x) {
      if (all(is.na(x))) return(TRUE)
      x <- x[!is.na(x)]
      length(x) == 1L || all(diff(x) > 0) 
    }) 
    
    if (!all(ok)) {
      stop(
        "/nsymbology control:\n",
        "to specify the colour correctly, the thresholds must be ordered from\n",
        "best to worst status;\n",
        "if this has already been done, then this error suggests there are\n",
        "inconsistencies in the threshold values applied to this assessment",
        call. = FALSE
      )
    }   
    
    # when good_status equates to low concentrations, negative t_diff is good
    #   since t_diff = clLY - threshold < 0
    # when good status equates to high concentrations, positive t_diff is good
    # to make colour calculation 'simple' change sign on t_diff when 
    #    goodStatus == high
    
    t_diff <- paste0(thresholds, "_diff")
    
    wk <- summary[t_diff]
    wk[] <- lapply(wk, "*", dplyr::if_else(good_status == "low", 1, -1))
    
    summary$colour <- apply(wk, 1, function(x) {
      
      if (all(is.na(x))) return(control$no_threshold)
      
      thresholds <- thresholds[!is.na(x)]
      x <- x[!is.na(x)]
      
      if (any(x < 0)) {
        id <- thresholds[which.max(x < 0)]
        control$colour[[id]]$below
      } else {  
        id <- thresholds[length(x)]
        control$colour[[id]]$above
      }
      
    })  

  } else {
    
    summary$colour <- control$no_threshold
    
  }
  
  
  # adjust shape and colour for nonparametric test method = "none" 
  
  if (control$adjust_nonparam & !is.null(control$colour)) {
    
    # get the variables that contain the non-parametric test result 
    # note these might not exist if e.g. only biological effects have been
    # assessed

    t_below <- paste0(thresholds, "_below")
    
    ok <- t_below %in% names(summary)
    
    if (any(ok)) {
      
      wk <- dplyr::select(summary, dplyr::any_of(t_below))
      
      # adjust good_status = high as before
      
      wk[] <- lapply(wk, function(x) {
        id <- !is.na(x) & good_status == "high"
        if (any(id)) {
          x[id] <- dplyr::if_else(x[id] == "below", "above", "below")
        }
        x
      })
      
      wk <- apply(wk, 1, function(x) {
        
        if (all(is.na(x))) return(NA_character_)
        
        thresholds <- thresholds[!is.na(x)]
        x <- x[!is.na(x)]
        
        if (any(x == "below")) {
          id <- thresholds[which.max(x == "below")]
          control$colour[[id]]$below
        } else {
          id <- thresholds[length(x)]
          control$colour[[id]]$above
        }
      })  
      
      # only apply this where method = "none" - i.e. no parametric model
      
      id <- summary$method %in% "none" & !is.na(wk)
      summary$colour[id] <- wk[id]
      summary$shape[id] <- control$shape$mean
    }
    
  }

  out <- summary[c("series", "shape", "colour")]

  names(out)[2:3] <- c(control$names$shape, control$names$colour)
    
  out
}


symbology_default_cntrl <- function(control, info) {

  default <- list(
    # trends
    shape = list(
      none = "small_open_circle", 
      mean = "small_filled_circle", 
      flat = "large_filled_circle", 
      up = "upward_triangle", 
      down = "downward_triangle"
    ),
    change = "recent", 
    alpha = 0.05,

    # status
    colour = NULL,
    no_threshold = "black",
    
    # other
    adjust_nonparam = TRUE,
    names = list(
      colour = "colour",
      shape = "shape"
    )
  )
  
  control <- modifyList(default, control, keep.null = TRUE)

  
  # trend  

  ok <- length(control$change) == 1L && 
    control$change %in% c("recent", "overall")
  if (!ok) {
    stop(
      "\nsymbology_control - change:\n", 
      "change must be either 'recent' or 'overall'", 
      call. = FALSE
    )
  }
  
  ok <- length(control$shape) == 5L &&
    all(c("none", "mean", "flat", "up", "down") %in% names(control$shape))
  if (!ok ) {
    stop(
      "\nsymbology_control - shape:\n", 
      "shape must be a list of 5 elements with names 'none', 'mean', 'flat'\n", 
      "'up', 'down'; each element of the list must be a character string", 
      call. = FALSE
    )
  }

  
  if (control$alpha <= 0 || control$alpha > 0.5) {
    stop(
      "\nsymbology_control - alpha:\n", 
      "alpha must be numeric and between 0 and 0.5; \n", 
      "typically it would be 0.1, 0.05 or 0.01", 
      call. = FALSE
    )
  }  
  
  
  # status

  ok <- length(control$no_threshold) == 1L && is.character(control$no_threshold)
  if (!ok) {
    stop(
      "\nsymbology_control - no_threshold:\n", 
      "no_threshold must be a character string giving the colour used for time\n",
      "series with no thresholds; to omit these time series from the colour\n", 
      "symbology, use 'no_threshold = NA_character_'", 
      call. = FALSE
    )
  }
  
  
  if (!is.null(control$colour)) {
    
    if (is.null(info$.threshold_group)) {
      stop(
        "\nsymbology_control - colour:\n", 
        "colour has been specified, but no thresholds were used", 
        call. = FALSE
      )
    }
    
    ok <- all(names(control$colour) %in% info$.threshold_group)
    if (!ok) {
      stop(
        "\nsymbology_control - colour:\n", 
        "the thresholds are not all recognised;\n",
        "if the thresholds have been grouped (using argument ", 
        "'threshold_group'),\nthen check the ",
        "thresholds specified in colour match the group names",
        call. = FALSE
      )
    }
    
    ok <- sapply(
      control$colour, 
      function(x) identical(sort(names(x)), c("above", "below"))
    )
    if (!all(ok)) {
      stop(
        "\nsymbology_control - colour:\n", 
        "each threshold must have two colours specified as a vector with \n",
        "names `below` and `above` (in either order)",
        call. = FALSE
      )
    }

    
    not_ok <- control$no_threshold %in% unlist(control$colour)
    if (not_ok) {
      stop(
        "\nsymbology_control - no_threshold and colour:\n", 
        "no_threshold must be different to the colours specified for each \n", 
        "threshold", 
        call. = FALSE
      )
    }

  }


  control
}



group_thresholds <- function(summary, threshold_groups, thresholds) {

  # NB not all thresholds need to be specified in threshold_groups; only those
  # where some grouping is going on
    
  # lots of error trapping to catch anything stupid!

  if (is.null(thresholds)) {
    stop(
      "threshold_groups specified but no thresholds used in the assessment", 
      call. = FALSE
    )
  }    

  if (is.null(names(threshold_groups))) {
    stop(
      "threshold_groups must be a named list", 
      call. = FALSE
    )
  }

  id <- names(threshold_groups) 
  ok <- !duplicated(id)
  if (!all(ok)) {
    stop(
      "cannot give different threshold_groups the same name",
      call. = FALSE
    )
  }
    
  id <- unlist(threshold_groups)
  ok <- !duplicated(id)
  if (!all(ok)) {
    id <- id[!ok]
    id <- sort(unique(id))
    stop(
      "these thresholds are specified more than once in threshold_groups:\n",
      paste(id, collapse = "; "),
      call. = FALSE
    )
  }
  
  id <- unlist(threshold_groups)
  ok <- id %in% thresholds
  if (!all(ok)) {
    id <- id[!ok]
    stop(
      "these thresholds are specified in threshold_groups but not ",
      "used in the assessment:\n", 
      paste(id, collapse = ", "), 
      call. = FALSE
    )
  }
  
  
  # check that no threshold is used as a group name for a 'different' group
  # i.e. one that does not contain that threshold
  
  ok <- sapply(
    names(threshold_groups), 
    function(id) {
      !id %in% thresholds || id %in% threshold_groups[[id]]
    }  
  )
    
  if (!all(ok)) {
    id <- names(threshold_groups)[!ok]
    stop(
      "these threshold groups do not contain the threshold of the same name:\n ", 
      paste(id, collapse = ", "), 
      call. = FALSE
    )
  }
  
  
  # add ungrouped thresholds to threshold_groups to keep everything simple!

  ok <- thresholds %in% unlist(threshold_groups)
  if (!all(ok)) {
    extra <- thresholds[!ok]
    names(extra) <- extra
    extra <- as.list(extra)

    threshold_groups <- c(threshold_groups, extra)
  }


  threshold_groups <- threshold_groups[sort(names(threshold_groups))]

  
  # set up type variables for each threshold

  id <- paste0(thresholds, "_type")
  summary[id] <- lapply(thresholds, function(x) {
    x_id <- paste0(x, "_value")
    dplyr::if_else(is.na(summary[[x_id]]), NA_character_, x)
  })


  suffix <- c("_type", "_value", "_diff", "_achieved", "_below") 
  
  for (group_id in names(threshold_groups)) {
    
    threshold_id <- threshold_groups[[group_id]]
    
    # if only one threshold, then the only thing to do is rename it if it
    # has a different name to the group name
    
    if (length(threshold_id) == 1L && threshold_id != group_id) {
      old_names <- paste0(threshold_id, suffix) 
      new_names <- paste0(group_id, suffix)
      
      pos <- match(old_names, names(summary))
      names(summary)[pos] <- new_names
    }
    
    
    if (length(threshold_id) > 1L) {

      # check multiple thresholds in the group haven't been applied to the same
      # time series
      
      id <- paste0(threshold_id, "_type")
      ok <- apply(summary[id], 1, function(x) {sum(!is.na(x)) <= 1L})
      if (!all(ok)) {
        series_id <- summary$series[which.min(ok)]
        stop(
          "cannot create the ", group_id, " threshold group because more than ", 
          "one component threshold\n", 
          "has been applied to the same series. For example, see ",
          "series:\n", series_id, 
          call. = FALSE
        )
      }
      
      in_id <- paste0(threshold_id, "_type")
      out_id <- paste0(group_id, "_type")
      summary[out_id] <- apply(summary[in_id], 1, group_thresholds_engine)
      
      in_id <- paste0(threshold_id, "_value")
      out_id <- paste0(group_id, "_value")
      summary[out_id] <- apply(summary[in_id], 1, group_thresholds_engine)
      
      in_id <- paste0(threshold_id, "_diff")
      out_id <- paste0(group_id, "_diff")
      summary[out_id] <- apply(summary[in_id], 1, group_thresholds_engine)
      
      in_id <- paste0(threshold_id, "_achieved")
      out_id <- paste0(group_id, "_achieved")
      summary[out_id] <- apply(summary[in_id], 1, group_thresholds_engine)
      
      in_id <- paste0(threshold_id, "_below")
      out_id <- paste0(group_id, "_below")
      summary[out_id] <- apply(summary[in_id], 1, group_thresholds_engine)

      # remove redundant columns
      
      id <- setdiff(threshold_id, group_id)
      if (length(id) >= 1L) {
        id <- paste0(rep(id, each = length(suffix)), suffix)
        summary <- dplyr::select(summary, - dplyr::all_of(id))
      }
      
    }
  }    
  
  # reorder column names 
  
  id <- paste0(rep(names(threshold_groups), each = length(suffix)), suffix)
  
  summary <- dplyr::relocate(
    summary, 
    dplyr::all_of(id), 
    .after = climit_last_year
  )

  list(summary = summary, thresholds = names(threshold_groups))
}
  
  
group_thresholds_engine <- function(x) {
  if (all(is.na(x))) {
    out <- switch(
      class(x), 
      numeric = NA_real_, 
      character = NA_character_
    )
    return(out)
  }
  x <- x[!is.na(x)]
  if (length(x) > 1L) {
    stop("multiple values not allowed")
  }
  x
}


# html reports ----

#' Reports the assessment of individual time series
#'
#' Generates a series of html reports with, for each time series, meta data, 
#' plots of the data with the fitted assessment model, statistical summaries, 
#' and a simple interpretation of the fitted model.
#'
#' @param assessment_obj An assessment object resulting from a call to
#'   run_assessment
#' @param subset An optional vector specifying which timeseries are to be
#'   reported. An expression will be evaluated in the timeSeries component of
#'   assessment_obj; use 'series' to identify individual timeseries.
#' @param output_dir The output directory for the assessment plots (possibly
#'   supplied using 'file.path'). The default is the working directory. The
#'   output directory must already exist.
#' @param output_file An alterntive file name to override the default. This is  
#'   currently only implemented for a single report. If not supplied, the .html
#'   extension will be added. 
#' @param max_report The maximum number of reports that will be generated.
#'   Defaults to 100. Each report is about 1MB in size and takes a few seconds 
#'   to run, so this prevents a ridiculous number of reports being created. 
#'
#' @returns A series of html files with, for each time series, meta data, 
#'   plots of the data with the fitted assessment model, statistical summaries, 
#'   and a simple interpretation of the fitted model.
#'
#'
#' @export
report_assessment <- function(
    assessment_obj, 
    subset = NULL, 
    output_dir = ".",
    output_file = NULL, 
    max_report = 100L) {
  
  # reporting_functions.R
  
  if (!dir.exists(output_dir)) {
    stop(
      "\nThe output directory '", output_dir, "' does not exist.\n", 
      "Create it or check the information supplied to argument 'output_dir'",
      " is correct.",
      call. = FALSE
    )
  }
  
  if (!is.null(output_file) & length(output_file) > 1) {
    stop(
      "\n`output_file` can currently only be a single character string for",
      " renaming a single\nreport.", 
      call. = FALSE
    )
  }

  
  info <- assessment_obj$info
  timeSeries <- assessment_obj$timeSeries   
  
  # set up time series information:
  # - merge with station information
  # - add in additional useful variables 
  # - subset if necessary

  timeSeries <- tibble::rownames_to_column(timeSeries, "series")
  
  timeSeries <- dplyr::left_join(
    timeSeries, 
    assessment_obj$stations, 
    by = "station_code"
  )
  
  timeSeries$group <- ctsm_get_info(
    info$determinand, 
    timeSeries$determinand, 
    "group", 
    info$compartment,
    sep = "_"
  )
  
  timeSeries$distribution <- ctsm_get_info(
    info$determinand, 
    timeSeries$determinand, 
    "distribution"
  )
  
  # if (info$compartment == "water") {
  #   timeSeries$matrix <- "WT"
  # }
  
  timeSeries <- apply_subset(timeSeries, subset, parent.frame())

  if (nrow(timeSeries) == 0L) {
    warning(
      "no timeseries were selected - nothing has been reported"
    )
    return(invisible())
  }

  series_id <- row.names(timeSeries)
  

  
  
  # ensure number of series does not exceed max_report
  
  n_series <- length(series_id)
  
  if (n_series > max_report) {
    stop(
      "\nYou have asked for ", n_series, " reports which exceeds the number ", 
      "allowed.\n", 
      "To continue increase the report limit with the 'max_report' argument.\n", 
      "Be aware that each report will be larger than 1MB.",
      call. = FALSE
    )
  }
    

  # if output_file supplied, ensure there is only one series
  
  if (!is.null(output_file) & n_series > 1) {
    stop(
      "\n`output_file` can currently only be used to rename a single report", 
      " and ", n_series, " reports have\nbeen requested", 
      call. = FALSE
    )
  }
  
  
  # report on each time series
  
  lapply(series_id, function(id) {

    # get file name
    # if not supplied, use id and add country and station name for easier 
    # identification

    if (!is.null(output_file)) {
      
      output_id = output_file
      
    } else {
    
      series <- timeSeries[id, ]
      
      output_id <- sub(
        series$station_code,
        paste(series$station_code, series$country, series$station_name), 
        id,
        fixed=TRUE
      )
      
      # get rid of any slashes that might have crept in 
      
      output_id <- gsub(" / ", " ", output_id, fixed = TRUE)
      output_id <- gsub("/", " ", output_id, fixed = TRUE)
      
      output_id <- gsub(" \ ", " ", output_id, fixed = TRUE)
      output_id <- gsub("\\", " ", output_id, fixed = TRUE)

      # and any % e.g. %DNATAIL!
      
      output_id <- gsub("%", "", output_id, fixed = TRUE)
    }
          
    package_dir = system.file(package = "harsat")
    template_dir = file.path(package_dir, "markdown")
    report_file <- file.path(template_dir, "report_assessment.Rmd")

    rmarkdown::render(
      report_file, 
      params = list(
        assessment_object = assessment_obj, 
        series = id
      ),
      output_file = output_id, 
      output_dir = output_dir
    )
  })
    
  invisible() 
}



# OHAT ----

write_hat <- function(assessments) {
  
  output_path <- file.path("output", "example_OSPAR", "hat")
  
  assessment_id <- names(assessments)
  
  wk <- lapply(assessment_id, function(id) {

    assessment <- assessments[[id]]

    out <- read.csv(
      file.path("output", "example_OSPAR", paste0(id, "_summary.csv")), 
      na.strings = "", 
      fileEncoding = "UTF-8-BOM"
    )
    
    var_id <- c(
      "series", 
      assessment$info$region$names, 
      "country",           
      "station_code", "station_name", "station_longname", 
      "station_latitude", "station_longitude", 
      "species", 
      "determinand_group", "determinand", 
      "filtration", "matrix", "basis", "unit", 
      "sex", "method_analysis", "subseries", 
      "shape", "colour", "shape_env", "colour_env", 
      "shape_health", "colour_health"
    )
    
    out <- dplyr::select(out, dplyr::any_of(var_id))
    
    out <- dplyr::rename(
      out, 
      # region = "OSPAR_region",
      # subregion = "OSPAR_subregion",
      latitude = "station_latitude",
      longitude = "station_longitude",
      measurement_type = "determinand_group"
    )
    
    out <- dplyr::mutate(
      out,
      measurement_type = factor(
        measurement_type, 
        levels = c(
          "Metals", "Organotins", 
          "PAH parent compounds", "PAH alkylated compounds", "PAH metabolites", 
          "Polybrominated diphenyl ethers", "Organobromines (other)", 
          "Organofluorines", 
          "Polychlorinated biphenyls", "Dioxins", "Organochlorines (other)",
          "Pesticides", 
          "Imposex", "Biological effects (other)"
        )
      ), 
      # determinand = factor(determinand, levels = determinand_order)
      determinand = factor(determinand)
    )
    
    if (any(is.na(out$determinand))) {
      stop("missing some determinand - investigate")
    }
    
    out$measurement <- ctsm_get_info(
      assessment$info$determinand, out$determinand, "common_name"
    )
    
    if (id %in% "biota") {
      out <- dplyr::mutate(
        out,
        species_latin = species,
        family = ctsm_get_info(
          assessment$info$species, species_latin, "species_group"
        ), 
        species = ctsm_get_info(
          assessment$info$species, species_latin, "common_name"
        ),
        family = dplyr::recode(
          family, 
          Crustacean = "Shellfish", 
          Bivalve = "Shellfish",
          Gastropod = "Shellfish"
        )
      )
    }
    
    if (id %in% "water") {
      out <- dplyr::mutate(
        out,
        menu1_title = "Filtration",
        menu1_entry = stringr::str_to_sentence(filtration),
        menu2_title = NA_character_,
        menu2_entry = NA_character_
      )
    }
    
    if (id %in% "sediment") {
      out <- dplyr::mutate(
        out,
        menu1_title = "Sediment fraction",
        menu1_entry = dplyr::recode(
          matrix, 
          SED20 = "<20 micron",
          SED63 = "<63 micron",
          SEDTOT = "Total sediment"
        ),
        menu2_title = NA_character_,
        menu2_entry = NA_character_
      )
    }
    
    if (id %in% "biota") {
      
      out <- dplyr::mutate(
        out,
        menu1_title = dplyr::case_when(
          measurement_type %in% "PAH metabolites" ~ "Chemical analysis",
          determinand %in% "EROD" ~ "Tissue",
          !(measurement_type %in% c("Imposex", "Biological effects (other)")) ~ 
            "Tissue",
          TRUE ~ NA_character_
        ), 
        menu1_entry = dplyr::case_when(
          measurement_type %in% "PAH metabolites" ~ method_analysis,
          determinand %in% "EROD" ~ matrix,
          !(measurement_type %in% c("Imposex", "Biological effects (other)")) ~ 
            matrix,
          TRUE ~ NA_character_
        ), 
        menu2_title = dplyr::case_when(
          determinand %in% "EROD" ~ "Sex",
          !(measurement_type %in% 
              c("PAH metabolites", "Imposex", "Biological effects (other)")) ~ 
            "Animal grouping",
          TRUE ~ NA_character_
        ), 
        menu2_entry = dplyr::case_when(
          determinand %in% "EROD" ~ sex,
          !(measurement_type %in% 
              c("PAH metabolites", "Imposex", "Biological effects (other)")) ~ 
            subseries,
          TRUE ~ NA_character_
        ), 
        .matrix = ctsm_get_info(assessment$info$matrix, matrix, "name"), 
        .matrix = stringr::str_to_sentence(.matrix),  
        .matrix = dplyr::recode(
          .matrix,
          "Erythrocytes (red blood cells in vertebrates)" = "Red blood cells",
          "Egg homogenate of yolk and albumin" = "Egg yolk & albumin", 
          "Liver s9 fraction" = "Liver S9 fraction"
        ), 
        menu1_entry = dplyr::if_else(
          menu1_title %in% "Tissue", 
          .matrix, 
          menu1_entry
        ),
        menu2_entry = dplyr::case_when(
          menu2_title %in% "Sex" ~ dplyr::recode(
            menu2_entry, "F" = "Female", "M" = "Male"
          ),
          menu2_title %in% "Animal grouping" & is.na(menu2_entry) ~ "All",
          TRUE ~ menu2_entry
        ),
        menu2_entry = sub("_", " ", menu2_entry, fixed = TRUE),
        .matrix = NULL
      )
      
      
      # only keep menu2_title and menu2_entry values for mammal_group
      # if there is a mammal time series for that determinand
      
      det_id <- out |> 
        dplyr::filter(
          menu2_title %in% "Animal grouping",
          ! menu2_entry %in% "All"
        ) |> 
        dplyr::pull(determinand) |> 
        as.character() |> 
        unique()
      
      det_id <- c(det_id, "EROD")
      
      out <- dplyr::mutate(
        out, 
        .change = determinand %in% det_id,
        menu2_title = dplyr::if_else(.change, menu2_title, NA_character_),
        menu2_entry = dplyr::if_else(.change, menu2_entry, NA_character_),
        .change = NULL
      )
      
    }
    
    
    col_id = c(
      "series", 
      assessment$info$region$names, 
      "country", 
      "station_code", "station_name", "station_longname", "latitude", "longitude", 
      "family", "species_latin", "species",
      "measurement_type", "determinand", "measurement", 
      "filtration", "matrix", "basis", "unit", 
      "sex", "method_analysis", "subseries", 
      "menu1_title", "menu1_entry", "menu2_title", "menu2_entry", 
      "shape", "colour", "shape_env", "colour_env", "shape_health", "colour_health"
    )
    
    col_id <- intersect(col_id, names(out))
    
    out <- out[col_id]
    
    
    # order shape and colour so that less important time series are plotted first 
    # (to avoid masking)
    
    if (id %in% c("biota", "water", "sediment")) {
      
      out <- out |> 
        dplyr::mutate(
          .shape = dplyr::recode(
            shape, 
            upward_triangle = "up", 
            downward_triangle = "down",
            large_filled_circle = "flat",
            small_filled_circle = "mean",
            small_open_circle = "other"
          )
        ) |> 
        tidyr::unite(".ord", .shape, colour, remove = FALSE) |> 
        dplyr::mutate(
          .shape = NULL,
          .ord = factor(
            .ord,
            levels = c(
              paste("other", c("black", "blue", "orange", "green", "red"), sep = "_"),
              "mean_black", "flat_black", 
              paste("mean", c("blue", "orange", "green", "red"), sep = "_"),
              paste("flat", c("blue", "orange", "green", "red"), sep = "_"),
              paste("down", c("black", "blue", "orange", "green", "red"), sep = "_"),
              paste("up", c("black", "blue", "orange", "green", "red"), sep = "_")
            ), 
            ordered = TRUE
          )
        )
      
      # out <- dplyr::arrange(
      #   out, measurement_type, determinand, region, subregion, .ord, station_code)
      
      out <- dplyr::arrange(
        out, measurement_type, determinand, .ord, station_code)
      
      out <- dplyr::mutate(out, .ord = NULL)
      
      out$order <- 1:nrow(out)
      
    }
    
    
    if (id %in% "biota_env") {
      
      out <- out |> 
        dplyr::mutate(
          .shape = dplyr::recode(
            shape_env, 
            upward_triangle = "up", 
            downward_triangle = "down",
            large_filled_circle = "flat",
            small_filled_circle = "mean",
            small_open_circle = "other"
          )
        ) |> 
        tidyr::unite(".ord", .shape, colour_env, remove = FALSE) |> 
        dplyr::mutate(
          .shape = NULL,
          .ord = factor(
            .ord,
            levels = c(
              paste("other", c("black", "blue", "orange", "green", "red"), sep = "_"),
              "mean_black", "flat_black", 
              paste("mean", c("blue", "orange", "green", "red"), sep = "_"),
              paste("flat", c("blue", "orange", "green", "red"), sep = "_"),
              paste("down", c("black", "blue", "orange", "green", "red"), sep = "_"),
              paste("up", c("black", "blue", "orange", "green", "red"), sep = "_")
            ), 
            ordered = TRUE
          )
        )
      
      # out <- dplyr::arrange(
      #   out, measurement_type, determinand, region, subregion, .ord, station_code)
      
      out <- dplyr::arrange(
        out, measurement_type, determinand, .ord, station_code)

      out <- dplyr::mutate(out, .ord = NULL)
      
      out$order_env <- 1:nrow(out)
      
      
      out <- out |> 
        dplyr::mutate(
          .shape = dplyr::recode(
            shape_health, 
            upward_triangle = "up", 
            downward_triangle = "down",
            large_filled_circle = "flat",
            small_filled_circle = "mean",
            small_open_circle = "other"
          )
        ) |> 
        tidyr::unite(".ord", .shape, colour_health, remove = FALSE) |> 
        dplyr::mutate(
          .shape = NULL,
          .ord = factor(
            .ord,
            levels = c(
              paste("other", c("black", "blue", "orange", "green", "red"), sep = "_"),
              "mean_black", "flat_black", 
              paste("mean", c("blue", "orange", "green", "red"), sep = "_"),
              paste("flat", c("blue", "orange", "green", "red"), sep = "_"),
              paste("down", c("black", "blue", "orange", "green", "red"), sep = "_"),
              paste("up", c("black", "blue", "orange", "green", "red"), sep = "_")
            ), 
            ordered = TRUE
          )
        )
      
      # out <- dplyr::arrange(
      #   out, measurement_type, determinand, region, subregion, .ord, station_code)
      
      out <- dplyr::arrange(
        out, measurement_type, determinand, .ord, station_code)
      
      out <- dplyr::mutate(out, .ord = NULL)
      
      out$order_health <- 1:nrow(out)
      
      # out <- dplyr::arrange(
      #   out, measurement_type, determinand, region, subregion, order_env, station_code)
      
      out <- dplyr::arrange(
        out, measurement_type, determinand, order_env, station_code)

      out <- dplyr::relocate(out, order_env, .after = colour_env)
    }
    
    
    outfile <- file.path(output_path, paste0(id, "_summary.csv"))
    
    readr::write_excel_csv(out, outfile, na = "")
    
    out
  })
  
  names(wk) <- stringr::str_to_sentence(assessment_id)  

  
  # OHAT schema
  
  # measurement_type, menu1_title, menu2_title, health, order
  
  # identify health determinands with no AC
  
  # wk_id <- wk$Biota |>
  #   group_by(determinand) |>
  #   summarise(all_black = all(colour_health %in% "black"), .groups = "drop_last") |>
  #   dplyr::filter(!all_black) |>
  #   pull(determinand) |>
  #   as.character()
   
  # get animal grouping in 'correct' order
  
  if ("biota" %in% assessment_id) {
    wk_ag <- unique(wk$Biota$menu2_entry) |> na.omit() |> c()
    wk_ag <- sort(wk_ag)
    wk_ag <- c(setdiff(wk_ag, "All"), "All")
    wk$Biota <- dplyr::mutate(
      wk$Biota, 
      menu2_entry = factor(menu2_entry, levels = wk_ag)
    )
  }
    
      
  wk2 <- wk |> 
    lapply(function(ls) {
      if ("colour_health" %in% names(ls)) {
        ls <- dplyr::mutate(
          ls, 
          health = if_else(determinand %in% wk_id, "TRUE", NA_character_)
        )
      }
      id <- c("measurement_type", "measurement", "menu1_title", "menu2_title", "health")
      ls |> 
        dplyr::distinct(dplyr::across(any_of(id))) |> 
        dplyr::mutate_if(is.factor, as.character) |> 
        dplyr::mutate(order = 1:dplyr::n())
    }) |> 
    dplyr::bind_rows(.id = "compartment")
  
  readr::write_excel_csv(
    wk2, 
    file.path(output_path, "schema_1.csv"), 
    na = ""
  )
  
  # additionally determinand, menu1_entry, menu2_entry
  
  wk2 <- wk |> 
    lapply(function(ls) {
      if ("colour_health" %in% names(ls)) {
        ls <- dplyr::mutate(ls, health = if_else(determinand %in% wk_id, "TRUE", NA_character_))
      }
      id <- c(
        "measurement_type", "determinand", "measurement", "menu1_title", "menu1_entry",
        "menu2_title", "menu2_entry", "health")
      ls |> 
        dplyr::distinct(across(any_of(id))) |> 
        dplyr::arrange(
          measurement_type, determinand, 
          menu1_title, menu1_entry, menu2_title, menu2_entry
        ) |> 
        dplyr::mutate(
          determinand = NULL, 
          order = 1:dplyr::n()) |> 
        dplyr::mutate_if(is.factor, as.character)
    }) |> 
    dplyr::bind_rows(.id = "compartment")
  
  readr::write_excel_csv(
    wk2, 
    file.path(output_path, "schema_2.csv"), 
    na = ""
  )
  
  
  # biota gives species groups
  # redefining family as a factor allows us to ge the ordering as desired
  # - trouble with using foracts::fct_relevel is that it issues a warning if
  # a level doesn't exist in the data

  if ("biota" %in% assessment_id) {
    
    wk2 <- wk$Biota |> 
      dplyr::distinct(measurement_type, family) |> 
      dplyr::mutate(
        family = factor(
          as.character(family), 
          c("Shellfish", "Fish", "Bird", "Mammal")
        )
      )  |> 
      dplyr::arrange(measurement_type, family) |> 
      dplyr::mutate_if(is.factor, as.character) |> 
      dplyr::mutate(order = 1:dplyr::n())
    
    readr::write_excel_csv(
      wk2, 
      file.path(output_path, "schema_3.csv"), 
      na = ""
    )

  }    
  
  
  # legends and categories
  
  # determinand_order is used to order determinands how you want them to appear
  # the code below is a hack to get things to work
  
  determinand_order <- lapply(assessments, function(x) {
    unique(x$data$determinand)
  })
  
  determinand_order <- unlist(determinand_order)
  determinand_order <- sort(unique(determinand_order))

  determinand_groups <- list(
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
  
    
  
  wk <- ctsm_OHAT_legends(
    assessments = assessments,
    determinandGroups = determinand_groups,
    determinands = determinand_order,
    symbology = list(
      biota = list(
        below = c(
          "BAC" = "blue", "NRC" = "blue", "EAC" = "green", "FEQG" = "green",
          "LRC" = "green", "QSsp" = "green"
        ),
        above = c(
          "BAC" = "orange", "NRC" = "orange", "EAC" = "red", "FEQG" = "red",
          "LRC" = "red", "QSsp" = "red"
        ),
        none = "black"
      ),
      sediment = list(
        below = c(
          "BAC" = "blue", "ERL" = "green", "EAC" = "green", "EQS" = "green", 
          "FEQG" = "green"
        ),
        above = c(
          "BAC" = "orange", "ERL" = "red", "EAC" = "red", "EQS" = "red", 
          "FEQG" = "red"
        ),
        none = "black"
      ),
      water = list(
        below = c("EQS" = "green"), 
        above = c("EQS" = "red"), 
        none = "black"
      )
    ),    
    regionalGroups = list(
      biota = c(
        "Metals", "PAH parent compounds", "PAH metabolites", "Polychlorinated biphenyls", 
        "Polybrominated diphenyl ethers", "Imposex"
      ),
      sediment = c(
        "Metals", "PAH parent compounds", "Polychlorinated biphenyls", 
        "Polybrominated diphenyl ethers", "Organotins"
      )
    ),
    distanceGroups = list(
      biota = c("Metals", "PAH parent compounds"),
      sediment = c("Metals", "PAH parent compounds")
    ),
    path = output_path
  )
  
  
  # wk_health <- ctsm_OHAT_legends(
  #   assessments = list(biota = biota_assessment),
  #   determinandGroups = determinand_groups,
  #   determinands = determinand_order, 
  #   symbology = list(
  #     biota = list(
  #       below = c("MPC" = "green", "QShh" = "green"),
  #       above = c("MPC" = "red", "QShh" = "red"),
  #       none = "black"
  #     )
  #   ),
  #   regionalGroups = list(
  #     Biota = c(
  #       "Metals", "PAH parent compounds", "PAH metabolites", "Polychlorinated biphenyls", 
  #       "Polybrominated diphenyl ethers", "Imposex"
  #     ),
  #     Sediment = c(
  #       "Metals", "PAH parent compounds", "Polychlorinated biphenyls", 
  #       "Polybrominated diphenyl ethers", "Organotins"
  #     )
  #   ),
  #   distanceGroups = list(
  #     Biota = c("Metals", "PAH parent compounds"),
  #     Sediment = c("Metals", "PAH parent compounds")
  #   ),
  #   path = output_path
  # )
  
  
  # ad-hoc corrections
  
  wk$legends <- wk$legends |> 
    dplyr::mutate(
      .biota_HG = Compartment == "Biota" & Determinand_code == "HG",
      Label = dplyr::case_when(
        .biota_HG & Label == "below BAC"  ~ "below BAC or NRC",
        .biota_HG & Label == "below QSsp" ~ "below QSsp or LRC",
        .biota_HG & Label == "above QSsp" ~ "above QSsp or LRC",
        TRUE                              ~ Label
      ),
      Tooltip = dplyr::if_else(
        .biota_HG & Label == "below BAC or NRC", 
        "background or no-risk concentrations", 
        Tooltip
      )
    ) |>
    dplyr::select(- .biota_HG)
  
  
  # wk$legends <- dplyr::bind_rows(
  #   list("Environmental" = wk$legends, "Human health" = wk_health$legends),
  #   .id = "Threshold"
  # )
  
  # wk$legends <- wk$legends |> 
  #   dplyr::mutate(
  #     Determinand_code = factor(Determinand_code) |> forcats::fct_inorder()
  #   ) |> 
  #   dplyr::arrange(Compartment, Determinand_code, Threshold) 
  
  wk$legends <- wk$legends |> 
    dplyr::mutate(
      Determinand_code = factor(Determinand_code) |> forcats::fct_inorder()
    ) |> 
    dplyr::arrange(Compartment, Determinand_code) 

  wk$legends <- dplyr::select(wk$legends, Compartment, everything())
  
  wk$legends$order <- 1:nrow(wk$legends)
  
  names(wk$legends) <- tolower(names(wk$legends))
  
  readr::write_excel_csv(
    wk$legends, 
    file.path(output_path, "legends.csv"), 
    na = ""
  )
  
  names(wk$help) <- tolower(names(wk$help))
  
  readr::write_excel_csv(
    wk$help, 
    file.path(output_path, "help_and_information.csv"), 
    na = ""
  )
  
  invisible()
  
}



#' @export
ctsm_OHAT_legends <- function(
  assessments, determinandGroups, determinands, symbology, 
  regionalGroups = NULL, distanceGroups = NULL, path) {

  # silence non-standard evaluation warnings
  info <- NULL

  out <- sapply(names(assessments), simplify = FALSE, USE.NAMES = TRUE, FUN = function(media) {

    assessment <- assessments[[media]]
    classColour <- symbology[[media]]
    regionalGroups <- regionalGroups[[media]]
    distanceGroups <- distanceGroups[[media]]
    
    legends <- ctsm.web.AC(assessment, classColour)
        
    determinands <- intersect(determinands, rownames(legends))   # need this to keep ordering right
    legends <- legends[determinands, , drop = FALSE]   

    compartment <- assessment$info$compartment
    group <- ctsm_get_info(
      assessment$info$determinand, determinands, "group", compartment, sep = "_"
    )
    web_group <- factor(
      group, 
      levels = determinandGroups$levels, 
      labels = determinandGroups$labels
    )
    web_group <- web_group[, drop = TRUE]
    
    goodStatus <- ctsm_get_info(assessment$info$determinand, determinands, "good_status")
    goodStatus <- as.character(goodStatus)

    legendName <- apply(legends, 1, function(i) paste(colnames(legends)[i], collapse = " "))
    legendName <- paste(compartment, group, goodStatus, legendName, sep = " ")

    legends <- data.frame(legends, legendName, group, web_group, goodStatus, stringsAsFactors = FALSE)

    ctsm_OHAT_add_legends(legends, classColour, regionalGroups, distanceGroups, assessment$info)
  })
  
  legends <- lapply(out, "[[", "legends") |> 
    dplyr::bind_rows(.id = "Compartment")
  
  help <- lapply(out, "[[", "help") |> 
    dplyr::bind_rows(.id = "Compartment")
  
  list(legends = legends, help = help)
}  
  



ctsm_OHAT_add_legends <- function(legends, classColour, regionalGroups, distanceGroups, info) {
  
  # get shape for good status = low and good status = high
  # for the latter, just need to sway the interpretation of a downward and upward trend

  recent.period <- paste("last", info$recent.trend, "years")

  standard_shape <- list(
    list(
      Legend = "trend", 
      Colour = classColour[["none"]],      
      Shape = "downward_triangle", 
      Label = "downward trend", 
      Tooltip = paste("concentrations decreased in", recent.period)
    ),
    list(
      Legend = "trend", 
      Colour = classColour[["none"]],      
      Shape = "upward_triangle", 
      Label = "upward trend", 
      Tooltip = paste("concentrations increased in", recent.period)
    ),
    list(
      Legend = "trend", 
      Colour = classColour[["none"]],      
      Shape = "large_filled_circle", 
      Label = "no trend", 
      Tooltip = paste("concentrations stable over", recent.period)
    ),
    list(
      Legend = "trend", 
      Colour = classColour[["none"]],      
      Shape = "small_filled_circle", 
      Label = "status assessment only", 
      Tooltip = "not enough data to assess trends"
    ),
    list(
      Legend = "trend", 
      Colour = classColour[["none"]],      
      Shape = "small_open_circle", 
      Label = "informal status assessment", 
      Tooltip = "only 1-2 years of data"
    )
  )

  standard_shape <- standard_shape |> 
    dplyr::bind_rows() |> 
    as.data.frame()
  
  standard_shape <- list(low = standard_shape, high = standard_shape)
  
  standard_shape$high$Tooltip[1:2] <- paste("levels got", c("worse", "better"), "in", recent.period)

  AC_explanation <- list(
    low = c("below BAC" = "near background concentrations",
            "above BAC" = "above background concentrations",
            "below ERL" = "few adverse effects on marine life",
            "above ERL" = "adverse effects on marine life",
            "below EAC" = "few adverse effects on marine life",
            "above EAC" = "adverse effects on marine life",
            "below FEQG" = "few adverse effects on marine life",
            "above FEQG" = "adverse effects on marine life",
            "below EQS" = "few adverse effects on marine life",
            "above EQS" = "adverse effects on marine life",
            "below EQS.OSPAR" = "few adverse effects on marine life",
            "above EQS.OSPAR" = "adverse effects on marine life",
            "below ERM" = "some adverse effects on marine life",
            "above ERM" = "many adverse effects on marine life",
            "below MPC" = "safe to eat",
            "above MPC" = "do not eat me",
            "below HQS" = "safe to eat",
            "above HQS" = "do not eat me",
            "no assessment criteria" = "assessment criteria under development"),
    high = c("above BAC" = "near background levels",
             "below BAC" = "below background levels",
             "above EAC" = "few adverse effects on marine life",
             "below EAC" = "adverse effects on marine life",
             "no assessment criteria" = "assessment criteria under development")
  )

  # split up legends data frame into AC identification and legendName

  legendName <- legends$legendName
  group <- as.character(legends$group)
  goodStatus <- legends$goodStatus
  regional <- legends$web_group %in% regionalGroups
  distance <- legends$web_group %in% distanceGroups
  legends <- legends[!(names(legends) %in% c("legendName", "group", "web_group", "goodStatus"))]
  

  legend_info <- lapply(1:nrow(legends), function(i) {

    # get names of each AC and the appropriate colour
    # add on an extra row for 'above' category - this needs to be tidied up, because 'above' is a legacy
    #   when only contaminants were being considered
    # deal with the cases where there are no AC at all
    # deal with the case where there is a BAC and no EAC for some determinands and an EAC for others

    ACid <- colnames(legends)[unlist(legends[i,])]

    AC.none <- "none" %in% ACid           # tests for stations with no AC
    AC.onlyBAC <- "BAC_only" %in% ACid     # tests for stations with BAC and no EAC
    ACnames <- subset(ACid, !(ACid %in% c("none", "BAC_only")))      # gets all AC


    if (length(ACnames) > 0) {                   # indicates some AC available
      AClast <- tail(ACnames, 1)
      ACsymbol <- classColour[["below"]][ACnames]
      ACnames <- c(ACnames, AClast)
      ACsymbol <- c(ACsymbol, classColour[["above"]][AClast])
      AClabel <- paste(c(rep("below", length(ACnames) - 1), "above"), ACnames)
    } else {
      ACsymbol <- AClabel <- character(0)
    }

    if (AC.onlyBAC & !("above BAC" %in% AClabel)) {
      ACsymbol <- append(ACsymbol, classColour[["above"]]["BAC"], 1)
      AClabel <- append(AClabel, "above BAC", 1)
      ACnames <- append(ACnames, "BAC", 1)
    }

    if (AC.none) {
      ACsymbol <- c(ACsymbol, classColour[["none"]])
      AClabel <- c(AClabel, "no assessment criteria")
      ACnames <- c(ACnames, "none")
    }


    statusID <- goodStatus[i]

    if (statusID == "high") {
      belowID <- grepl("below", AClabel)
      aboveID <- grepl("above", AClabel)
      AClabel[belowID] <- sub("below", "above", AClabel[belowID])
      AClabel[aboveID] <- sub("above", "below", AClabel[aboveID])
    }

    AC_explanation <- AC_explanation[[statusID]][AClabel]
    
    out <- data.frame(
      Legend = "status", 
      Colour = ACsymbol,
      Shape = "large_filled_circle",
      Label = AClabel, 
      Tooltip = AC_explanation
    )
    
    out <- dplyr::mutate_if(out, is.factor, as.character)
    
    out <- dplyr::bind_rows(out, standard_shape[[statusID]])
    
    out
  })
    

  # add in regional assessment (if present), distance to BC (if present), methods, AC and FAQ help files

  help_info <- lapply(1:nrow(legends), function(i) {
    
    if (regional[i]) {
      info_regional <- list(
        Legend = "info", 
        File = tolower(paste0("regional_assessment_", info$compartment, "_", group[i], ".html")),
        Label = "Regional assessment",
        Order = 1
      )
    } else {
      info_regional <- NULL
    }
    
    if (distance[i]) {
      info_distance <- list(
        Legend = "info", 
        File = tolower(paste0("distance_LC_", info$compartment, "_", group[i], ".html")),
        Label = "Distance to background",
        Order = 1 + regional[i]
      )
    } else {
      info_distance <- NULL
    }

    not_contaminant <- group[i] %in% c("Metabolites", "Effects", "Imposex")
    group_txt <- if (not_contaminant) group[i] else "contaminants"

    help_info <- list(
      list(
        Legend = "help",
        File = tolower(paste0("help_methods_", info$compartment, "_", group_txt, ".html")), 
        Label = "Assessment methodology", 
        Order = 1
      ), 
      list(
        Legend = "help",
        File = tolower(paste0("help_ac_", info$compartment, "_", group_txt, ".html")), 
        Label = "Assessment criteria", 
        Order = 2
      ),
      list(
        Legend = "help",
        File = "help_faq.html",
        Label = "Frequently asked questions", 
        Order = 3
      )
    )
    
    out <- dplyr::bind_rows(info_regional, info_distance, help_info)
    
    out <- as.data.frame(out)
    
    out
  })
    
  names(legend_info) <- names(help_info) <- row.names(legends)
  
  legend_info <- dplyr::bind_rows(legend_info, .id = "Determinand_code")
  help_info <- dplyr::bind_rows(help_info, .id = "Determinand_code")
  
  list(legends = legend_info, help = help_info)
}
