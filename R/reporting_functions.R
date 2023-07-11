# mapping functions ----

ctsm.Mercator <- function(x) atanh(sin(x * pi / 180))

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

#' @export
ctsm_subset_assessment <- function(assessment_obj, subset) {
  
  # reporting_functions.R
  # subsets an assessment object by filtering on the timeSeries component
  
  timeSeries <- assessment_obj$timeSeries
  
  ok <- eval(substitute(subset), timeSeries, parent.frame())
  timeSeries <- timeSeries[ok, ]
  
  assessment_obj$timeSeries <- timeSeries

  
  # update other components to be consistent
  
  id <- row.names(timeSeries)
  
  assessment_obj$assessment <- assessment_obj$assessment[id]
  
  ok <- assessment_obj$data$seriesID %in% id
  assessment_obj$data <- assessment_obj$data[ok, ]
  
  ok <- row.names(assessment_obj$stations) %in% timeSeries$station 
  assessment_obj$stations <- assessment_obj$stations[ok, ]
  
  
  # drop redundant factor levels
  
  id <- c("timeSeries", "data", "stations")
  assessment_obj[id] <- lapply(assessment_obj[id], droplevels)
  
  assessment_obj
}  


#' @export
ctsm_summary_overview <- function(
    assessment, timeSeries, info, classColour, fullSummary = FALSE) {
  
  # reporting_functions.R
  
  # gets shape and colour for each time series
  
  # first get list of assessment summaries 
  # summary structures differ between detGroups, so need to be careful
  
  # order assessment so that it is compatible with timeSeries - 
  # need to resolve this
  
  assessment <- assessment[row.names(timeSeries)]
  
  summaryList <- sapply(assessment, "[[", "summary", simplify = FALSE)
  
  if (any(is.null(summaryList)) | (length(summaryList) != nrow(timeSeries))) {
    stop("coding error - contact HARSAT development team")
  }
    
  
  # get combined summary names across detGroups
  
  summaryNames <- unique(do.call("c", lapply(summaryList, names)))
  
  # create enlarged structure to hold summaries
  
  summaryObject <- do.call(
    "data.frame", 
    sapply(summaryNames, USE.NAMES = TRUE, simplify = FALSE, FUN = function (i) NA)
  )
  
  # now write each summary to the enlarged structure
  
  out <- do.call("rbind", sapply(summaryList, USE.NAMES = TRUE, simplify = FALSE, FUN = function(x)
  {
    if (is.null(x)) return(summaryObject)
    out <- summaryObject
    out[names(x)] <- x
    out
  }))
  
  
  # get shape and colour of plotting symbols

  out <- ctsm_symbology_OSPAR(out, info, timeSeries, classColour)

  if (fullSummary) out else out[c("shape", "colour")]
}


#' @export
ctsm_symbology_OSPAR <- function(summary, info, timeSeries, classColour, alpha = 0.05) {
  
  # reporting_functions.R
  
  summary$shape <- with(summary, {
    
    shape <- character(nrow(summary))
    
    # trend symbols: a trend is estimated if pltrend is present - default shape is a large filled circle
    
    trendFit <- !is.na(pltrend)
    shape[trendFit] <- "large_filled_circle"
    
    # show a significant trend based on prtrend 
    # NB prtrend might not exist even if pltrend does, because there are too few years of data in the
    # recent window
    
    # isImposex <- timeSeries$detGroup %in% "Imposex"
    
    isTrend <- !is.na(prtrend) & prtrend < alpha
    upTrend <- isTrend & rtrend > 0
    downTrend <- isTrend & rtrend < 0
    
    shape[downTrend] <- "downward_triangle"
    shape[upTrend] <- "upward_triangle"
    
    
    # status based on upper confidence limit in last year, but for imposex, clLY available only with
    # 1 or 2 years if individual data, so use nyfit instead
    
    # statusFit <- !trendFit & !is.na(clLY) 
    statusFit <- !trendFit & nyfit >= 3
    shape[statusFit] <- "small_filled_circle"
    
    shape[!trendFit & !statusFit] <- "small_open_circle"
    
    shape
  })


  # colour based on 
  # - upper confidence limit (nyfit > 2 or, for VDS, individual measurements) 
  # - meanly (nyfit <= 2)
  
  # get the names of the variables which contain the difference between the 
  # meanly and the AC
  
  if (!is.null(info$AC)) {
    
    ACdiff <- paste(info$AC, "diff", sep = "")
    
    
    # when goodStatus is indicated by low concentrations, negative ACdiff is good
    # since ACdiff = clLY - AC < 0
    # when indicated by high concentrations, positive ACdiff is good
    # to make colour calculation 'simple' change sign on ACdiff when goodStatus == high
    
    goodStatus <- ctsm_get_info(info$determinand, timeSeries$determinand, "good_status")
    
    wk <- summary[ACdiff]
    wk[] <- lapply(wk, "*", ifelse(goodStatus == "low", 1, -1))
    
    summary$colour <- apply(wk, 1, function(x) {
      
      if (all(is.na(x))) return(classColour$none)
      
      AC <- info$AC[!is.na(x)]
      x <- x[!is.na(x)]
      
      if (any(x < 0)) classColour$below[AC[which.max(x < 0)]]
      else classColour$above[AC[length(x)]]
    })  
    
    # need to adjust for null summaries
    
    summary <- within(summary, colour[is.na(shape)] <- NA)
    
  } else {
    
    summary$colour <- classColour$none
    
  }
  
  
  # adjust shape and colour for nonparametric test if nyfit <= 2 and nyall > 2
  
  # get the names of the variables which contain the result of the non-parametric test for each AC
  
  if (!is.null(info$AC)) {
    
    ACbelow <- paste(info$AC, "below", sep = "")
    
    wk <- summary[ACbelow]
    wk[] <- lapply(wk, function(x) {
      
      ok <- !is.na(x) & goodStatus == "high"
      if (any(ok))
        x[ok] <- ifelse(x[ok] == "below", "above", "below")
      x
    })
    
    wk <- apply(wk, 1, function(x) {
      
      if (all(is.na(x))) return(NA)
      
      AC <- info$AC[!is.na(x)]
      x <- x[!is.na(x)]
      
      if (any(x == "below")) classColour$below[AC[which.max(x == "below")]]
      else classColour$above[AC[length(x)]]
    })  
    
    id <- with(summary, nyfit <= 2 & nyall > 2 & !is.na(wk))
    summary$colour[id] <- wk[id]
    summary$shape[id] <- "small_filled_circle"
    
  }
  
  summary
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
  
  # identity all AC that are relevant to the overall assessment
  # more AC might be included in the assessment for each timeseries - this
  # is a legacy issue that needs to be resolved - has arisen in looking
  # at both environmental and health criteria
  
  AC_id <- assessment_ob$info$AC
  stopifnot(AC_id %in% names(classification[["below"]]))

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
#' * the status assessments (if there any thresholds
#' * (optionally) a symbology summarising the trend (shape) and status (colour)
#' of each time series
#'
#' @param assessment_obj An assessment object resulting from a call to
#'   run_assessment.
#' @param output_file The name of the output csv file. If using NULL, the file
#'   will be called biota_summary.csv, sediment_summary.csv or water_summary.csv
#'   as appropriate. By default the file will be written to the working
#'   directory. If a file name is provided, a path to the output file can also
#'   be provided (e.g. using `file.path`). The output_dir option can also be
#'   used to specify the output file directory.
#' @param output_dir The output directory for `output_file`. The default is the
#'   working directory. Any file path provided in `output_file`, will be
#'   appended to `output_dir`. The resulting output directory must already
#'   exist.
#' @param export Logical. TRUE (the default) writes the summary table to a csv
#'   file. FALSE returns the summary table as an R object (and does not write to
#'   a csv file).
#' @param determinandGroups
#' @param classColour
#' @param collapse_AC
#'
#' @returns
#'
#' @export
write_summary_table <- function(
  assessment_obj, output_file = NULL, output_dir = ".", export = TRUE, 
  determinandGroups = NULL, classColour = NULL, collapse_AC = NULL) {

  # reporting_functions.R
  
  library(dplyr)
  
  assessment <- assessment_obj$assessment
  timeSeries <- assessment_obj$timeSeries
  info <- assessment_obj$info
  
  
  # output information
  # check valid extension and path
  # merge output_file and output_dir to give final output destination
  
  if (export) {

    # get default output_file 
    
    if (is.null(output_file)) {
      output_file <- paste0(info$compartment, "_summary.csv")
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

  }
  
    
  # assessment criteria must have an appropriate class colour
    
  is_AC <- !is.null(info$AC)
  
  if (is_AC) {
    AC <- info$AC
    stopifnot(
      !is.null(classColour),
      AC %in% names(classColour$below), 
      AC %in% names(classColour$above),
      "none" %in% names(classColour)
    )
  }


  # assessment criteria in collapse_AC must be present in the data
  # can turn some of these errors into warnings by creating extra variables 
  # later on - have raised an issue 
  
  if (!is.null(collapse_AC)) {
    
    if (is.null(names(collapse_AC))) {
      stop(
        "collapse_AC must be a names list of valid assessment criterion", 
        call. = FALSE
      )
    }
    
    if (is_AC) {

      id <- unlist(collapse_AC)
      if (any(duplicated(id))) {
        stop(
          "cannot specify the same assessment criterion more than once in ", 
          "collapse_AC", 
          call. = FALSE
        )
      }
      
      ok <- id %in% AC
      if (!all(ok)) {
        id <- id[!ok]
        stop(
          "these assessment criteria are specified in collapse_AC but were not ",
          "used in the\n", 
          "in the assessment: ",
          paste(id, collapse = ", "), 
          call. = FALSE
        )
      }
      
      other_AC <- !AC %in% unlist(collapse_AC)
      if (any(other_AC)) {
        id <- AC[other_AC]
        if (any(id %in% names(collapse_AC))) {
          stop(
            "cannot use the name of another assessment criterion in the names of \n", 
            "collapse_AC",
            .call = FALSE
          )
        }
      }
      
      
    } else {

      stop(
        "collapse_AC specified, but no assessment criteria were used in the assessment", 
        call. = FALSE
      )
      
    }

  }
      

  ## augment timeSeries structure

  # get determinand group
  # NB detGroup is currently a character if determinandGroups is null and a 
  # factor otherwise - the use of the factor assists with ordering, but probably
  # not needed here as this is primarily an OHAT requirement - have raised
  # an issue to tidy this up
  
  timeSeries <- rownames_to_column(timeSeries, ".series")
  
  timeSeries$detGroup <- ctsm_get_info(
    info$determinand, 
    timeSeries$determinand, 
    "group", 
    info$compartment, 
    sep = "_"
  )

  if (!is.null(determinandGroups)) {
    
    if (!all(timeSeries$detGroup %in% determinandGroups$levels)) {
      stop('some determinand groups present in data, but not in groups argument')
    }
    
    timeSeries$detGroup <- factor(
      timeSeries$detGroup, 
      levels = determinandGroups$levels, 
      labels = determinandGroups$labels, 
      ordered = TRUE
    )

    timeSeries$detGroup <-   timeSeries$detGroup[, drop = TRUE]
  }  
    
  
  # merge stations with timeSeries
  
  timeSeries <- left_join(
    timeSeries, 
    assessment_obj$stations,
    by = "station_code"
  )
  
  timeSeries <- column_to_rownames(timeSeries, ".series")
  

  ## get summary from assessment 
  
  summary <- ctsm_summary_overview(
    assessment, timeSeries, info, classColour, fullSummary = TRUE
  )
  
  summary <- cbind(timeSeries, summary)
    
  summary$series <- row.names(summary)
    

  # double check no legacy data 
  
  if (any(is.na(summary$shape))) {
    stop('some legacy data have crept through')
  }
    

  ## tidy up output
  
  # reorder variables 
  
  wk <- c(
    "series", 
    info$region_id,  
    "country", "CMA", 
    "station_code", "station_name", "station_longname", 
    "station_latitude", "station_longitude", "station_type", "waterbody_type", 
    "determinand", "detGroup", "species", "filtered",
    "submedia", "matrix", "basis", "unit", "sex", "method_analysis", "subseries", 
    "shape", "colour"
  ) 
  
  summary <- summary[c(wk[wk %in% names(summary)], setdiff(names(summary), wk))]

  sortID <- intersect(
    c(info$region_id, "country", "CMA", "station_name", 
      "species", "detGroup", "determinand", "matrix"), 
    names(summary)
  )
  summary <- summary[do.call(order, summary[sortID]), ]
  
  
  # rename variables
  
  summary <- rename(
    summary, 
    determinand_group = detGroup, 
    n_year_all = nyall,
    n_year_fit = nyfit,
    n_year_positive = nypos,
    first_year_all = firstYearAll,
    first_year_fit = firstYearFit,
    last_year = lastyear,
    p_linear_trend = pltrend,
    linear_trend = ltrend,
    p_recent_trend = prtrend,
    recent_trend = rtrend,
    detectable_trend = dtrend,
    mean_last_year = meanLY,
    climit_last_year = clLY
  )
  
  if ("class" %in% names(summary)) {
    summary <- rename(summary, imposex_class = class)
  }
  
  names(summary) <- gsub("diff$", "_diff", names(summary))
  names(summary) <- gsub("achieved$", "_achieved", names(summary))
  names(summary) <- gsub("below$", "_below", names(summary))
  
  
  ## simplify (collapse) AC output if required
  
  if (is_AC && !is.null(collapse_AC)) {

    # augment collapse_AC with unspecified AC
    
    other_AC <- !AC %in% unlist(collapse_AC)

    if (any(other_AC)) {
      other_AC <- AC[other_AC]
      
      extra <- as.list(other_AC)
      names(extra) <- other_AC
      
      collapse_AC <- c(collapse_AC, extra)
      
      collapse_AC <- collapse_AC[sort(names(collapse_AC))]
    }
    

    # set up type and value variables for each AC

    id <- paste0(AC, "_type")
    summary[id] <- lapply(AC, function(x) {
      ifelse(!is.na(summary[[x]]), x, NA_character_)
    })
    
    id <- match(AC, names(summary))
    names(summary)[id] <- paste0(AC, "_value")
    

    var_id <- c("type", "value", "diff", "achieved", "below") 

        
    for (group_id in names(collapse_AC)) {
      
      AC_id <- collapse_AC[[group_id]]
      
      # if only one AC and it has a different name, then rename
      
      if (length(AC_id) == 1 && group_id != AC_id) {
        pos <- match(paste(AC_id, var_id, sep = "_"), names(summary))
        names(summary)[pos] <- paste(group_id, var_id, sep = "_")
      }
      
      if (length(AC_id) > 1) {
      
        in_id <- paste0(AC_id, "_type")
        out_id <- paste0(group_id, "_type")
        summary[out_id] <- apply(summary[in_id], 1, ctsm_collapse_AC, type = "character")
        
        in_id <- paste0(AC_id, "_value")
        out_id <- paste0(group_id, "_value")
        summary[out_id] <- apply(summary[in_id], 1, ctsm_collapse_AC, type = "real")
        
        in_id <- paste0(AC_id, "_diff")
        out_id <- paste0(group_id, "_diff")
        summary[out_id] <- apply(summary[in_id], 1, ctsm_collapse_AC, type = "real")
        
        in_id <- paste0(AC_id, "_achieved")
        out_id <- paste0(group_id, "_achieved")
        summary[out_id] <- apply(summary[in_id], 1, ctsm_collapse_AC, type = "real")
        
        in_id <- paste0(AC_id, "_below")
        out_id <- paste0(group_id, "_below")
        summary[out_id] <- apply(summary[in_id], 1, ctsm_collapse_AC, type = "character")
      }
      
      # remove unwanted columns
      
      id <- setdiff(AC_id, group_id)
      id <- paste(rep(id, each = 5), var_id, sep = "_")
    
      id <- setdiff(names(summary), id)
    
      summary <- summary[id]
    }    
    
    # reorder column names 
    
    id <- paste(rep(names(collapse_AC), each = 5), var_id, sep = "_")
    
    summary <- dplyr::relocate(
      summary, 
      tidyselect::all_of(id), 
      .after = climit_last_year
    )
    
  }
  
  
  # rename region variables if required
  
  if (!identical(info$region_id, info$region_names)) {
    pos <- match(info$region_id, names(summary))
    names(summary)[pos] <- info$region_names
  }
  

  # write summary to output_file or return summary object
    
  if (export) {
    readr::write_excel_csv(summary, output_file, na = "")
    return(invisible())
  } else {
    return(summary)
  }
    
}


ctsm_collapse_AC <- function(x, type = c("real", "character")) {
  type = match.arg(type)
  if (all(is.na(x))) {
    out <- switch(
      type, 
      real = NA_real_, 
      character = NA_character_
    )
    return(out)
  }
  x <- x[!is.na(x)]
  if (length(x) > 1) stop("multiple values not allowed")
  x
}




# OHAT ----

#' @export
ctsm_OHAT_legends <- function(
  assessments, determinandGroups, regionalGroups = NULL, distanceGroups = NULL, path) {

  out <- sapply(names(assessments), simplify = FALSE, USE.NAMES = TRUE, FUN = function(media) {

    assessment.ob <- assessments[[media]]
    assessment <- assessment.ob$assessment
    classColour <- assessment.ob$classColour
    determinands <- assessment.ob$determinands
    regionalGroups <- regionalGroups[[media]]
    distanceGroups <- distanceGroups[[media]]
    
    legends <- ctsm.web.AC(assessment, classColour)
        
    determinands <- intersect(determinands, rownames(legends))   # need this to keep ordering right
    legends <- legends[determinands, , drop = FALSE]   

    compartment <- assessment$info$compartment
    group <- ctsm_get_info(
      info$determinand, determinands, "group", compartment, sep = "_"
    )
    web_group <- factor(
      group, 
      levels = determinandGroups$levels, 
      labels = determinandGroups$labels
    )
    web_group <- wk_group[, drop = TRUE]
    
    goodStatus <- ctsm_get_info(info$determinand, determinands, "good_status")
    goodStatus <- as.character(goodStatus)

    legendName <- apply(legends, 1, function(i) paste(colnames(legends)[i], collapse = " "))
    legendName <- paste(compartment, group, goodStatus, legendName, sep = " ")

    legends <- data.frame(legends, legendName, group, web_group, goodStatus, stringsAsFactors = FALSE)

    ctsm_OHAT_add_legends(legends, classColour, regionalGroups, distanceGroups, assessment$info)
  })
  
  legends <- lapply(out, "[[", "legends") %>% 
    bind_rows(.id = "Compartment")
  
  help <- lapply(out, "[[", "help") %>% 
    bind_rows(.id = "Compartment")
  
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

  standard_shape <- standard_shape %>% 
    bind_rows() %>% 
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
    
    out <- mutate_if(out, is.factor, as.character)
    
    out <- bind_rows(out, standard_shape[[statusID]])
    
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
    
    out <- bind_rows(info_regional, info_distance, help_info)
    
    out <- as.data.frame(out)
    
    out
  })
    
  names(legend_info) <- names(help_info) <- row.names(legends)
  
  legend_info <- bind_rows(legend_info, .id = "Determinand_code")
  help_info <- bind_rows(help_info, .id = "Determinand_code")
  
  list(legends = legend_info, help = help_info)
}
