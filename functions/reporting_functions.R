# web construction ----

ctsm_web_initialise <- function(
  assessmentObject, 
  classColour = NULL, 
  determinandGroups) {

  # reporting_functions.R


  # check assessment criteria have an appropriate class colour
  
  AC <- assessmentObject$info$AC
  if (!is.null(AC)) 
    stopifnot(
      AC %in% names(classColour$below), 
      AC %in% names(classColour$above),
      "none" %in% names(classColour)
    )
  
  
  # create submedia (levelTwo) and levelSix and levelSeven variables 
  # (should be done externally) and 'species' (levelThree) for sediment
  
  compartment <- assessmentObject$info$compartment
  purpose <- assessmentObject$info$purpose
  
  assessmentObject$timeSeries <- within(
    assessmentObject$timeSeries, 
    switch(
      compartment, 
      biota = {
        level2name <- ctsm_get_info("species", species, "family")[drop = TRUE]
        levels(level2name) <- list(
          "Fish" = "Fish", 
          "Shellfish" = c("Bivalvia", "Gastropoda"), 
          "Crustacea" = "Crustacea", 
          "Bird" = "Bird",
          "Mammal" = "Mammal")
        level2element <- "Family"
        level3name <- ctsm_get_info("species", species, "common.name")
        level3element <- "Species"
      }, 
      sediment = {                                                       
        # null values 
        level2name <- "Sediment"   
        level2element <- as.character(NA)
        level3name <- "Sediment"   
        level3element <- as.character(NA)
      }, 
      water = {                                                       
        # null values 
        switch(
          purpose, 
          HELCOM = {
            level2name <- ifelse(filtered, "filtered", "unfiltered")
            level2element <- "Filtration"
          }, {
            level2name <- "Water"
            level3element <- "Water"
          }
        )
        level3name <- "Water"   
        level3element <- as.character(NA)
      }
    ))
  
  

  # create detGroup variable in each timeSeries structure

  assessmentObject$timeSeries <- within(
    assessmentObject$timeSeries, {
      detGroup <- ctsm_get_info(
        "determinand", determinand, "group", compartment, sep = "_"
      )
      if (!all(detGroup %in% determinandGroups$levels)) 
        stop('some determinand groups present in data, but not in groups argument')
      detGroup <- factor(
        detGroup, 
        levels = determinandGroups$levels, labels = determinandGroups$labels, 
        ordered = TRUE)
      detGroup <- detGroup[, drop = TRUE]
  })
  
  assessmentObject <- ctsm_web_setup(assessmentObject)
  
  list(
    assessment = assessmentObject, 
    # determinands = determinands, 
    classColour = classColour
  )
}


ctsm_web_setup <- function(assessmentObject) {
  
  # extract timeSeries and stations and get more useful names
  
  timeSeries <- assessmentObject$timeSeries
  timeSeries <- swap.names(timeSeries, "station", "stationID")
  timeSeries <- within(timeSeries, stationID <- as.character(stationID))
  
  stations <- assessmentObject$stations
  stations <- swap.names(stations, "name", "stationName")
  
  
  # convert station factors to characters for output 
  
  stations <- tibble::rownames_to_column(stations, ".rownames")
  stations <- dplyr::mutate(stations, dplyr::across(where(is.factor), as.character))
  stations <- tibble::column_to_rownames(stations, ".rownames")
  
  
  # check lower case stationIDs are unique (for unique filenames)
  
  wk <- tolower(row.names(stations))
  if (any(duplicated(wk))) {
    stop(
      "stationID not unique when converted to lower case ", 
      "- might cause difficulties with filenames"
    )
  }
  

  # merge timeSeries and stations  
  
  timeSeries <- cbind(timeSeries, stations[timeSeries$stationID, ])

  assessmentObject$timeSeries <- timeSeries
  assessmentObject
}



# mapping functions ----

ctsm.Mercator <- function(x) atanh(sin(x * pi / 180))

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


ctsm.web.overview <- function(assessmentObject, classColour, fullSummary = FALSE) {
  # gets shape and colour for each time series
  
  assessment <- assessmentObject$assessment
  info <- assessmentObject$info
  timeSeries <- assessmentObject$timeSeries
  
  # first get list of assessment summaries - some are null, so need to be careful
  # also summary structures differ between detGroups, so need to be doubly careful
  
  summaryList <- sapply(row.names(timeSeries), simplify = FALSE, USE.NAMES = TRUE, 
                        FUN = function(i) assessment[[i]]$summary)
  
  
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
  
  
  # get shape of plotting symbols
  
  out$shape <- with(out, {
    
    shape <- character(nrow(out))
    
    # trend symbols: a trend is estimated if pltrend is present - default shape is a large filled circle
    
    trendFit <- !is.na(pltrend)
    shape[trendFit] <- "large_filled_circle"
    
    # show a significant trend based on prtrend 
    # NB prtrend might not exist even if pltrend does, because there are too few years of data in the
    # recent window
    
    # isImposex <- timeSeries$detGroup %in% "Imposex"
    
    isTrend <- !is.na(prtrend) & prtrend < 0.05
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
    
    goodStatus <- ctsm_get_info("determinand", timeSeries$determinand, "good_status")
    
    wk <- out[ACdiff]
    wk[] <- lapply(wk, "*", ifelse(goodStatus == "low", 1, -1))
    
    out$colour <- apply(wk, 1, function(x) {
      
      if (all(is.na(x))) return(classColour$none)
      
      AC <- info$AC[!is.na(x)]
      x <- x[!is.na(x)]
      
      if (any(x < 0)) classColour$below[AC[which.max(x < 0)]]
      else classColour$above[AC[length(x)]]
    })  
    
    # need to adjust for null summaries
    
    out <- within(out, colour[is.na(shape)] <- NA)
    
  } else {
    
    out$colour <- classColour$none
    
  }
  
  
  # adjust shape and colour for nonparametric test if nyfit <= 2 and nyall > 2
  
  # get the names of the variables which contain the result of the non-parametric test for each AC
  
  if (!is.null(info$AC)) {
    
    ACbelow <- paste(info$AC, "below", sep = "")
    
    wk <- out[ACbelow]
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
    
    id <- with(out, nyfit <= 2 & nyall > 2 & !is.na(wk))
    out$colour[id] <- wk[id]
    out$shape[id] <- "small_filled_circle"
    
  }
  
  if (fullSummary) out else out[c("shape", "colour")]
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

ctsm_summary_table <- function(
  assessments, determinandGroups, path = ".", export = TRUE,  
  include_all = FALSE, collapse_AC = TRUE) {

  # reporting functions
  # get summary files
  
  # path is for output

  library(dplyr)

  out <- sapply(names(assessments), simplify = FALSE, FUN = function(x) {
  
    assessment <- assessments[[x]]$assessment
    classColour <- assessments[[x]]$classColour
    purpose <- assessment$info$purpose
    compartment <- assessment$info$compartment

    summary <- ctsm.web.overview(assessment, classColour, fullSummary = TRUE)

    summary <- cbind(assessment$timeSeries, summary)
    
    summary$series <- row.names(summary)
    

    # double check no legacy data and exclude those columns that aren't useful

    if (any(is.na(summary$shape))) 
      stop('some legacy data have crept through')
    
    if (!include_all) {
      notok <- match(
        c("stationID", 
          "level6element", "level6name", "level7element", "level7name", 
          "level2element", "level2name", "level3element", "level3name", 
          "fileName", "offshore"), 
        names(summary)
      )
      notok <- na.omit(notok)
      if (length(notok)) summary <- summary[, -notok]
    }
      
    # reorder to get most useful output
  
    wk <- c(
      "series", 
      assessment$info$region_id,  
      "country", "CMA", 
      "code", "station", "stationName", "latitude", "longitude", "MSTAT", "WLTYP", 
      "determinand", "detGroup", "species", "filtered",
      "submedia", "matrix", "basis", "unit", "sex", "metoa", "AMAP_group", "shape", "colour"
    ) 
  
    summary <- summary[c(wk[wk %in% names(summary)], setdiff(names(summary), wk))]

    sortID <- intersect(
      c(assessment$info$region_id, "country", "CMA", "station", 
        "species", "detGroup", "determinand", "matrix"), 
      names(summary)
    )
    summary <- summary[do.call(order, summary[sortID]), ]
    
      
    
    # now streamline to get colour coded fit
    
    wk <- c(
      assessment$info$region_id, 
      "country", "CMA", 
      "station", "stationName", "latitude", "longitude", "determinand", "detGroup", 
      "species", "matrix", "sex", "metoa", "AMAP_group", "filtered", "shape", "colour", 
      "nyfit", "nypos"
    )
  
    overview <- summary[wk[wk %in% names(summary)]]

    overview <- within(overview, determinand2 <- as.character(determinand))
    
    # ad hoc fix for biota - need to generalise

    paste_overview <- function(det2, var) ifelse(is.na(var), det2, paste(det2, var))
    
    if ("matrix" %in% names(overview))
      overview <- within(overview, determinand2 <- paste_overview(determinand2, matrix))

    if ("sex" %in% names(overview))
      overview <- within(overview, determinand2 <- paste_overview(determinand2, sex))
    
    if ("metoa" %in% names(overview))
      overview <- within(overview, determinand2 <- paste_overview(determinand2, metoa))
    
    if ("AMAP_group" %in% names(overview))
      overview <- within(overview, determinand2 <- paste_overview(determinand2, AMAP_group))
    
    overview <- within(overview, {
      shape[shape %in% "downward_triangle"] <- "down"
      shape[shape %in% "upward_triangle"] <- "up"
      shape[shape %in% c("large_filled_circle", "small_filled_circle", "small_open_circle")] <- 
        "flat"
      summary <- paste(nyfit, colour, shape)
    })

    overview <- within(overview, id <- as.character(station))
    
    if ("country" %in% names(overview)) overview <- within(overview, id <- paste(country, id))
    
    if ("species" %in% names(overview)) overview <- within(overview, id <- paste(id, species))
    
    if ("filtered" %in% names(overview)) overview <- within(overview, id <- paste(id, filtered))
    
    
    station.var <- c(
      assessment$info$region_id, 
      "country", "CMA", 
      "station", "stationName", "latitude", "longitude", "species", "filtered"
    )
    station.var <- station.var[station.var %in% names(summary)]
    
    overview <- overview[c(station.var, "determinand2", "summary", "id")] 

    overview <- reshape(overview, direction = "wide", v.names = "summary", 
              timevar = "determinand2", idvar = "id")

    overview <- overview[-match("id", names(overview))]

    names(overview) <- gsub("summary.", "", names(overview), fixed = TRUE)
    
    
    # order both rows and columns

    sortID <- intersect(
      c(assessment$info$region_id, "CMA", "station"), 
      names(overview)
    )
    overview <- overview[do.call(order, overview[sortID]), ]
    
    wk <- setdiff(names(overview), station.var)
    wk.det <- sapply(strsplit(wk, " "), "[[", 1)
    wk.group <- ctsm_get_info(
      "determinand", wk.det, "group", assessment$info$compartment, sep = "_"
    )
    wk.group <- factor(
      wk.group, 
      levels = c(
        "Metals", "Organotins", "PAH_parent", "PAH_alkylated", "PAH_metabolites", 
        "PBDEs", "Organobromines", "Organofluorines", "Chlorobiphenyls", 
        "Dioxins", "Organochlorines", "Imposex", "Effects"
      )
    )  
    wk.group <- droplevels(wk.group)
    
    overview <- overview[c(station.var, wk[order(wk.group, wk.det)])]

    if (!export) return(list(summary = summary, overview = overview))
    
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
    
    if ("class" %in% names(summary))
      summary <- rename(summary, imposex_class = class)
    
    names(summary) <- gsub("diff$", "_diff", names(summary))
    names(summary) <- gsub("achieved$", "_achieved", names(summary))
    names(summary) <- gsub("below$", "_below", names(summary))
    
    
    # put all BAC type results together and all EAC type results together
    # this is developmental and, if there are alternative BAC types (AMAP NRC), 
    # they will have to be edited externally for now
    
    if (collapse_AC) {

      AC <- assessment$info$AC
      
      # set up basic type and value variables for each AC
      # BAC_type is what we want to keep
      # EAC_type needs to be collapsed across the remaining EAC related variables
      
      id <- paste0(AC, "_type")
      summary[id] <- lapply(AC, function(x) {
        ifelse(!is.na(summary[[x]]), x, NA_character_)
      })
      
      id <- match(AC, names(summary))
      names(summary)[id] <- paste0(AC, "_value")
      
      
      EAC_id <- setdiff(AC, "BAC")
      
      if (length(EAC_id) >= 1) {

        collapse <- function(x, type = c("real", "character")) {
          type = match.arg(type)
          if (all(is.na(x))) {
            return(switch(type, real = NA_real_, character = NA_character_))
          }
          x <- x[!is.na(x)]
          if (length(x) > 1) stop("multiple EAC values not allowed")
          x
        }
        
        id <- paste0(EAC_id, "_type")
        summary["EAC_type"] <- apply(summary[id], 1, collapse, type = "character")

        id <- paste0(EAC_id, "_value")
        summary["EAC_value"] <- apply(summary[id], 1, collapse)
        
        id <- paste0(EAC_id, "_diff")
        summary["EAC_diff"] <- apply(summary[id], 1, collapse)
        
        id <- paste0(EAC_id, "_achieved")
        summary["EAC_achieved"] <- apply(summary[id], 1, collapse)
        
        id <- paste0(EAC_id, "_below")
        summary["EAC_below"] <- apply(summary[id], 1, collapse, type = "character")
      }
      
      # remove unwanted columns

      id <- setdiff(EAC_id, "EAC")
      id <- paste(
        rep(id, each = 5), 
        c("type", "value", "diff", "achieved", "below"), 
        sep = "_"
      )
      
      id <- setdiff(names(summary), id)
      
      summary <- summary[id]
      
      # reorder so that BAC_type and EAC type are first in their class
      
      id <- paste(
        rep(c("BAC", "EAC"), each = 5), 
        c("type", "value", "diff", "achieved", "below"), 
        sep = "_"
      )
      
      summary <- dplyr::relocate(
        summary, 
        tidyselect::any_of(id), 
        .after = climit_last_year
      )
    }
    

    # rename some variables
    
    summary <- rename(
      summary, 
      station_code = code,
      station_name = station,
      station_long_name = stationName,
    )
    
    if ("AMAP_group" %in% names(summary)) {
      summary <- rename(summary, mammal_group = AMAP_group)
    }
    
    if (!identical(assessment$info$region_id, assessment$info$region_names)) {
      pos <- match(assessment$info$region_id, names(summary))
      names(summary)[pos] <- assessment$info$region_names
    }
    

    outfile <- file.path(path, paste0(tolower(x), "_summary.csv"))
     
    readr::write_excel_csv(summary, outfile, na = "")
  })  
  
  if (export) invisible() else out
}


# OHAT ----

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
      "determinand", determinands, "group", compartment, sep = "_"
    )
    web_group <- factor(group, levels = determinandGroups$levels, labels = determinandGroups$labels)[, drop = TRUE]
    
    goodStatus <- ctsm_get_info("determinand", determinands, "good_status")
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
