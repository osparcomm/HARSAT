# main assessment functions

# Set up functions ----

#' Assess timeseries for trends and status
#'
#' Fits a model to each timeseries, test for any temporal trend and compare with
#' thresholds. Need to add a lot more in details.
#'
#' @param ctsm_ob A HARSAT object resulting from a call to create_timeSeries
#' @param subset An optional vector specifying which timeseries are to be
#'   assessed. Might be used if the assessment is to be done in chunks because
#'   of size, or when refitting a timeseries model which has not converged. An
#'   expression will be evaluated in the timeSeries component of ctsm_ob; use
#'   'series' to identify individual timeseries.
#' @param AC A character vector identifying the thresholds to be used in status
#'   assessments. These should be in the threshold reference table. Defaults to
#'   NULL; i.e. no thresholds are used.
#' @param get_AC_fn An optional function that overrides get_AC_default. See
#'   details (which need to be written).
#' @param recent_trend An integer giving the number of years which are used in
#'   the assessment of recent trends. For example, a value of 20 (the default)
#'   consider trends in the last twenty year.
#' @param parallel A logical which determines whether to use parallel
#'   computation; default = FALSE.
#' @param extra_data `r lifecycle::badge("experimental")` A named list used to
#'   pass additional data to specific assessment routines. At present it is 
#'   only used for imposex assessments, where it passes two data frames called
#'   `VDS_estimates` and `VDS_confidence_limits`. Defaults to NULL, This 
#'   argument will be generalised in the near future, so expect it to change.
#' @param control `r lifecycle::badge("experimental")` A list of control 
#'   parameters that allow the user to modify the way the assessment is run. 
#'   At present, these only include parameters involved in post-hoc power 
#'   calculations, but it is intended to move other structures such as 
#'   `recent_trend` here. See details (which need to be written).
#' @param ... Extra arguments which are passed to assessment_engine. See
#'   details (which need to be written).
#' @export
run_assessment <- function(
  ctsm_ob, 
  subset = NULL, 
  AC = NULL, 
  get_AC_fn = NULL, 
  recent_trend = 20L, 
  parallel = FALSE,
  extra_data = NULL,
  control = list(),
  ...) {
  
  # location: assessment_functions.R

  # set up assessment structure
  
  ctsm_ob$call <- match.call()

  ctsm_ob$info$recent.trend <- recent_trend

  ctsm_ob$info$AC <- AC
  
  ctsm_ob$info$get_AC_fn <- get_AC_fn
  if (!is.null(AC) && is.null(ctsm_ob$info$get_AC_fn)) {
    ctsm_ob$info$get_AC_fn <- get_AC[[ctsm_ob$info$compartment]]
  }
  
  if (any(ctsm_ob$data$group %in% "Imposex")) {
    if (is.null(extra_data)) {
      stop("`extra_data` must be supplied for imposex assessments")
    }
    
    ok <- c("VDS_estimates", "VDS_confidence_limits") %in% names(extra_data)
    if (!all(ok)) {
      stop(
        "argument extra_data must be a list with components ", 
        "VDS_estimates and VDS_confidence_limits"
      )
    }
  }
  
  ctsm_ob$info$extra_data <- extra_data

  
  # update control information
  
  cntrl <- run_control_default()
  
  cntrl <- run_control_modify(cntrl, control)
  
  if (any(names(cntrl) %in% names(ctsm_ob$info))) {
    id <- names(cntrl)
    id <- id[id %in% names(ctsm_ob$info)]
    warning(
      "\n conflict between components of ctsm_ob$info and control parameters ", 
      "- results may be unexpected:\n ",
      paste(id, collapse = ", "), "\n",
      call. = FALSE, immediate. = TRUE)
  }
  
  ctsm_ob$info <- append(ctsm_ob$info, cntrl)
  

  ctsm_ob$assessment <- vector(mode = "list", length = nrow(ctsm_ob$timeSeries))
  names(ctsm_ob$assessment) <- row.names(ctsm_ob$timeSeries)

  ctsm_ob$call.data <- NULL

  
  # identify which series are to be assessed in this run - defaults to all
  
  series_id <- row.names(ctsm_ob$timeSeries)
  
  if (!is.null(substitute(subset))) {
    timeSeries <- tibble::rownames_to_column(ctsm_ob$timeSeries, "series")
    ok <- eval(substitute(subset), timeSeries, parent.frame())
    series_id <- timeSeries[ok, "series"]
  }

  
  # assess timeseries
  
  out <- assessment_engine(
    ctsm_ob, 
    series_id,
    parallel = parallel, 
    ...
  )
  
  
  # organise output
  
  ctsm_ob$assessment[names(out)] <- out
  
  ctsm_ob
}


#' Update timeseries assessments
#'
#' Refits models for particular timeseries, or fits new models when an
#' assessment is being done in chunks.
#'
#' @param ctsm_ob A HARSAT object resulting from a call to run_assessment
#' @param subset A vector specifying which timeseries assessements are to be
#'   updated or fit for the first time. An expression will be evaluated in the
#'   timeSeries component of ctsm_ob; use 'series' to identify individual
#'   timeseries.
#' @param parallel A logical which determines whether to use parallel
#'   computation; default = FALSE.
#' @param ... Extra arguments which are passed to assessment_engine.  See
#'   details (which need to be written).
#' @export
update_assessment <- function(ctsm_ob, subset = NULL, parallel = FALSE, ...) {
  
  # location: assessment_functions.R
  
  # error trapping
  
  # ctsm_ob must be a valid assessment object
  
  if (!"assessment" %in% names(ctsm_ob) || 
      !identical(names(ctsm_ob$assessment), row.names(ctsm_ob$timeSeries))) {
    stop(
      "'ctsm_ob' is not a valid harsat assessment object as it is not the ",
      "result of a call\n  to run_assessment.")
  }
  
  # AC cannot be passed as an argument, as it already specified in 
  # run_assessment (easily copied over by mistake)
  
  if ("AC" %in% names(match.call())) {
    stop(
      "'AC' cannot be an argument of update_assessment as the thresholds were ",
      "specified in\n  an earlier call to run_assessment and cannot be changed."
    )
  }
  

  # identify which series are to be assessed in this run
  
  if (is.null(substitute(subset))) {
    stop("subset is NULL: nothing to update")
  }
  
  timeSeries <- tibble::rownames_to_column(ctsm_ob$timeSeries, "series")
  ok <- eval(substitute(subset), timeSeries, parent.frame())
  series_id <- timeSeries[ok, "series"]

  
  # assess timeseries
  
  out <- assessment_engine(
    ctsm_ob, 
    series_id,
    parallel = parallel, 
    ...
  )
  
  
  # organise output
  
  ctsm_ob$assessment[names(out)] <- out
  
  ctsm_ob
}



assessment_engine <- function(ctsm.ob, series_id, parallel = FALSE, ...) {

  # silence non-standard evaluation warnings
  .data <- index <- NULL

  # location: assessment_functions.R
  
  # assess each time series in turn
  
  info <- ctsm.ob$info

  timeSeries <- ctsm.ob$timeSeries[series_id, ]
  
  data <- dplyr::filter(ctsm.ob$data, .data$seriesID %in% row.names(timeSeries))
  data <- droplevels(data)

  stations <- tibble::column_to_rownames(ctsm.ob$stations, "station_code")
  stations <- stations[unique(timeSeries$station_code), ]

  
  
  # set up parallel processing information 
  
  if (parallel) {

    cat("Setting up clusters: be patient, it can take a while!")
        
    n_cores <- parallel::detectCores()
    cluster_id <- parallel::makeCluster(n_cores - 1, outfile = "")
    on.exit(parallel::stopCluster(cluster_id))

    is_imposex <- "Imposex" %in% data$group
    export_objects <- parallel_objects(is_imposex)
    package_environment <- environment(sys.function(sys.nframe()))
    parallel::clusterExport(cluster_id, export_objects, envir = package_environment)

  } else {
    cluster_id <- NULL
  }


  # split data into each series
  
  data <- split(data, data$seriesID)

  assessment <- pbapply::pblapply(data, ..., cl = cluster_id, FUN = function(x, ...) {
  
    # get info about the time series
    
    seriesID <- x$seriesID[1]
    seriesInfo <- sapply(
      timeSeries[seriesID,], 
      function(i) if (is.numeric(i)) i else as.character(i), 
      USE.NAMES = TRUE, 
      simplify = FALSE
    )
    seriesInfo <- seriesInfo[vapply(seriesInfo, function(x) !is.na(x), NA)]
    
    cat(
      "\nassessing series: ", 
      paste(names(seriesInfo), seriesInfo, collapse = "; "), 
      "\n"
    )
    
    
    # return if no data to assess
    
    if (all(is.na(x$concentration))) return()
    

    # add station information to seriesInfo
    
    seriesInfo <- c(seriesInfo, stations[seriesInfo$station_code, ])


    # construct annual index
    
    determinand <- seriesInfo$determinand
    
    var_id <- c("year", "concentration")
    
    var_id <- if (determinand %in% c("VDS", "IMPS", "INTS")) {
      append(var_id, c("n_individual", "%FEMALEPOP"))
    } else {
      append(var_id, c("censoring", "uncertainty"))
    }
    
    if (determinand %in% "MNC") {
      var_id <- append(var_id, "MNC-QC-NR")
    }

    if (determinand %in% "%DNATAIL") {
      var_id <- append(var_id, "CMT-QC-NR")
    }
    
    x <- x[!is.na(x$concentration), var_id]
    
    row.names(x) <- NULL

    annualIndex <- get_index(determinand, x, info)
    
    # get assessment concentrations - need to extract some key variables first 
    # could streamline in future
    
    if ("AC" %in% names(info)) {
      rt_id <- c("thresholds", "determinand", "species")
      rt_id <- rt_id[rt_id %in% names(info)]
      AC <- info$get_AC_fn(
        data = as.data.frame(seriesInfo), 
        AC = info$AC,
        rt = info[rt_id]
      )
      AC <- unlist(AC)
    } else { 
      AC <- NULL
    }
    
    # initialise output from function
    
    out <- list(fullData = x, annualIndex = annualIndex, AC = AC)
    
    
    # do assessment
    
    if (determinand %in% c("VDS", "IMPS", "INTS")) {
      
      species <- seriesInfo$species
      station_code <- seriesInfo$station_code

      theta_id <- switch(
        info$purpose, 
        # CSEMP = stations[station_code, "CMA"],
        OSPAR = paste(seriesInfo$country, stations[station_code, "ospar_subregion"]),
        HELCOM = seriesInfo$country
      )

      theta_id <- paste(theta_id, species)
      
      is_theta <- theta_id %in% names(info$extra_data$VDS_estimates)

      theta <- NULL
      if (is_theta) {
        theta <- info$extra_data$VDS_estimates[[theta_id]]
      }
      
      
      # if any individual data, need to augment annual indices with confidence 
      # intervals
      
      indiID <- with(x, tapply(n_individual, year, function(y) all(y == 1)))
      
      if (any(indiID) & is_theta) {

        out$annualIndex[c("lower", "upper")] <- NA
        
        clID <- paste(station_code, names(indiID)[indiID], species)
        out$annualIndex[indiID, c("lower", "upper")] <- 
          info$extra_data$VDS_confidence_limits[clID, c("lower", "upper")]
        
        # adjust when indices are zero or max (K = #cutpoints)
        # had to add an observation with a 1 or n-1 to get a fit
        
        out$annualIndex <- within(out$annualIndex, {
          lower[index == 0] <- 0
          upper[index == theta$K] <- theta$K
        })  
      }  
      
      out <- c(out, assess_imposex(
        data = x, 
        annualIndex = out$annualIndex, 
        AC = AC, 
        recent.years = info$recent_years,
        determinand = determinand, 
        species = species,
        station_code = station_code,
        theta = theta, 
        max.year = info$max_year, 
        info.imposex = info$imposex, 
        recent.trend = info$recent.trend)
      )
      
      out$convergence <- 0L
      
      out
    }
    else {

      # ad-hoc fix for missing uncertainties for SFG: 
      # trap for any lognormal or normal distributed data that have missing 
      #   uncertainties

      distribution <- ctsm_get_info(info$determinand, determinand, "distribution")
      good.status <- ctsm_get_info(info$determinand, determinand, "good_status")

      if (determinand %in% "SFG" && any(is.na(x$uncertainty))) {      
        cat("  warning: ad-hoc fix to estimate missing uncertainties\n")
        pos <- is.na(x$uncertainty)
        x$uncertainty[pos] <- 0.1
      }

      if (any(is.na(x$uncertainty)) && distribution %in% c("lognormal", "normal")) {
        stop("missing uncertainties not allowed")
      }
      
            
      args.list <- list(
        data = x, 
        annualIndex = annualIndex,
        AC = AC, 
        recent.years = info$recent_years, 
        determinand = determinand, 
        max.year = info$max_year, 
        recent.trend = info$recent.trend, 
        distribution = distribution, 
        good.status = good.status,
        power = info$power
      )
      
      args.list <- c(args.list, list(...))
      fit <- try(do.call("assess_lmm", args.list))
      if (is.character(fit) & length(fit) == 1L) {
        out <- c(out, error = fit)
      } else {
        out <- c(out, fit)
      }
      
      out$convergence <- check_convergence_lmm(out)
      
      out
    }
    
  })
  
  assessment
}



#' Default control parameters for `run_assessment`
#'
#' `r lifecycle::badge("experimental")` Default parameters that control the way
#'  the assessment is run. Presently only includes parameters for post-hoc 
#'  power, but it is intended to move `recent_trend` here, along with the 
#'  arguments that control the calculation of numerical derivatives.  
#'
#' @returns A list with the following components:
#' * `power` A list with the following components (all expressed as 
#'   percentages):  
#'   - `target_power` default = 90%
#'   - `target_trend` default = 5%
#'   - `size` default = 5%  
#'   The power calculations are currently only applied to log-normally 
#'   distributed data, which is why the trend is expressed as a percentage.
#'
run_control_default <- function() {

  # location: assessment_functions.R
  
  power = list(
    target_power = 90,
    target_trend = 5,
    size = 5
  )
  
  list(power = power)  
}





#' Modifies control parameters for `run_assessment`
#'
#' Undates default control parameters with user specification and does basis
#' error checking.
#'
#' @param run_control_default List of default control parameters produced by a 
#'   call to `run_control_default`
#' @param control List of replacement control parameters; defaults to an empty 
#'   list. 
#'
#' @returns List of updated control parameters  
#' 
run_control_modify <- function(control_default, control = list()) {
  
  # location: assessment_functions.R
  
  control <- modifyList(control_default, control, keep.null = TRUE)

  if (control$power$target_power <= control$power$size) {
    stop(
      "error in target_power component of control$power:\n", 
      "target_power must be greater than size"
    )
  }
    
  if (control$power$size <= 0 | control$power$size >= 100) {
    stop(
      "error in size component of control$power:\n", 
      "size (%) must be greater than 0 and less than 100"
    )
  }
  
  if (control$power$target_power >= 100) {
    stop(
      "error in target_power component of control$power:\n", 
      "target_power must be less that 100%"
    )
  }

  if (control$power$target_trend <= -100) {
    stop(
      "error in target_trend component of control$power:\n", 
      "target_trend must be greater than -100%"
    )
  }
  
  control
}





parallel_objects <- function(imposex = FALSE) {
  
  # assessment_functions.R
  # objects required for clusterExport

  package.environment <- environment(sys.function(sys.nframe()))
  
  out <- c(
    "negTwiceLogLik", 
    "check_convergence_lmm",
    objects(package.environment, pattern = "^assess_*"),
    objects(package.environment, pattern = "^get*"),
    objects(package.environment, pattern = "^ctsm*")
  )

  if (imposex) {
    out <- c(
      out, 
      "imposex.assess.index", "imposex_class", 
      "imposex.family", "cuts6.varmean",  
      "imposex_assess_clm", "imposex.clm.fit", "imposex.clm.X", 
      "imposex.clm.X.change", "imposex.clm.loglik.calc", "imposex.VDS.p.calc", 
      "imposex.clm.predict", "imposex.clm.cl", "imposex.clm.contrast"
    )
  }
  
  out
}


# Annual indices ----

# get.index <- function(determinand, data, info) {
#   
#   group <- ctsm_get_info(
#     info$determinand, 
#     determinand, 
#     "group", 
#     info$compartment, 
#     sep = "_"
#   )
#   
#   function_id <- paste("get.index", info$compartment, group, sep = ".")
#   
#   do.call(function_id, list(data = data, determinand = determinand, info = info), envir = sys.frame()) 
# }


# get.index.default <- function(data, determinand, info) {
# 
#   data <- data[!is.na(data$concentration), ]
#   
#   # median (log) concentrations with flag to denote if less thans used in their 
#   # construction
#   
#   distribution <- ctsm_get_info(info$determinand, determinand, "distribution")
#   
#   data$response <- switch(
#     distribution, 
#     lognormal = log(data$concentration), 
#     data$concentration
#   )
#   
#   index <- tapply(data$response, data$year, median, na.rm = TRUE)
#   
#   censoring <- by(data, data$year, function(x) {
#     x <- x[order(x$concentration),]
#     n <- nrow(x)
#     n0 <- ceiling(n / 2)
#     any(x$censoring[n0:n] %in% c("<", "D", "Q"))
#   })
#   censoring <- sapply(censoring, function(i) i)	
#   
#   data <- data.frame(year = as.numeric(names(index)), index, censoring, row.names = NULL)
#   
#   data <- within(data, censoring <- factor(ifelse(censoring, "<", "")))
#   
#   data
# }



get_index <- function(determinand, data, info) {
  
  distribution <- ctsm_get_info(info$determinand, determinand, "distribution")

  good_status <- ctsm_get_info(info$determinand, determinand, "good_status")
  
  
  # crude trap for unexpected censoring values
  
  censoring_trap <- FALSE
  
  if (any(data$censoring != "") && good_status == "high") {
    censoring_trap <- TRUE
  } 
  
  if (any(data$censoring != "") && !distribution %in% c("normal", "lognormal")) {
    censoring_trap <- TRUE
  } 
  
  if (censoring_trap) {
    warning(
      "unexpected censoring values found for determinand ", 
      determinand, "\n",
      "results might be compromised\n", 
      "please consult harsat development team", 
      call. = FALSE
    )
  }

  
  if (distribution %in% "lognormal") {
    get_index_median(data, log = TRUE)
  } else if (distribution %in% c("normal", "survival")) {
    get_index_median(data, log = FALSE)
  } else if (distribution %in% c("multinomial", "quasibinomial")) {
    get_index_imposex(data, determinand, info)
  }else if (distribution %in% c("beta", "negativebinomial")) {
    get_index_weighted_mean(data, determinand) 
  } else {
    stop(
      "distribution ", distribution, " not recognised\n", 
      "please consult harsat development team", 
      call. = FALSE
    )
  }
    
}


get_index_median <- function(data, log = TRUE) {
  
  data <- data[!is.na(data$concentration), ]
  
  # median (log) concentrations with flag to denote if less thans used in their 
  # construction
  
  data$response <- data$concentration
  if (log) {
    data$response <- log(data$response)
  }
  
  index <- tapply(data$response, data$year, median, na.rm = TRUE)
  
  censoring <- by(data, data$year, function(x) {
    x <- x[order(x$concentration),]
    n <- nrow(x)
    n0 <- ceiling(n / 2)
    any(x$censoring[n0:n] %in% c("<", "D", "Q"))
  })
  censoring <- sapply(censoring, function(i) i)	
  censoring <- ifelse(censoring, "<", "")

  data.frame(
    year = as.numeric(names(index)), 
    index, 
    censoring, 
    row.names = NULL
  )
  
}


get_index_weighted_mean <- function(data, determinand) {
  
  # only currently used by MNC and %DNATAIL where we get a weighted average
  # based on the number of cells

  if (!determinand %in% c("MNC", "%DNATAIL")) {
    stop(
      "get_index_weighted_mean can not be used with determinand ", 
      determinand, "\n", 
      "currently only coded for MNC or %DNATAIL\n", 
      "please contact harsat development team", 
      call. = FALSE
    )
  }
    
  data <- data[!is.na(data$concentration), ]
  
  nCell_id <- switch(
    determinand,
    MNC = "MNC-QC-NR",
    "%DNATAIL" = "CMT-QC-NR"
  )
  
  out <- by(data, data$year, function(x) {
    names(x)[match(nCell_id, names(x))] <- "nCell"
    data.frame(
      year = with(x, unique(year)),
      index = with(x, sum(concentration * nCell) / sum(nCell)),
      nCell = with(x, sum(nCell))
    )
  })
  
  out <- do.call("rbind", out)
  
  out$censoring <- rep("", nrow(out))
  
  out <- out[c("year", "index", "censoring", "nCell")]
  
  out$row.names <- NULL
  
  return(out)
}



# Mixed model assessment functions ----

assess_lmm <- function(
    data, annualIndex, AC, recent.years, determinand, max.year, 
    recent.trend = 20, distribution, good.status, choose_model, 
    power, ...) {

  # silence non-standard evaluation warnings
  year <- NULL

  # choose_model forces exit with a particular model: 2 = linear, 3 = smooth on 2df etc, with an 
  # error if that model doesn't exist 
  # this should only be used with extreme care - it is provided for those very rare casese where a 
  # model has been 'over fit' and it is not possible to get standard errors


  # order data

  data <- data[order(data$year), ]
  

  # get total number of years and first year in the time series
  # need to do this here because early years can get stripped out of the data set 
  
  nYearFull <- length(unique(data$year)) 
  firstYearFull <- min(data$year)
  
  
  # deal with data sets that have crept in by mistake and have no recent data
  
  if (max(data$year) < min(recent.years)) return (NULL)
  

  data$response <- switch( 
    distribution,
    lognormal = log(data$concentration), 
    data$concentration
  )
  
  
  # auxiliary variables
    
  if (distribution %in% "lognormal") {
    data$sdAnalytic <- data$uncertainty / data$concentration
  } else if (distribution %in% "normal") {
    data$sdAnalytic <- data$uncertainty
  } else if (determinand %in% "%DNATAIL") {
    data$weight <- data[["CMT-QC-NR"]]
  } else if (determinand %in% "MNC") {
    data$offset <- data[["MNC-QC-NR"]]
  } 

  
  # non-parametric classification (looks at last e.g. ten years of data as defined by recent.trend)

  if (distribution %in% c("lognormal", "normal")) {
    below.result <- lapply(AC, function(i) ctsm.test.below(
      annualIndex$year, 
      annualIndex$index, 
      value = switch(distribution, lognormal = log(i), normal = i), 
      min.year = max.year - recent.trend + 1,
      below = switch(good.status, low = TRUE, high = FALSE)
    ))
  }
  

  # splits data into e.g. 6 year blocks and checks at least one observation in each block - 
  # ensures there aren't humungous gaps in the time series

  reporting.window <- length(recent.years)
  
  year.group <- cut(
    data$year, seq(max.year, min(data$year) - reporting.window, by = - reporting.window))
  
  in.group <- table(year.group) > 0
  
  if (!all(in.group)) {
  
    ok.group <- cumprod(rev(in.group)) == 1
    ok.group <- names(which(ok.group))
    
    id <- year.group %in% ok.group
    data <- data[id, ]
  }

  
  # split off to assess other distributions - need to harmonise this
  
  if (!distribution %in% c("lognormal", "normal")) {
    wk_fn <- paste0("assess_", distribution)
    output <- do.call(
      wk_fn, 
      args = list(
        data = data, 
        annualIndex = annualIndex, 
        AC = AC,
        recent.years = recent.years, 
        determinand = determinand, 
        max.year = max.year,
        recent.trend = recent.trend,
        nYearFull = nYearFull, 
        firstYearFull = firstYearFull
      )
    )
    output$convergence <- 0L
    return(output)
  }

  
  # strip off early years with 'less-thans' - ensure that at least 50% of years
  # have at least one real observation - but keep at least two years of data
  # also, if have 5 or more years of positives (so that a trend is fitted), then remove all 
  # less-thans at the start of the time series
  
  data <- ctsm.remove.early.lessThans(data)
  
  
  # decide whether there are sufficient years to model data
  # appropriate type of fit depends on number of years

  nYear <- length(unique(data$year))

  if (nYear > 2) {
    
    # fit a mean level for (starters) 
    # if get bound convergence, strip off early years to see if that resolves the problem - essentially
    # seeing if an outlying year is responsible in a very crude way
    
    fitMean <- ctsm.lmm.fit(data, dfSmooth = 0, ...)
    
    while (!fitMean$convergence$bound$fixed | !fitMean$convergence$bound$random) {
        
      data <- subset(data, year > min(year))
      data <- ctsm.remove.early.lessThans(data)
      
      nYear <- length(unique(data$year))
      if (nYear == 2) break()

      fitMean <- ctsm.lmm.fit(data, dfSmooth = 0, ...)
    }      
  }  

  # calculate number of years with positives - determines type of fit:
  # nYear <= 2 none
  # nYearPos <= 4 mean (might be some years with as few as two positives)
  # nyearPos >= 5 linear or smooth
  # nYearPos >= 7, 10, 15 try smooths on 2, 3, 4 df
  
  nYearPos <- with(data, sum(tapply(censoring, year, function(x) any(x == ""))))
  

  # initialise output
  
  output <- list(data = data)
  
    
  # now do the remaining fits 
  
  if (nYear <= 2) 
    
    output$method <- "none"
  
  else {
    
    fits <- list(fitMean)
    
    if (nYearPos >= 5) 
      fits[[2]] <- ctsm.lmm.fit(
        data, dfSmooth = 1, varComp.start = with(fits[[1]], coefficients[idRandom, ]), ...)

    if (nYearPos >= 7) 
      fits[[3]] <- ctsm.lmm.fit(
        data, dfSmooth = 2, varComp.start = with(fits[[2]], coefficients[idRandom, ]), ...)
    
    if (nYearPos >= 10) 
      fits[[4]] <- ctsm.lmm.fit(
        data, dfSmooth = 3, varComp.start = with(fits[[3]], coefficients[idRandom, ]), ...)
    
    if (nYearPos >= 15) 
      fits[[5]] <- ctsm.lmm.fit(
        data, dfSmooth = 4, varComp.start = with(fits[[4]], coefficients[idRandom, ]), ...)
    
    output$anova <- do.call(
      rbind, lapply(fits, function(x) as.data.frame(x[c("twiceLogLik", "AIC", "AICc")])))

    row.names(output$anova) <- 
      c("mean", "linear", "smooth (df = 2)", "smooth (df = 3)", "smooth (df = 4)")[1:nrow(output$anova)]
  
  
    # choose best model: 
    # nYearPos <= 4 gives mean model
    # otherwise choose linear or smooth model with minimum AICc
    # if choose_model specified, take model with corresponding degrees of freedom (for use in emergency)
    
    if (missing(choose_model)) 
      bestFit <- if (nYearPos <= 4) 1 else max(2, which.min(output$anova$AICc))
    else 
      bestFit <- choose_model
    
    if (bestFit > length(fits))
      stop("model choice does not exist in ctsm.anyyear.lmm", call. = FALSE)

    fit <- ctsm.lmm.hess(fits[[bestFit]], ...)

    output <- within(output, {
      method <- if (bestFit == 1) "mean" else if (bestFit == 2) "linear" else "smooth" 
      if (method == "smooth") dfSmooth <- bestFit - 1
      pred <- fit$pred
      coefficients <- fit$coefficients
    })
  }

  
  # use coefficients to estimate the variance of the mean response observed each year
  # assuming equal number of samples and equal analytical quality
  # number of samples per year is nrow(data) / nYear
  
  sdAnalytic <- with(data, median(sdAnalytic))
  sdSample <- if (nrow(data) > nYear) output$coefficients["sdSample", "est"] else 0
  sdYear <- output$coefficients["sdYear", "est"]
  nSample <- nrow(data) / nYear
  
  sigma <- sqrt(sdYear ^ 2 + (sdSample ^ 2 + sdAnalytic ^ 2) / nSample)
  
  output$sd_components <- c(
    sd_analytic = sdAnalytic, sd_sample = sdSample, sd_year = sdYear, n_sample = nSample, 
    sd_index = sigma
  )
  
  
  # get estimated change in concentration over whole time series and in the most recent 
  # e.g. twenty years of monitoring (truncate when data missing and only compute if at least five years)
  # in that period
  # NB p value from contrast is NOT the same as from likelihood ratio test even if method = "linear"

  if (output$method %in% c("linear", "smooth")) {

    contrast.whole <- ctsm.lmm.contrast(fit, start = min(data$year), end = max(data$year))
    row.names(contrast.whole) <- "whole"
    
    start_recent <- max(max.year - recent.trend + 1, min(data$year))
    if (sum(unique(data$year) >= start_recent - 0.5) >= 5) {
      contrast.recent <- ctsm.lmm.contrast(fit, start = start_recent, end = max(data$year))
      row.names(contrast.recent) <- "recent"
      contrast.whole <- rbind(contrast.whole, contrast.recent)
    }		
    
    output$contrasts <- contrast.whole
  }
     
  # compare fitted concentration in final year against assessment concentrations
  # ACs may be log-transformed to get onto the same scale as the index
  # some determinands have good status when high, some when low
     
  if (output$method %in% c("mean", "linear", "smooth")) {

    output$reference.values <- lapply(AC, function(i) {
      ctsm.lmm.refvalue(
        fit, 
        yearID = max(data$year), 
        refvalue = switch(distribution, lognormal = log(i), normal = i), 
        lower.tail = switch(good.status, low = TRUE, high = FALSE)
      )
    })
    
    output$reference.values <- do.call("rbind", output$reference.values)
  }


  # compute power statistics (other than dtrend)
  # need to extend this to normally distributed data at some point
  
  if (distribution == "lognormal") {
    output$power <- ctsm_lmm_power(
      output, 
      target_power = power$target_power,
      target_trend = power$target_trend,
      size = power$size
    )
  }
  
  
  # construct summary output -
  
  output$summary <- data.frame(
    nyall = nYearFull, nyfit = nYear, nypos = nYearPos, 
    firstYearAll = firstYearFull, firstYearFit = min(data$year), lastyear = max(data$year), 
    p_nonlinear = NA, p_linear = NA, p_overall = NA, pltrend = NA, ltrend = NA, prtrend = NA, 
    rtrend = NA, dtrend = NA, meanLY = NA, clLY = NA)
  
  
  output$summary <- within(output$summary, {
    
    if (output$method == "smooth") {
      
      p_nonlinear <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["linear", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })

      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })

      p_overall <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })
      
    }
      
    if (output$method == "linear") {
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- p_linear
      
    }
    
    if (output$method %in% c("linear", "smooth")) {
      
      # pltrend
      # method = "linear" use p_linear (from likelihood ratio test)  
      # method = "smooth" use p from the Wald test in contrasts 
      # for linear model, likelihood ratio test is a better test (fewer 
      #   approximations) than the Wald test
      # for smooth model, would be better to go into profile likelihood 
      #   territory (future enhancement) 
      
      # prtrend
      # same approach; however p_linear could be misleading when the years at 
      #   the end of the time series are all censored values and a flat model is 
      #   fitted; the estimate of rtrend is shrunk to reflect this, but p_linear 
      #   might be misleadingly significant; something to think about in the 
      #   future
      # however, there is a pathological case when all the fitted values in the 
      #   recent period have the same value; rtrend is zero, and yet can still be 
      #   significant based on p_linear even though there are no data to support 
      #   this; in this case use p from the Wald test (which is unity)
      
      if (output$method == "linear") {
        pltrend <- p_linear
      } else {
        pltrend <- output$contrasts["whole", "p"]
      }
      
      ltrend <- with(output$contrasts["whole", ], estimate / (end - start))
      
      if ("recent" %in% row.names(output$contrasts)) {
        
        if (
          output$method == "linear" & 
          max(data$year[data$censoring %in% ""]) > start_recent
        ) {
          prtrend <- p_linear
        } else {
          prtrend <- output$contrasts["recent", "p"]
        }          

        rtrend <- with(output$contrasts["recent", ], estimate / (end - start))
      }
    }
                              
    # if parametric model cannot be fitted, use maximum index in last two monitoring years 
    # (if low values are good) or minimum index (if low values are bad) for crude extra data

    if (output$method == "none") 
      meanLY <- local({
        index <- tail(annualIndex$index, nYear)
        switch(
          good.status, 
          low = max(index), 
          high = min(index)
        )
      })
    else {
      meanLY <- tail(output$pred$fit, 1)
      clLY <- switch(
        good.status, 
        low = tail(output$pred$ci.upper, 1), 
        high = tail(output$pred$ci.lower, 1)
      )
      dtrend <- ctsm_dtrend(
        1:10, 
        sigma, 
        alpha = power$size / 100,
        power = power$target_power / 100
      )
    }
                             
                             
    # backtransform for log-normal distribution (and to give percentage trends)
    
    if (distribution == "lognormal") {
      ltrend <- ltrend * 100
      rtrend <- rtrend * 100
      dtrend <- dtrend * 100
      meanLY <- exp(meanLY)
      clLY <- exp(clLY)
    }
    
  })  
  
  if (!is.null(AC)) {
    output$summary <- data.frame(output$summary, do.call(cbind, lapply(names(AC), function(i) {
      
      value <- AC[i]
      diff <- with(output, if (method == "none") summary$meanLY - value else summary$clLY - value)
      
      # estimate number of years until meanLY reaches target - based on rtrend
      # might be already there but cl is too high
      
      maxYear <- max(data$year)
      bigYear <- 3000
      
      tillTarget <- with(output$summary, {
        
        if (good.status == "low") {
        
          if (is.na(value) || (meanLY >= value & is.na(rtrend)))
            NA
          else if (meanLY < value) 
            maxYear
          else if (rtrend >= 0)
            bigYear
          else {
            wk <- switch(
              distribution, 
              lognormal = 100 * (log(value) - log(meanLY)) / rtrend,
              (value - meanLY) / rtrend)
            wk <- round(wk + maxYear)
            min(wk, bigYear)
          }
          
        } else {
          
          if (is.na(value) || (meanLY <= value & is.na(rtrend)))
            NA
          else if (meanLY > value) 
            maxYear
          else if (rtrend <= 0)
            bigYear
          else {
            wk <- switch(
              distribution, 
              lognormal = 100 * (log(value) - log(meanLY)) / rtrend,
              (value - meanLY) / rtrend)
            wk <- round(wk + maxYear)
            min(wk, bigYear)
          }

        }
      })
        
      out <- data.frame(value, diff, tillTarget, below.result[[i]])
      names(out) <- paste0(i, c("", "diff", "achieved", "below"))
      out
    })))
  }
  
  output$summary <- within(output$summary, {
  
    # and round for ease of interpretation
    
    p_nonlinear <- round(p_nonlinear, 4)
    p_linear <- round(p_linear, 4)
    p_overall <- round(p_overall, 4)
    pltrend <- round(pltrend, 4)
    prtrend <- round(prtrend, 4)
    
    ltrend <- round(ltrend, 1)
    rtrend <- round(rtrend, 1)
    dtrend <- round(dtrend, 1)
  })
  
    
  rownames(output$summary) <- NULL
  output
}



ctsm.test.below <- function(year, index, value, min.year, below = TRUE) {
  
  if (is.na(value)) return (NA)
  
  ok <- year >= min.year
  index <- index[ok]
  n.year <- length(year[ok])
  
  if (n.year < 5) return (NA)     # will never get a significant result

  # now focus on most recent 5 years - used to use all available data but, 
  # with year-skipping stategies, ended up with time series being above the 
  # e.g. EAC due to a single observation more than ten years before the current 
  # monitoring year
  
  n.year <- 5
  index <- tail(index, n.year)

  n.bad <- if (below) sum(index >= value) else sum(index <= value)
  sig <- pbinom(n.bad, n.year, 0.5) < 0.05
  
  if (below)
    if (sig) "below" else "above"
  else
    if (sig) "above" else "below"
}



ctsm.remove.early.lessThans <- function(data) {

  check <- with(data, length(unique(year)) > 2 & any(censoring %in% c("<", "D", "Q")))
  if (!check) return(data)

  # identify years with at least one real
  
  wk <- with(data, tapply(censoring, year, function(x) any(x == "")))

  # if 5 or more years with reals, then exclude all completely less than years at the start of the time series
  
  if (sum(wk) >= 5) {
    wk <- wk[which.max(wk):length(wk)]
    data <- data[as.character(data$year) %in% names(wk), ]
  } 
  
  # reverse to find out, going backwards, which years give a span of data
  # in which 50% of years have at least one real
  
  ok <- round(100 * cumsum(rev(wk)) / seq(1, length(wk))) >= 50
  
  # ensure always include first two years of data (for an index) and reverse
  # back to get in correct order
  
  ok[1:2] <- TRUE
  ok <- rev(ok)
  
  # find which years satisfy constraint and simplify data set
  
  ok <- ok[which.max(ok):length(ok)]
  data[as.character(data$year) %in% names(ok), ]
}


ctsm.lmm.refvalue <- function(ctsm.ob, yearID, refvalue, ...) {

  ok <- ctsm.ob$pred$year %in% yearID 
  if (!any(ok)) stop("requested year not found in predicted values")
  pred <- ctsm.ob$pred[ok, ]
  
  fit <- pred$fit
  difference <- refvalue - fit
  se <- pred$se
  t.stat <- difference / se
  p.value = 1 - pt(t.stat, ctsm.ob$dfResid, ...)

  data.frame(year = yearID, fit, refvalue, difference, se, p = p.value)	
}  
  

ctsm.lmm.contrast <- function(ctsm.ob, start, end) {

  # error trapping
  
  if (length(start) > 1 | length(end) > 1) 
    stop('only a single contrast is allowed: start and end must both be scalars')  

  pos <- match(c(start, end), ctsm.ob$pred$year)
  if (any(is.na(pos))) stop('start or end year not found in predicted data')

  wk <- c(-1, 1)
  contrast <- t(wk) %*% ctsm.ob$pred$fit[pos]
  
  wk <- t(wk) %*% ctsm.ob$Xpred[pos, ]
  se.contrast <- sqrt(wk %*% ctsm.ob$vcov %*% t(wk))

  # catch pathological case where contrast = 0 and se.contrast = 0
  # this can happen if all the data between start and end are censored, so
  # a 'flat' model is fitted
  
  if (dplyr::near(contrast, 0L) & dplyr::near(se.contrast, 0L)) {
    p.contrast <- 1
  } else {
    t.stat <- contrast / se.contrast
    p.contrast <- 1 - pf(t.stat^2, 1, ctsm.ob$dfResid)
  }
  
  data.frame(start, end, estimate = contrast, se = se.contrast, p = p.contrast)
}


#' Check whether the assessments have converged
#'
#' @description Checks whether the assessments in a HARSAT assessment object
#'   have converged. Currently only does detailed checks for models with normal
#'   or lognormal errors (all chemical timeseries and some biological effects).
#'   Here it checks whether:
#'
#' * fixed effect estimates are away from their bounds
#' * random effect estimates are lower than their upper bounds; they can of
#'   course be equal to zero
#' * standard errors are present for model predictions and fixed effects
#'   estimates
#' * standard errors of the fixed effects estimates are realistic (very small standard errors indicate problems with
#'   the numerical differencing used to compute the standard errors)
#'
#'
#' @param assessment_ob A HARSAT assessment object resulting from a call to
#'   ctsm_assessment
#' @param save_result Saves the identifiers of the timeseries that have not
#'   converged; defaults to FALSE. When an assessment has been done in stages, the output also identifies those 
#'   timeseries that have not yet been assessed
#'
#' @returns A list of two character vectors: 
#' 
#' * not_converged identifies timeseries that have not converged
#' * not_assessed identifies timeseries that have not been assessed
#'
#' @export
check_assessment <- function(assessment_ob, save_result = FALSE) {

  assessment <- assessment_ob$assessment

  n <- length(assessment)
  
  
  # deal with situation where not all time series have been assessed 
    
  assessed <- sapply(assessment, function(x) !is.null(x))
  
  if (!any(assessed)) {
    warning(
      "No timeseries have been assessed\n", 
      call. = FALSE, immediate. = TRUE
    )
    return(invisible())
  }
  
  if (!all(assessed)) {
    warning(
      "Not all timeseries have been assessed; to see which ones, set argument\n",
      "'save_result' to TRUE\n", 
      call. = FALSE, immediate. = TRUE
    )
  }
  
  
  # restrict assessment to those time series that have been assessed
  
  assessment <- assessment[assessed]
  
  converged <- sapply(assessment, function(x) x$convergence == 0L)

  
  # identify time series that have not been assessed or have not converged  
  
  result <- list(
    not_converged = names(converged)[!converged],
    not_assessed = names(assessed[!assessed])
  )
  
  if (all(converged)) {
    cat("All assessment models have converged\n")
  } else {
    cat(
      "The following assessment models have not converged:\n",
      paste(result$not_converged, collapse = "\n"), 
      sep = ""
    )
  } 
  
  if (save_result) result else invisible()
}


#' Checks convergence of an assessment model
#'
#' @description Utilty function for use within assess_lmm. Checks whether a
#'   model has converged. Currently only checks assessments with normal (or
#'   lognormal) errors, where it considers whether:
#'
#' * fixed effect estimates are away from their bounds
#' * random effect estimates are lower than their upper bounds; they can of
#'   course be equal to zero
#' * standard errors are present for model predictions and fixed effects
#'   estimates
#' * standard errors of the fixed effects estimates are not unrealistically small;
#'   the tolerance is chosen to be much smaller than seen in typical OSPAR
#'   assessments which have converged
#'
#'   Model fits based on other distributions are assumed to have converged. Some
#'   checking is needed here in future
#'
#' @param assessment An assessment from assess_lmm
#' @param coeff_se_tol The tolerance for checking whether standard errors on the
#'   fixed effects estimates are unrealistically small; defaults to 0.001
#'
#' @returns An integer: 0 indicates convergence, 1 indicates an issue
check_convergence_lmm <- function(assessment, coeff_se_tol = 0.001) {

  # error in model fit
  
  if ("error" %in% names(assessment)) {
    return(1L)
  }
  
  # nothing to check if no model fitting (determined by presence of pred component)
  
  if (!"pred" %in% names(assessment)) {
    return(0L)
  }
  
  # all standard errors from pred component should be present
  
  if (any(is.na(assessment$pred$se))) {
    return(1L)
  }
  
  # fixed effect coefficients: standard errors should be present and not 
  #   ridiculously small
  # fixed effect coefficients should not be on bounds
  # random effect coefficients should not be on upper bounds
  
  coeff <- assessment$coefficients
  
  random_id <- grepl("sd", row.names(coeff))
  fixed <- coeff[!random_id, ]
  random <- coeff[random_id, ]
  
  if (any(is.na(fixed$se)) || min(fixed$se) < coeff_se_tol) {
    return(1L)
  }
  
  if (any(fixed$onBound)) {
    return(1L)
  }
  
  if (any(random$onBound & random$est > 0.0001)) {
    return(1L)
  }
  
  0L
}


# Power functions ----

ctsm_lmm_power <- function(assessment, target_power = 80, target_trend = 10, size = 5) {
  
  # intialise output
  
  id <-c(
    "dtrend_obs", "dtrend_seq", "dtrend_ten", 
    "nyear_seq", 
    "power_obs", "power_seq", "power_ten" 
  )
  
  out <- vector("list", 7) 
  names(out) <- id
  out[id] <- NA
  
  
  # get key data, and return if too few years to compute power
  
  year <- unique(assessment$data$year)

  sd <- assessment$sd_components["sd_index"]
  
  method <- assessment$method 
  
  if (method == "none") {
    return(out)
  }
  
  
  # detectable trend (on log scale) of 
  # 1 current time series over observed time span
  # 2 current time series with no gaps (e.g. annual monitoring from min_year to max_year)
  # 3 in ten years of sequential monitoring
  
  target_power <- target_power / 100
  
  if (method %in% c("linear", "smooth")) {
    out["dtrend_obs"] <- ctsm_dtrend(year, sd, power = target_power)
    out["dtrend_seq"] <- ctsm_dtrend(min(year):max(year), sd, power = target_power)
  }
  
  out["dtrend_ten"] <- ctsm_dtrend(1:10, sd, power = target_power)
  
  # back-transform to percentage annual (positive) change
  
  id <- c("dtrend_obs", "dtrend_seq", "dtrend_ten")
  
  out[id] <- lapply(out[id], function(y) round(100 * (exp(y) - 1), 1))
  
  
  # number of sequential years to detect the specified % change
  
  target_trend <- target_trend / 100
  
  out["nyear_seq"] <- ctsm_dyear(log(1 + target_trend), sd, power = target_power)
  
  
  # power to detect the specified % change with same options as dtrend
  
  if (method %in% c("linear", "smooth")) {
    out["power_obs"] <- ctsm_dpower(log(1 + target_trend), year, sd)
    out["power_seq"] <- ctsm_dpower(log(1 + target_trend), min(year):max(year), sd)
  }
  
  out["power_ten"] <- ctsm_dpower(log(1 + target_trend), 1:10, sd)
  
  # turn into percentages
  
  id <- c("power_obs", "power_seq", "power_ten") 

  out[id] <- lapply(out[id], function(y) round(100 * y))
  
  out
}



ctsm_dpower <- function(
    q, year, sigma, alpha = 0.05, sigma_type = c("index", "slope"), 
    alternative = c("two.sided", "less", "greater")) {
  
  # power of (log-)linear regression
  
  # q is the annual change (of the log-linear trend)
  # year is the vector of years that were sampled
  # sigma is the standard deviation:
  #   type = "index"; the standard deviation of the yearly index which is useful in 
  #     stylised situations when the index is measured with constant variance
  #   type = "slope"; the standard deviation of the estimator of beta, which allows for year-
  #     specific variances on the yearly indices; this will often be taken from a linear 
  #     regression anova table
  # alpha is the size of the test
  
  sigma_type <- match.arg(sigma_type)
  alternative = match.arg(alternative)
  
  stopifnot(
    length(year) >= 3, 
    !duplicated(year),
    sigma > 0, 
    alpha > 0,
    alpha < 1
  )
  
  df <- length(year) - 2
  
  delta <- switch(
    sigma_type,
    index = sqrt(sum((year - mean(year)) ^ 2)),
    slope = 1
  )
  
  delta <- q * delta / sigma
  
  switch(
    alternative, 
    two.sided = {     
      crit_val <- qf(1 - alpha, 1, df)
      delta <- delta ^ 2
      pf(crit_val, 1, df, delta, lower.tail = FALSE)
    },
    less = {
      crit_val <- qt(alpha, df)
      pt(crit_val, df, delta, lower.tail = TRUE)
    },
    greater = {
      crit_val <- qt(1 - alpha, df)
      pt(crit_val, df, delta, lower.tail = FALSE)
    }
  )
}


ctsm_dtrend <- function(
    year, sigma, alpha = 0.05, power = 0.8, sigma_type = c("index", "slope"), 
    alternative = c("two.sided", "less", "greater")) {
  
  sigma_type <- match.arg(sigma_type)
  alternative <- match.arg(alternative)
  
  # calculates lowest annual change detectable in time series 
  
  if (alternative == "less") {
    lower <- -90
    upper <- 0
    extendInt <- "no"
  } else {
    lower <- 0
    upper <- 10
    extendInt <- "upX"
  }
  
  uniroot(
    function(q) 
      ctsm_dpower(q, year, sigma, alpha, sigma_type, alternative) - power, 
    lower = lower, upper = upper, extendInt = extendInt
  )$root
}


ctsm_dyear <- function(
    q, sigma, alpha = 0.05, power = 0.8, sigma_type = c("index", "slope"), 
    alternative = c("two.sided", "less", "greater")) {
  
  sigma_type <- match.arg(sigma_type)
  alternative <- match.arg(alternative)
  
  # calculates number of years required to detect trend
  
  if ((alternative == "greater" & q < 0) | (alternative == "less" & q > 0))
    stop("sign of q incompatible with test alternative")
  
  achieved_power <- 0
  n_year <- 2
  while (achieved_power < power) {
    n_year <- n_year + 1
    achieved_power <- ctsm_dpower(q, 1:n_year, sigma, alpha, sigma_type, alternative)
  }
  
  n_year
}





# Other distributions ----

assess_survival <- function(
  data, annualIndex, AC, recent.years, determinand, max.year, recent.trend, 
  nYearFull, firstYearFull) {

  # silence non-standard evaluation warnings
  .data <- est <- lcl <- ucl <- p <- se <- NULL
  info <- NULL

  # assess survival data, often expressed as interval data 
  
  # check valid determinands 
  
  if (! determinand %in% c("LP", "NRR", "SURVT")) {
      stop("not yet coded for determinand: ", determinand)
  }
  
  
  # initialise output
  
  output <- list(data = data)

 
  # check suitable values and construct censoring information
  # event = 2 is left censored
  # event = 3 is interval censored 
  # time = lower bound
  # time2 = upper bound (or Inf if event = 2)
  # response can be either lower or upper bound!
  
  valid_values <- switch(
    determinand, 
    LP = c(0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50), 
    NRR = c(0, 15, 30, 60, 90, 120, 150, 180), 
    SURVT = seq(1, max(data$response))
  )  
  
  if (determinand %in% c("NRR", "LP")) {

    # response is the lower bound 
    
    data <- dplyr::mutate(
      data, 
      event = dplyr::if_else(.data$response == max(valid_values), 2, 3),
      time = .data$response,
      time2 = factor(
        .data$time, 
        levels = valid_values, 
        labels = as.character(c(valid_values[-1], Inf))
      ),
      time2 = as.numeric(as.character(.data$time2))
    )
    
  } else {
    
    # response is the upper bound 
    
    n_valid = length(valid_values)
    
    data <- dplyr::mutate(
      data, 
      event = 3,
      time2 = .data$response,
      time = factor(
        .data$time2, 
        levels = valid_values, 
        labels = as.character(c(0, valid_values[-n_valid]))
      ),
      time = as.numeric(as.character(.data$time))
    )
    
  }
    
  
  # lower bounds of zero are not compatible with survival distribution, so 
  #   replace with one tenth of upper bound
  
  data <- dplyr::mutate(
    data, 
    time = dplyr::if_else(.data$time == 0, .data$time2 / 10, .data$time)
  )
  
  
  # define survival distribution
  
  surv_dist <- switch(
    determinand, 
    NRR = "gamma", 
    SURVT = "gamma",
    NA
  )
  
  
  # get number of years and adjust year variable for stability
  # nYearPos is calculated for consistency with other responses
  
  nYear <- nYearPos <- length(unique(data$year))
  
  data$year_adj <- data$year - min(recent.years)
    
  
  # establish other info
  
  good_status <- ctsm_get_info(info$determinand, determinand, "good_status")
  

  # type of fit depends on number of years:
  # nYear <= 2 none
  # nYear <= 4 mean 
  # nyear >= 5 linear or smooth
  # nYear >= 7, 10, 15 try smooths on 2, 3, 4 df
  
  # have only currently coded for mean and linear - look at ctsm.anyyear.lmm for 
  # extensions to smoothers

  if (determinand %in% c("NRR", "SURVT") & nYear >= 7) {
    stop("time series too long: need to include code for smoothers")
  } 
    
  if (determinand %in% "LP" & nYear >= 3) {
    stop("first time parametric fit possible: need to select survival distribution")
  } 

      
  # now do fits 
  
  if (nYear <= 2) 
    
    output$method <- "none"
  
  else {
    
    fits <- list(mean = NULL)
    
    check_fit <- function(model_fit) {
      if (model_fit$opt$convergence != 0) {
        stop("non-convergence: investigate")
      } else {
        return(invisible)
      }
    }


    # mean model
    
    fits$mean <- flexsurv::flexsurvreg(
      Surv(time, time2, type = "interval2") ~ 1,
      dist = surv_dist, 
      data = data
    )
    
    check_fit(fits$mean)

    
    # linear model
    
    if (nYear >= 5) {
      fits$linear <- update(fits$mean, .~. + year_adj)
      check_fit(fits$linear)
    }  
      
    
    # full model for overdispersion test (need to incorporate random effects 
    # at some point)

    fits$full <- update(fits$mean, .~. + as.factor(year))
    check_fit(fits$full)


    # get basic anova 
    
    output$anova <- data.frame(
      p = sapply(fits, "[[", "npars"),
      twiceLogLik = sapply(fits, function(i) 2 * logLik(i)), 
      AIC = sapply(fits, AIC)
    )


    # choose best model based on AIC
    # options are currently 1 vs 3 if nYear <= 4, or 2 vs 3 if nYear >= 5

    bestFit <- if (nYear <= 4) 1 else 2
    
    if (output$anova["full", "AIC"] < output$anova[bestFit, "AIC"]) {
      cat("  warning: over-dispersed - consider incorporating random year effect\n")
      output$convergence <- "over-dispersed"
    } else {
      output$convergence <- 0
    }
      
    fit <- fits[[bestFit]]
    
    output$method <- 
      if (bestFit == 1) {
        "mean" 
      } else if (bestFit == 2) { 
        "linear" 
      } else {
        "smooth"
      }
      
    if (output$method == "smooth") {
      output$dfSmooth <- bestFit - 1
    }

    
    # predicted mean survival times with pointwise two-sided 90% confidence 
    #   limits
    # censoring means that the fitted values tend to be higher (NRR, LP) or 
    #   lower (SURVT) than the annual indices
    # store confidence level in output because needed for AC comparisons
    # confidence intervals and standard errors are simulation based, so need to 
    #   set seed for repeatability (hard-wired at present); standard predict 
    #   function doesn't allow user to increase B (number of simulations) so 
    #   use summary instead
        
    new_data <- data.frame(year = seq(min(data$year), max(data$year)))
    
    new_data$year_adj <- new_data$year - min(recent.years)

    output$conf_level_predictions <- 0.90
    
    output$seed <- 220526
    
    set.seed(output$seed)
    
    cat("  warning: random seed hard-wired in code\n")
    
    pred <- summary(
      fit, 
      new_data, 
      type = "mean",
      se = TRUE, 
      ci = TRUE, 
      cl = output$conf_level_predictions, 
      B = 100000, 
      tidy = TRUE
    )

    pred <- dplyr::rename(
      pred, 
      fit = est,
      ci.lower = lcl,
      ci.upper = ucl
    )
    
    pred <- pred[c("fit", "se", "ci.lower", "ci.upper")]
    
    output$pred <- cbind(
      new_data["year"], 
      pred[c("fit", "se", "ci.lower", "ci.upper")]
    )

    row.names(output$pred) <- NULL
    
    
    # model coefficients
        
    coefficients <- fit$res
    
    coefficients <- as.data.frame.matrix(coefficients)
    
    coefficients <- tibble::rownames_to_column(coefficients, ".term")
    
    coefficients <- dplyr::mutate(coefficients, t = est / se)

    if (nYear >= 5) {
      coefficients <- dplyr::mutate(
        coefficients, 
        p = 2 * pnorm(abs(t), lower.tail = FALSE),
        p = round(p, 4)
      )
    }

    coefficients <- tibble::column_to_rownames(coefficients, ".term")
    
    output$coefficients <- coefficients
  }
  
  
  # get estimated change in log mean survival over whole time series and in the 
  # most recent # e.g. twenty years of monitoring (truncate when data missing 
  # and only compute if at least five years in that period)
  # NB p value from contrast is NOT the same as from likelihood ratio test even 
  # if method = "linear"
  
  if (output$method %in% c("linear", "smooth")) {
    
    if (output$method == "smooth") {
      stop("need to update code")
    }
    
    contrast.whole <- assess_survival_contrast(
      output, 
      start = min(data$year), 
      end = max(data$year)
    )
    row.names(contrast.whole) <- "whole"
    
    start.year <- max(max.year - recent.trend + 1, min(data$year))
    if (sum(unique(data$year) >= start.year - 0.5) >= 5) {
      contrast.recent <- assess_survival_contrast(
        output, 
        start = start.year, 
        end = max(data$year)
      )
      row.names(contrast.recent) <- "recent"
      contrast.whole <- rbind(contrast.whole, contrast.recent)
    }		
    
    output$contrasts <- contrast.whole
  }

  
  # compare mean survival time in final year to assessment criteria
  # determinands typically have good status when high
  # results are presented on the log-scale because the confidence intervals
  #   are more symmetric; however there are still some slight inconsistencies 
  #   between the p-values produced here and the confidence limits from the 
  #   predict function
  # could get consistency by computing confidence intervals for different 
  #   coverages and finding the coverage that equates to the AC

  if (output$method %in% c("mean", "linear", "smooth")) {

    output$reference.values <- lapply(AC, function(i) {
      assess_survival_refvalue(
        output, 
        year_id = max(data$year), 
        refvalue = i,
        good_status = good_status
      )
    })
      
    output$reference.values <- do.call("rbind", output$reference.values)
  }
      
  
  # construct summary output -
  
  output$summary <- data.frame(
    nyall = nYearFull, nyfit = nYear, nypos = nYearPos, 
    firstYearAll = firstYearFull, firstYearFit = min(data$year), lastyear = max(data$year), 
    p_nonlinear = NA, p_linear = NA, p_overall = NA, pltrend = NA, ltrend = NA, prtrend = NA, 
    rtrend = NA, dtrend = NA, meanLY = NA, clLY = NA)

  
  output$summary <- within(output$summary, {
    
    if (output$method == "smooth") {
      
      p_nonlinear <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["linear", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })

    }
    
    if (output$method == "linear") {
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- p_linear
      
    }
    
    if (output$method %in% c("linear", "smooth")) {
      
      # for linear trend and recent trend, use pltrend (from likelihood ratio test) if 
      # method = "linear", because a better test 
      # really need to go into profile likelihood territory here!
      
      pltrend <- if (output$method == "linear") {
        p_linear
      } else {
        output$contrasts["whole", "p"]
      }
      
      ltrend <- with(output$contrasts["whole", ], estimate / (end - start))
      
      if ("recent" %in% row.names(output$contrasts)) {
        prtrend <- if (output$method == "linear") {
          p_linear
        } else {
          with(output$contrasts["recent", ], p)
        }
        rtrend <- with(output$contrasts["recent", ], estimate / (end - start))
      }
    }
    
    # if parametric model cannot be fitted, use maximum index in last two monitoring years 
    # (if low values are good) or minimum index (if low values are bad) for crude extra data
    
    if (output$method == "none") 
      meanLY <- local({
        index <- tail(annualIndex$index, nYear)
        switch(
          good_status, 
          low = max(index), 
          high = min(index)
        )
      })
    else {
      meanLY <- tail(output$pred$fit, 1)
      clLY <- switch(
        good_status, 
        low = tail(output$pred$ci.upper, 1), 
        high = tail(output$pred$ci.lower, 1)
      )
    }
    
    
    # turn trends into 'percentage trends'
    
    ltrend <- ltrend * 100
    rtrend <- rtrend * 100
  })  
  
  if (!is.null(AC)) {
    output$summary <- data.frame(output$summary, do.call(cbind, lapply(names(AC), function(i) {
      
      value <- AC[i]
      diff <- with(output, if (method == "none") summary$meanLY - value else summary$clLY - value)
      
      # estimate number of years until meanLY reaches target - based on rtrend
      # might be already there but cl is too high
      
      maxYear <- max(data$year)
      bigYear <- 3000
      
      tillTarget <- with(output$summary, {
        
        if (is.na(value) || (meanLY <= value & is.na(rtrend)))
          NA
        else if (meanLY >= value) 
          maxYear
        else if (rtrend <= 0)
          bigYear
        else {
          wk <- 100 * (log(value) - log(meanLY)) / rtrend
          wk <- round(wk + maxYear)
          min(wk, bigYear)
        }
      })
      
      out <- data.frame(value, diff, tillTarget, below.result = NA)
      names(out) <- paste0(i, c("", "diff", "achieved", "below"))
      out
    })))
  }
  
  output$summary <- within(output$summary, {
    
    # and round for ease of interpretation
    
    p_nonlinear <- round(p_nonlinear, 4)
    p_linear <- round(p_linear, 4)
    p_overall <- round(p_overall, 4)
    pltrend <- round(pltrend, 4)
    prtrend <- round(prtrend, 4)
    
    ltrend <- round(ltrend, 1)
    rtrend <- round(rtrend, 1)
    dtrend <- round(dtrend, 1)
  })
  
  
  rownames(output$summary) <- NULL
  output
}



assess_survival_contrast <- function(ctsm.ob, start, end) {

  # silence non-standard evaluation warnings
  fit <- NULL

  # based on ctsm.lmm.contrast - should be able to make it almost identical but
  # first need to get variance covariance matrix of fitted values
  
  # error trapping
  
  if (length(start) > 1 | length(end) > 1) 
    stop('only a single contrast is allowed: start and end must both be scalars')  
  
  pos <- match(c(start, end), ctsm.ob$pred$year)
  if (any(is.na(pos))) stop('start or end year not found in predicted data')
  
  # take logs of fitted values to get onto linear scale
  ctsm.ob$pred <- dplyr::mutate(ctsm.ob$pred, fit = log(fit))
  
  wk <- c(-1, 1)
  contrast <- t(wk) %*% ctsm.ob$pred$fit[pos]
  
  # until vcov is available, use standard error from model coefficient
  
  # wk <- t(wk) %*% ctsm.ob$Xpred[pos, ]
  # se.contrast <- sqrt(wk %*% ctsm.ob$vcov %*% t(wk))
  se.contrast <- ctsm.ob$coefficients["year_adj", "se"] * (end - start)
  
  t.stat <- contrast / se.contrast
  p.contrast <- 1 - pf(t.stat^2, 1, Inf)
  data.frame(start, end, estimate = contrast, se = se.contrast, p = p.contrast)
}



assess_survival_refvalue <- function(
  ctsm_ob, year_id, refvalue, good_status, ...) {

  ok <- ctsm_ob$pred$year %in% year_id 
  if (!any(ok)) {
    stop("requested year not found in predicted values")
  }
  pred <- ctsm_ob$pred[ok, ]
  
  # transform to the log-scale as confidence intervals more symmetric
  
  refvalue <- log(refvalue)
  
  id <- c("fit", "ci.lower", "ci.upper")
  pred[id] <- lapply(pred[id], log)
  
  fit <- pred$fit
  difference <- refvalue - fit
  
  # approximate standard error from the log-transformed confidence limits
  # base this on difference between lower confidence limit (if good status is
  #   high) or upper confidence limit (if low) so that we get consistency between
  #   confidence limits and p-values
  # conf_level of predictions are two-sided, but need a one-sided test here

  alpha <- (1 - ctsm_ob$conf_level_predictions) / 2
  
  if (good_status == "high") {
    ci = pred$ci.lower 
    lower_tail = TRUE 
  } else {
    ci = pred$ci.upper 
    lower_tail = FALSE 
  }
      
  se <- abs(ci - fit) / qnorm(1 - alpha)
  
  t_stat <- difference / se
  p_value = pnorm(t_stat, lower.tail = lower_tail)
  
  data.frame(year = year_id, fit, refvalue, difference, se, p = p_value)	
}




assess_beta <- function(
  data, annualIndex, AC, recent.years, determinand, max.year, recent.trend, 
  nYearFull, firstYearFull) {
  
  # silence non-standard evaluation warnings
  info <- weight <- NULL

  # percentage data that are not based on counts (e.g. comet assay) 
  
  # check valid determinands 
  
  if (! determinand %in% "%DNATAIL") {
    stop("not yet coded for determinand: ", determinand)
  }
  

  # initialise output
  
  output <- list(data = data)
  
  
  # check all values are valid
  # response currently expressed as a percentage
  
  data$response <- data$response / 100
  
  if (!all(data$response > 0 & data$response < 1)) {
    stop("invalid values for beta distribution data")
  }
  

  # set up weights - e.g. for comet assay these are the number of individuals
  # specified in CMT-QC-NR
  
  if (!("weight" %in% names(data))) {
    warning("  warning: weights not specified - assuming equal weights\n")
    data$weight <- 1
  }
  

  # get number of years and adjust year variable for stability
  # nYearPos is calculated for consistency with other responses
  
  nYear <- nYearPos <- length(unique(data$year))
  
  data$year_adj <- data$year - min(recent.years)
  
  data$year_fac <- factor(data$year)
  
  # establish other info
  
  good_status <- ctsm_get_info(info$determinand, determinand, "good_status")
  
  
  # type of fit depends on number of years:
  # nYear <= 2 none
  # nYear <= 4 mean 
  # nyear >= 5 linear or smooth
  # nYear >= 7, 10, 15 try smooths on 2, 3, 4 df
  
  # have only currently coded for mean and linear - look at ctsm.anyyear.lmm for 
  # extensions to smoothers
  
  if (nYear >= 7) {
    stop("time series too long: need to include code for smoothers")
  } 
  

  # do fits 
  
  if (nYear <= 2) 
    
    output$method <- "none"
  
  else {
    
    fits <- list(mean = NULL)
    
    # mean model
    
    fits$mean <- mgcv::gam(
      response ~ 1 + s(year_fac, bs = "re"), 
      weights = weight, 
      data = data, 
      family = "betar",
      method = "ML"
    )
    

    # linear model
    
    if (nYear >= 5) {
      fits$linear <- update(fits$mean, .~. + year_adj)
    }  
    
    
    # get basic anova 
    
    output$anova <- data.frame(
      p = sapply(fits, function(i) {
        edf <- i$edf
        id <- !grepl("year_fac", names(edf))
        edf <- edf[id]
        sum(edf) + 2
      }),
      twiceLogLik = sapply(fits, function(i) - 2 * i$gcv.ubre) 
    )
    
    output$anova$AIC <- - output$anova$twiceLogLik + 2 * output$anova$p
    
    
    # choose best model based on AIC

    bestFit <- if (nYear <= 4) 1 else max(2, which.min(output$anova$AIC))

    fit <- fits[[bestFit]]
    
    output$method <- 
      if (bestFit == 1) {
        "mean" 
      } else if (bestFit == 2) { 
        "linear" 
      } else {
        "smooth"
      }
    
    if (output$method == "smooth") {
      output$dfSmooth <- bestFit - 1
    }
    
    output$dfResid <- nYear - bestFit
    

    # predicted values with pointwise two-sided 90% confidence limits

    new_data <- data.frame(year = seq(min(data$year), max(data$year)))
    
    new_data$year_adj <- new_data$year - min(recent.years)
    
    new_data$year_fac <- factor(new_data$year)
    
    pred <- predict(fit, new_data, type = "lpmatrix")
    
    id <- !grepl("year_fac", dimnames(pred)[[2]])
    Xpred <- pred[, id, drop = FALSE]
    
    vcov <- vcov(fit)[id, id, drop = FALSE]
    coefficients <- coefficients(fit)[id]
    
    pred <- data.frame(
      year = new_data$year, 
      fit = c(Xpred %*% coefficients)
    )
    
    cov_fit <- Xpred %*% vcov %*% t(Xpred)
    pred$se <- sqrt(diag(cov_fit))

    output$conf_level_predictions <- 0.90
    
    alpha <- (1 - output$conf_level_predictions) / 2

    pred$ci.lower <- pred$fit + pred$se * qt(alpha, output$dfResid)
    pred$ci.upper <- pred$fit + pred$se * qt(1 - alpha, output$dfResid)
    
    output$pred <- pred
    
    row.names(output$pred) <- NULL
    
    
    # model coefficients
    coefficients <- data.frame(
      est = coefficients, 
      se = sqrt(diag(vcov))
    )
    
    coefficients$t <- coefficients$est / coefficients$se
    
    coefficients$p <- 2 * pt(
      abs(coefficients$t), 
      df = output$dfResid, 
      lower.tail = FALSE
    )
     
    coefficients$p = round(coefficients$p, 4)

    output$coefficients <- coefficients
  }
  
  
  # get estimated change in logit value over whole time series and in the 
  # most recent # e.g. twenty years of monitoring (truncate when data missing 
  # and only compute if at least five years in that period)
  # NB p value from contrast is NOT the same as from likelihood ratio test even 
  # if method = "linear"

  if (output$method %in% c("linear", "smooth")) {
    
    fit_info <- list(
      pred = output$pred, 
      dfResid = output$dfResid, 
      Xpred = Xpred,
      vcov = vcov
    )
    
    contrast.whole <- ctsm.lmm.contrast(fit_info, start = min(data$year), end = max(data$year))
    row.names(contrast.whole) <- "whole"
    
    start.year <- max(max.year - recent.trend + 1, min(data$year))
    if (sum(unique(data$year) >= start.year - 0.5) >= 5) {
      contrast.recent <- ctsm.lmm.contrast(fit_info, start = start.year, end = max(data$year))
      row.names(contrast.recent) <- "recent"
      contrast.whole <- rbind(contrast.whole, contrast.recent)
    }		
    
    output$contrasts <- contrast.whole
  }
        
  
  # compare mean value in final year to assessment criteria
  # results are presented on the logistic scale

  if (output$method %in% c("mean", "linear", "smooth")) {
    
    output$reference.values <- lapply(AC, function(i) {
      ctsm.lmm.refvalue(
        output, 
        yearID = max(data$year), 
        refvalue = qlogis(i / 100),
        lower.tail = switch(good_status, low = TRUE, high = FALSE)
      )
    })
    
    output$reference.values <- do.call("rbind", output$reference.values)
  }
  
  
  # construct summary output -
  
  output$summary <- data.frame(
    nyall = nYearFull, nyfit = nYear, nypos = nYearPos, 
    firstYearAll = firstYearFull, firstYearFit = min(data$year), lastyear = max(data$year), 
    p_nonlinear = NA, p_linear = NA, p_overall = NA, pltrend = NA, ltrend = NA, prtrend = NA, 
    rtrend = NA, dtrend = NA, meanLY = NA, clLY = NA)
  
  
  output$summary <- within(output$summary, {
    
    if (output$method == "smooth") {
      
      p_nonlinear <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["linear", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })
      
    }
    
    if (output$method == "linear") {
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- p_linear
      
    }
    
    if (output$method %in% c("linear", "smooth")) {
      
      # for linear trend and recent trend, use pltrend (from likelihood ratio test) if 
      # method = "linear", because a better test 
      # really need to go into profile likelihood territory here!
      
      pltrend <- if (output$method == "linear") {
        p_linear
      } else {
        output$contrasts["whole", "p"]
      }
      
      ltrend <- with(output$contrasts["whole", ], estimate / (end - start))
      
      if ("recent" %in% row.names(output$contrasts)) {
        prtrend <- if (output$method == "linear") {
          p_linear
        } else {
          with(output$contrasts["recent", ], p)
        }
        rtrend <- with(output$contrasts["recent", ], estimate / (end - start))
      }
    }
    
    # if parametric model cannot be fitted, use maximum index in last two monitoring years 
    # (if low values are good) or minimum index (if low values are bad) for crude extra data
    
    if (output$method == "none") 
      meanLY <- local({
        index <- tail(annualIndex$index, nYear)
        switch(
          good_status, 
          low = max(index), 
          high = min(index)
        )
      })
    else {
      meanLY <- tail(output$pred$fit, 1)
      meanLY <- 100 * plogis(meanLY)
      clLY <- switch(
        good_status, 
        low = tail(output$pred$ci.upper, 1), 
        high = tail(output$pred$ci.lower, 1)
      )
      clLY <- 100 * plogis(clLY)
    }
  })  
  
  if (!is.null(AC)) {
    output$summary <- data.frame(output$summary, do.call(cbind, lapply(names(AC), function(i) {
      
      value <- AC[i]
      diff <- with(output, if (method == "none") summary$meanLY - value else summary$clLY - value)
      
      # estimate number of years until meanLY reaches target - based on rtrend
      # might be already there but cl is too high
      
      maxYear <- max(data$year)
      bigYear <- 3000
      
      tillTarget <- with(output$summary, {
        
        if (good_status == "low") {
          
          if (is.na(value) || (meanLY >= value & is.na(rtrend)))
            NA
          else if (meanLY < value) 
            maxYear
          else if (rtrend >= 0)
            bigYear
          else {
            wk <- (qlogis(value / 100) - qlogis(meanLY / 100)) / rtrend
            wk <- round(wk + maxYear)
            min(wk, bigYear)
          }
          
        } else {
          
          if (is.na(value) || (meanLY <= value & is.na(rtrend)))
            NA
          else if (meanLY > value) 
            maxYear
          else if (rtrend <= 0)
            bigYear
          else {
            wk <- (qlogis(value / 100) - qlogis(meanLY / 100)) / rtrend
            wk <- round(wk + maxYear)
            min(wk, bigYear)
          }
          
        }
      })
      
      out <- data.frame(value, diff, tillTarget, below.result = NA)
      names(out) <- paste0(i, c("", "diff", "achieved", "below"))
      out
    })))
  }
  
  output$summary <- within(output$summary, {
    
    # and round for ease of interpretation
    
    p_nonlinear <- round(p_nonlinear, 4)
    p_linear <- round(p_linear, 4)
    p_overall <- round(p_overall, 4)
    pltrend <- round(pltrend, 4)
    prtrend <- round(prtrend, 4)
    
    ltrend <- round(ltrend, 3)
    rtrend <- round(rtrend, 3)
    dtrend <- round(dtrend, 3)
  })
  
  
  rownames(output$summary) <- NULL
  output
}



assess_negativebinomial <- function(
  data, annualIndex, AC, recent.years, determinand, max.year, recent.trend, 
  nYearFull, firstYearFull) {
  
  # silence non-standard evaluation warnings
  info <- weight <- NULL

  # over-dispersed count data (perhaps very low over-dispersed values from a 
  # binomial distribution, such an MNC) 
  
  # check valid determinands 
  
  if (! determinand %in% "MNC") {
    stop("not yet coded for determinand: ", determinand)
  }
  
  
  # initialise output
  
  output <- list(data = data)
  
  
  # set up offset - e.g. for MNC these are the number of individuals
  # specified in MNc-QC-NR
  
  if (!("offset" %in% names(data))) {
    data$offset <- 1
  }
  
  
  # get number of years and adjust year variable for stability
  # nYearPos is calculated for consistency with other responses
  
  nYear <- nYearPos <- length(unique(data$year))
  
  data$year_adj <- data$year - min(recent.years)
  
  data$year_fac <- factor(data$year)
  
  # establish other info
  
  good_status <- ctsm_get_info(info$determinand, determinand, "good_status")
  
  
  # type of fit depends on number of years:
  # nYear <= 2 none
  # nYear <= 4 mean 
  # nyear >= 5 linear or smooth
  # nYear >= 7, 10, 15 try smooths on 2, 3, 4 df
  
  # have only currently coded for mean and linear - look at ctsm.anyyear.lmm for 
  # extensions to smoothers
  
  if (nYear >= 3) {
    stop("time series too long: need to include code for smoothers")
  } 
  
  
  # do fits 
  
  if (nYear <= 2) 
    
    output$method <- "none"
  
  else {
    
    fits <- list(mean = NULL)
    
    # mean model
    
    fits$mean <- mgcv::gam(
      response ~ 1 + s(year_fac, bs = "re"), 
      weights = weight, 
      data = data, 
      family = "betar",
      method = "ML"
    )
    
    
    # linear model
    
    if (nYear >= 5) {
      fits$linear <- update(fits$mean, .~. + year_adj)
    }  
    
    
    # get basic anova 
    
    output$anova <- data.frame(
      p = sapply(fits, function(i) {
        edf <- i$edf
        id <- !grepl("year_fac", names(edf))
        edf <- edf[id]
        sum(edf) + 2
      }),
      twiceLogLik = sapply(fits, function(i) - 2 * i$gcv.ubre) 
    )
    
    output$anova$AIC <- - output$anova$twiceLogLik + 2 * output$anova$p
    
    
    # choose best model based on AIC
    
    bestFit <- if (nYear <= 4) 1 else max(2, which.min(output$anova$AIC))
    
    fit <- fits[[bestFit]]
    
    output$method <- 
      if (bestFit == 1) {
        "mean" 
      } else if (bestFit == 2) { 
        "linear" 
      } else {
        "smooth"
      }
    
    if (output$method == "smooth") {
      output$dfSmooth <- bestFit - 1
    }
    
    output$dfResid <- nYear - bestFit
    
    
    # predicted values with pointwise two-sided 90% confidence limits
    
    new_data <- data.frame(year = seq(min(data$year), max(data$year)))
    
    new_data$year_adj <- new_data$year - min(recent.years)
    
    new_data$year_fac <- factor(new_data$year)
    
    pred <- predict(fit, new_data, type = "lpmatrix")
    
    id <- !grepl("year_fac", dimnames(pred)[[2]])
    Xpred <- pred[, id, drop = FALSE]
    
    vcov <- vcov(fit)[id, id, drop = FALSE]
    coefficients <- coefficients(fit)[id]
    
    pred <- data.frame(
      year = new_data$year, 
      fit = c(Xpred %*% coefficients)
    )
    
    cov_fit <- Xpred %*% vcov %*% t(Xpred)
    pred$se <- sqrt(diag(cov_fit))
    
    output$conf_level_predictions <- 0.90
    
    alpha <- (1 - output$conf_level_predictions) / 2
    
    pred$ci.lower <- pred$fit + pred$se * qt(alpha, output$dfResid)
    pred$ci.upper <- pred$fit + pred$se * qt(1 - alpha, output$dfResid)
    
    output$pred <- pred
    
    row.names(output$pred) <- NULL
    
    
    # model coefficients
    coefficients <- data.frame(
      est = coefficients, 
      se = sqrt(diag(vcov))
    )
    
    coefficients$t <- coefficients$est / coefficients$se
    
    coefficients$p <- 2 * pt(
      abs(coefficients$t), 
      df = output$dfResid, 
      lower.tail = FALSE
    )
    
    coefficients$p = round(coefficients$p, 4)
    
    output$coefficients <- coefficients
  }
  
  
  # get estimated change in logit value over whole time series and in the 
  # most recent # e.g. twenty years of monitoring (truncate when data missing 
  # and only compute if at least five years in that period)
  # NB p value from contrast is NOT the same as from likelihood ratio test even 
  # if method = "linear"
  
  if (output$method %in% c("linear", "smooth")) {
    
    fit_info <- list(
      pred = output$pred, 
      dfResid = output$dfResid, 
      Xpred = Xpred,
      vcov = vcov
    )
    
    contrast.whole <- ctsm.lmm.contrast(fit_info, start = min(data$year), end = max(data$year))
    row.names(contrast.whole) <- "whole"
    
    start.year <- max(max.year - recent.trend + 1, min(data$year))
    if (sum(unique(data$year) >= start.year - 0.5) >= 5) {
      contrast.recent <- ctsm.lmm.contrast(fit_info, start = start.year, end = max(data$year))
      row.names(contrast.recent) <- "recent"
      contrast.whole <- rbind(contrast.whole, contrast.recent)
    }		
    
    output$contrasts <- contrast.whole
  }
  
  
  # compare mean value in final year to assessment criteria
  # results are presented on the logistic scale
  
  if (output$method %in% c("mean", "linear", "smooth")) {
    
    output$reference.values <- lapply(AC, function(i) {
      ctsm.lmm.refvalue(
        output, 
        yearID = max(data$year), 
        refvalue = qlogis(i / 100),
        lower.tail = switch(good_status, low = TRUE, high = FALSE)
      )
    })
    
    output$reference.values <- do.call("rbind", output$reference.values)
  }
  
  
  # construct summary output -
  
  output$summary <- data.frame(
    nyall = nYearFull, nyfit = nYear, nypos = nYearPos, 
    firstYearAll = firstYearFull, firstYearFit = min(data$year), lastyear = max(data$year), 
    p_nonlinear = NA, p_linear = NA, p_overall = NA, pltrend = NA, ltrend = NA, prtrend = NA, 
    rtrend = NA, dtrend = NA, meanLY = NA, clLY = NA)
  
  
  output$summary <- within(output$summary, {
    
    if (output$method == "smooth") {
      
      p_nonlinear <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["linear", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- with(output, {
        smoothID <- paste0("smooth (df = ", dfSmooth, ")")
        diff <- anova[smoothID, "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, dfSmooth - 1, lower.tail = FALSE)
      })
      
    }
    
    if (output$method == "linear") {
      
      p_linear <- with(output, {
        diff <- anova["linear", "twiceLogLik"] - anova["mean", "twiceLogLik"]
        pchisq(diff, 1, lower.tail = FALSE)
      })
      
      p_overall <- p_linear
      
    }
    
    if (output$method %in% c("linear", "smooth")) {
      
      # for linear trend and recent trend, use pltrend (from likelihood ratio test) if 
      # method = "linear", because a better test 
      # really need to go into profile likelihood territory here!
      
      pltrend <- if (output$method == "linear") {
        p_linear
      } else {
        output$contrasts["whole", "p"]
      }
      
      ltrend <- with(output$contrasts["whole", ], estimate / (end - start))
      
      if ("recent" %in% row.names(output$contrasts)) {
        prtrend <- if (output$method == "linear") {
          p_linear
        } else {
          with(output$contrasts["recent", ], p)
        }
        rtrend <- with(output$contrasts["recent", ], estimate / (end - start))
      }
    }
    
    # if parametric model cannot be fitted, use maximum index in last two monitoring years 
    # (if low values are good) or minimum index (if low values are bad) for crude extra data
    
    if (output$method == "none") 
      meanLY <- local({
        index <- tail(annualIndex$index, nYear)
        switch(
          good_status, 
          low = max(index), 
          high = min(index)
        )
      })
    else {
      meanLY <- tail(output$pred$fit, 1)
      meanLY <- 100 * plogis(meanLY)
      clLY <- switch(
        good_status, 
        low = tail(output$pred$ci.upper, 1), 
        high = tail(output$pred$ci.lower, 1)
      )
      clLY <- 100 * plogis(clLY)
    }
  })  
  
  if (!is.null(AC)) {
    output$summary <- data.frame(output$summary, do.call(cbind, lapply(names(AC), function(i) {
      
      value <- AC[i]
      diff <- with(output, if (method == "none") summary$meanLY - value else summary$clLY - value)
      
      # estimate number of years until meanLY reaches target - based on rtrend
      # might be already there but cl is too high
      
      maxYear <- max(data$year)
      bigYear <- 3000
      
      tillTarget <- with(output$summary, {
        
        if (good_status == "low") {
          
          if (is.na(value) || (meanLY >= value & is.na(rtrend)))
            NA
          else if (meanLY < value) 
            maxYear
          else if (rtrend >= 0)
            bigYear
          else {
            wk <- (qlogis(value / 100) - qlogis(meanLY / 100)) / rtrend
            wk <- round(wk + maxYear)
            min(wk, bigYear)
          }
          
        } else {
          
          if (is.na(value) || (meanLY <= value & is.na(rtrend)))
            NA
          else if (meanLY > value) 
            maxYear
          else if (rtrend <= 0)
            bigYear
          else {
            wk <- (qlogis(value / 100) - qlogis(meanLY / 100)) / rtrend
            wk <- round(wk + maxYear)
            min(wk, bigYear)
          }
          
        }
      })
      
      out <- data.frame(value, diff, tillTarget, below.result = NA)
      names(out) <- paste0(i, c("", "diff", "achieved", "below"))
      out
    })))
  }
  
  output$summary <- within(output$summary, {
    
    # and round for ease of interpretation
    
    p_nonlinear <- round(p_nonlinear, 4)
    p_linear <- round(p_linear, 4)
    p_overall <- round(p_overall, 4)
    pltrend <- round(pltrend, 4)
    prtrend <- round(prtrend, 4)
    
    ltrend <- round(ltrend, 1)
    rtrend <- round(rtrend, 1)
    dtrend <- round(dtrend, 1)
  })
  
  
  rownames(output$summary) <- NULL
  output
}


