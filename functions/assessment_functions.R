# main assessment functions

# Set up functions ----

ctsm.assessment.setup <- function(ctsm.ob, AC = NULL, recent.trend = 20) {
  
  # set up assessment object (only needed for safety because trend assessment routines
  # are now much slower and don't want things to crash near the end, losing everything)

  ctsm.ob$call <- match.call()

  ctsm.ob <- within(ctsm.ob, {
    
    info$AC <- AC
    info$recent.trend <- recent.trend

    assessment <- vector(mode = "list", length = nrow(timeSeries))
    names(assessment) <- row.names(timeSeries)

  })
  
  ctsm.ob$call.data <- ctsm.ob$QA <- NULL
  ctsm.ob

}



ctsm.assessment <- function(
  ctsm.ob, determinandID, seriesID, get.assessment.criteria = get.AC, clusterID = NULL, ...) {

  # assess each time series
  
  require("pbapply")
  require("dplyr")
  
  info <- ctsm.ob$info

  timeSeries <- ctsm.ob$timeSeries
  
  if (! missing(determinandID)) 
    timeSeries <- droplevels(subset(timeSeries, determinand %in% determinandID))

  if (! missing(seriesID)) 
    timeSeries <- droplevels(timeSeries[seriesID, ])
  
  data <- droplevels(subset(ctsm.ob$data, seriesID %in% row.names(timeSeries)))

  stations <- ctsm.ob$stations
  stations <- stations[row.names(stations) %in% timeSeries$station, ]
  
  data <- split(data, data$seriesID)
  
  assessment <- pblapply(data, ..., cl = clusterID, FUN = function(x, ...) {
  
    # get info about the time series
    
    seriesID <- as.character(x$seriesID[1])
    seriesInfo <- sapply(timeSeries[seriesID,], function(i) if (is.numeric(i)) i else as.character(i), 
                         USE.NAMES = TRUE, simplify = FALSE)
    seriesInfo <- seriesInfo[vapply(seriesInfo, function(x) !is.na(x), NA)]
    
    print(paste("assessing series:", paste(names(seriesInfo), seriesInfo, collapse = "; ")))
    
    
    # return if no data to assess
    
    if (all(is.na(x$concentration))) return()
    
    if ("country" %in% names(stations))
      seriesInfo$country <- as.character(stations[seriesInfo$station, "country"])


    # construct annual index
    
    determinand <- seriesInfo$determinand
    
    var_id <- c("year", "concentration")
    
    var_id <- if (determinand %in% c("VDS", "IMPS", "INTS")) {
      append(var_id, c("noinp", "%FEMALEPOP"))
    } else {
      append(var_id, c("qflag", "uncertainty"))
    }
    
    if (determinand %in% "MNC") {
      var_id <- append(var_id, "MNC-QC-NR")
    }

    if (determinand %in% "%DNATAIL") {
      var_id <- append(var_id, "CMT-QC-NR")
    }
    
    x <- x[!is.na(x$concentration), var_id]
    
    row.names(x) <- NULL

    annualIndex <- get.index(info$compartment, determinand, x)
    

    # get assessment concentrations - need to extract some key variables first, could streamline in future
    
    if ("AC" %in% names(info))
      AC <- unlist(get.AC(info$compartment, determinand, seriesInfo, info$AC))
    else 
      AC <- NULL
    
    # initialise output from function
    
    out <- list(fullData = x, annualIndex = annualIndex, AC = AC)
    
    
    # do assessment
    
    if (determinand %in% c("VDS", "IMPS", "INTS")) {
      
      species <- seriesInfo$species
      station <- seriesInfo$station

      thetaID <- switch(
        info$purpose, 
        CSEMP = as.character(stations[station, "CMA"]),
        OSPAR = paste(seriesInfo$country, stations[station, "OSPAR_subregion"]),
        HELCOM = seriesInfo$country
      )

      thetaID <- paste(thetaID, species)
      
      

      # if any individual data, need to augment annual indices with confidence 
      # intervals
      
      indiID <- with(x, tapply(noinp, year, function(y) all(y == 1)))
      
      if (any(indiID) & thetaID %in% names(biota.VDS.estimates)) {

        out$annualIndex[c("lower", "upper")] <- NA
        
        clID <- paste(station, names(indiID)[indiID], species)
        out$annualIndex[indiID, c("lower", "upper")] <- 
          biota.VDS.cl[clID, c("lower", "upper")]
        
        # adjust when indices are zero or max (had to add an observation with a 
        # 1 or n-1 to get a fit)
        
        ntheta <- biota.VDS.estimates[[thetaID]]$K
        theta <- biota.VDS.estimates[[thetaID]]$par
        theta <- theta[as.character(0:(ntheta-1))]

        out$annualIndex <- within(out$annualIndex, {
          lower[index == 0] <- 0
          upper[index == ntheta] <- ntheta
        })  
      }  
      
      c(out, assess_imposex(
        data = x, 
        annualIndex = out$annualIndex, 
        AC = AC, 
        recent.years = info$recent_years,
        determinand = determinand, 
        species = species,
        station = station,
        thetaID = thetaID, 
        max.year = info$max_year, 
        recent.trend = info$recent.trend)
      )
    }
    else {

      # ad-hoc fix for missing uncertainties for SFG: 
      # trap for any lognormal or normal distributed data that have missing 
      #   uncertainties

      distribution <- get.info("determinand", determinand, "distribution")

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
        recent.trend = info$recent.trend)
      
      args.list <- c(args.list, list(...))
      fit <- try(do.call("ctsm.anyyear.lmm", args.list))
      if (is.character(fit) & length(fit) == 1) 
        return(c(out, error = fit))
      else
        return(c(out, fit))			
    }
    
  })
  
  assessment
}



# Annual indices ----

get.index <- function(compartment, determinand, data) {
  group <- get.info("determinand", determinand, "group", compartment)
  function_id <- paste("get.index", compartment, group, sep = ".")
  do.call(function_id, list(data = data, determinand = determinand)) 
}


get.index.default <- function(data, determinand) {

  data <- data[!is.na(data$concentration), ]
  
  # median (log) concentrations with flag to denote if less thans used in their 
  # construction
  
  distribution <- get.info("determinand", determinand, "distribution")
  
  data$response <- switch(
    distribution, 
    lognormal = log(data$concentration), 
    data$concentration
  )
  
  index <- tapply(data$response, data$year, median, na.rm = TRUE)
  
  qflag <- by(data, data$year, function(x) {
    x <- x[order(x$concentration),]
    n <- nrow(x)
    n0 <- ceiling(n / 2)
    any(x$qflag[n0:n] %in% c("<", "D", "Q"))
  })
  qflag <- sapply(qflag, function(i) i)	
  
  data <- data.frame(year = as.numeric(names(index)), index, qflag, row.names = NULL)
  
  data <- within(data, qflag <- factor(ifelse(qflag, "<", "")))
  
  data
}

get.index.sediment.Metals <- get.index.default
get.index.sediment.PAH_parent <- get.index.default
get.index.sediment.PAH_alkylated <- get.index.default
get.index.sediment.Chlorobiphenyls <- get.index.default
get.index.sediment.PBDEs <- get.index.default
get.index.sediment.Organobromines <- get.index.default
get.index.sediment.Organotins <- get.index.default
get.index.sediment.Organochlorines <- get.index.default
get.index.sediment.Organofluorines <- get.index.default
get.index.sediment.Dioxins <- get.index.default

get.index.biota.Metals <- get.index.default
get.index.biota.PAH_parent <- get.index.default
get.index.biota.PAH_alkylated <- get.index.default
get.index.biota.Chlorobiphenyls <- get.index.default
get.index.biota.PBDEs <- get.index.default
get.index.biota.Organobromines <- get.index.default
get.index.biota.Organotins <- get.index.default
get.index.biota.Organochlorines <- get.index.default
get.index.biota.Organofluorines <- get.index.default
get.index.biota.Dioxins <- get.index.default
get.index.biota.Metabolites <- get.index.default

get.index.biota.Effects <- function(data, determinand) {

  # median value apart from MNC and %DNATAIL where we get a weighted average
  # based on the number of cells

  # have put in a trap for qflags for anything that is not lognormally distributed
  # or that has good.status = high
  
  distribution <- get.info("determinand", determinand, "distribution")
  
  good_status <- get.info("determinand", determinand, "good.status")
  
  
  qflag_trap <- FALSE
  
  if (any(data$qflag != "") && good_status == "high") {
    qflag_trap <- TRUE
  } 
  
  if (any(data$qflag != "") && !distribution %in% c("normal", "lognormal")) {
    qflag_trap <- TRUE
  } 
  
  if (qflag_trap) {
    stop("surprising qflags: need to investigate")
  }
  
  
  # default for everything other than %DNATAIL and MNC: median (log) concentration
  
  if (!determinand %in% c("%DNATAIL", "MNC")) {
    out <- get.index.default(data, determinand)
    return(out)
  }
  
  
  # %DNATAIL and MNC: mean 'concentration' weighted by number of cells
  
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
  
  out$qflag <- factor(rep("", nrow(out)), levels = c("", "<"))
  
  out <- out[c("year", "index", "qflag", "nCell")]
    
  out$row.names <- NULL
  
  return(out)
}

get.index.water.Metals <- get.index.default
get.index.water.PAH_parent <- get.index.default
get.index.water.PAH_alkylated <- get.index.default
get.index.water.Chlorobiphenyls <- get.index.default
get.index.water.PBDEs <- get.index.default
get.index.water.Organobromines <- get.index.default
get.index.water.Organotins <- get.index.default
get.index.water.Organochlorines <- get.index.default
get.index.water.Organofluorines <- get.index.default
get.index.water.Dioxins <- get.index.default
get.index.water.Pesticides <- get.index.default


# Mixed model assessment functions ----

ctsm.anyyear.lmm <- function(data, annualIndex, AC, recent.years, determinand, max.year, 
                             recent.trend = 20, choose_model, ...) {

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
  
  
  # sort out assumed distribution (lognormal or normal) and whether high or low indicates good status
  
  distribution <- get.info("determinand", determinand, "distribution")
  good.status <- get.info("determinand", determinand, "good.status")
  

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
    wk_fn <- paste0("ctsm_assess_", distribution)
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
  
  nYearPos <- with(data, sum(tapply(qflag, year, function(x) any(x == ""))))
  

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
    
    start.year <- max(max.year - recent.trend + 1, min(data$year))
    if (sum(unique(data$year) >= start.year - 0.5) >= 5) {
      contrast.recent <- ctsm.lmm.contrast(fit, start = start.year, end = max(data$year))
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
        year = max(data$year), 
        refvalue = switch(distribution, lognormal = log(i), normal = i), 
        lower.tail = switch(good.status, low = TRUE, high = FALSE)
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
      
      pltrend <- if (output$method == "linear") p_linear else with(output$contrasts["whole", ], p)
      ltrend <- with(output$contrasts["whole", ], estimate / (end - start))
      
      if ("recent" %in% row.names(output$contrasts)) {
        prtrend <- if (output$method == "linear") p_linear else with(output$contrasts["recent", ], p)
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
      dtrend <- ctsm.dtrend(1:10, sigma, power = 0.9)
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

  check <- with(data, length(unique(year)) > 2 & any(qflag %in% c("<", "D", "Q")))
  if (!check) return(data)

  # identify years with at least one real
  
  wk <- with(data, tapply(qflag, year, function(x) any(x == "")))

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

  t.stat <- contrast / se.contrast
  p.contrast <- 1 - pf(t.stat^2, 1, ctsm.ob$dfResid)
  data.frame(start, end, estimate = contrast, se = se.contrast, p = p.contrast)
}


ctsm_check_convergence <- function(assessment_ob, coeff_se_tol = 0.001) {
  ok <- sapply(assessment_ob, function(x) {
    
    # trap errors in convergence 
    
    if ("error" %in% names(x)) 
      return(FALSE)

    # nothing to check if no model fitting (determined by presence of pred component)
    
    if (! "pred" %in% names(x))
      return(TRUE)
    
    # all standard errors from pred component should be present
      
    if (any(is.na(x$pred$se)))
      return(FALSE)
    
    # fixed effect coefficients: standard errors should be there and not ridiculously small
    # fixed effect coefficients should not be on bounds
    # random effect coefficients should not be on upper bounds
    
    coeff <- x$coefficients
    
    random_id <- grepl("sd", row.names(coeff))
    fixed <- coeff[!random_id, ]
    random <- coeff[random_id, ]

    if (any(is.na(fixed$se)) || min(fixed$se) < coeff_se_tol)
      return(FALSE)
    
    if (any(fixed$onBound))
      return(FALSE)
    
    if (any(random$onBound & random$est > 0.0001))
      return(FALSE)
    
    TRUE
  })
  
  names(ok[!ok])
}

ctsm.power <- function(
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


ctsm.dtrend <- function(
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
      ctsm.power(q, year, sigma, alpha, sigma_type, alternative) - power, 
    lower = lower, upper = upper, extendInt = extendInt
  )$root
}


ctsm.dyear <- function(
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
    achieved_power <- ctsm.power(q, 1:n_year, sigma, alpha, sigma_type, alternative)
  }
  
  n_year
}


# Other distributions ----

ctsm_assess_survival <- function(
  data, annualIndex, AC, recent.years, determinand, max.year, recent.trend, 
  nYearFull, firstYearFull) {

  # assess survival data, often expressed as interval data 
  
  library("flexsurv")


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
    time = if_else(.data$time == 0, .data$time2 / 10, .data$time)
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
  
  good_status <- get.info("determinand", determinand, "good.status")
  

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
    
    fits$mean <- flexsurvreg(
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
    
    contrast.whole <- ctsm_assess_survival_contrast(
      output, 
      start = min(data$year), 
      end = max(data$year)
    )
    row.names(contrast.whole) <- "whole"
    
    start.year <- max(max.year - recent.trend + 1, min(data$year))
    if (sum(unique(data$year) >= start.year - 0.5) >= 5) {
      contrast.recent <- ctsm_assess_survival_contrast(
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
      ctsm_assess_survival_refvalue(
        output, 
        year = max(data$year), 
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



ctsm_assess_survival_contrast <- function(ctsm.ob, start, end) {

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



ctsm_assess_survival_refvalue <- function(
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




ctsm_assess_beta <- function(
  data, annualIndex, AC, recent.years, determinand, max.year, recent.trend, 
  nYearFull, firstYearFull) {
  
  # percentage data that are not based on counts (e.g. comet assay) 
  
  library("mgcv")
  
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
  
  good_status <- get.info("determinand", determinand, "good.status")
  
  
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
    
    fits$mean <- gam(
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
        year = max(data$year), 
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



ctsm_assess_negativebinomial <- function(
  data, annualIndex, AC, recent.years, determinand, max.year, recent.trend, 
  nYearFull, firstYearFull) {
  
  # over-dispersed count data (perhaps very low over-dispersed values from a 
  # binomial distribution, such an MNC) 
  
  library("mgcv")
  
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
  
  good_status <- get.info("determinand", determinand, "good.status")
  
  
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
    
    fits$mean <- gam(
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
        year = max(data$year), 
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


# Post analysis functions ----

ctsm_post_analysis_power <- function(assessment_obj, target_power = 0.8) {
  
  lapply(assessment_obj$assessment, function(x) {
    
    # intialise output
    
    id <-c(paste("dtrend", 1:3, sep = "_"), "dyear", paste("dpower", 1:3, sep = "_"))
    
    x$summary[id] <- NA
    
    
    # get key data, and return if too few years to compute power
    
    year <- unique(x$data$year)
    n_year <- length(year) 
    
    sd <- x$sd_components["sd_index"]
    
    if (n_year < 3)
      return(x)
    
    
    # detectable trend (on log scale) of 
    # 1 current time series over observed time span
    # 2 current time series with no gaps (e.g. annual monitoring from min_year to max_year)
    # 3 in ten years of sequential monitoring
    
    if (n_year >= 5) {
      x$summary$dtrend_1 <- ctsm.dtrend(year, sd, power = target_power)
      x$summary$dtrend_2 <- ctsm.dtrend(min(year):max(year), sd, power = target_power)
    }
    
    x$summary$dtrend_3 <- ctsm.dtrend(1:10, sd, power = target_power)
    
    # back-transform to percentage annual (positive) change
    
    id <- paste("dtrend", 1:3, sep = "_")
    
    x$summary[id] <- lapply(x$summary[id], function(y) round(100 * (exp(y) - 1), 1))
    
    
    # number of sequential years to detect a 10% trend
    
    x$summary$dyear <- ctsm.dyear(log(1 + 0.1), sd, power = target_power)
    
    
    # power to detect an annual 10% change with same options as dtrend
    
    if (n_year >= 5) {
      x$summary$dpower_1 <- ctsm.power(log(1 + 0.1), year, sd)
      x$summary$dpower_2 <- ctsm.power(log(1 + 0.1), min(year):max(year), sd)
    }
    
    x$summary$dpower_3 <- ctsm.power(log(1 + 0.1), 1:10, sd)
    
    # turn into percentages
    
    id <- paste("dpower", 1:3, sep = "_")
    
    x$summary[id] <- lapply(x$summary[id], function(y) round(100 * y))
    
    x
  })  
}

