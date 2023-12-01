# functions to assess imposex data

get.index.biota.Imposex <- function(data, determinand, info) {

  # silence non-standard evaluation warnings
  n_individual <- NULL

  data$nfemale <- round(data$n_individual * data[["%FEMALEPOP"]] / 100)
  if (any(is.na(data$nfemale)))
  {
    warning("missing FEMALEPOP - inferring 50%", immediate. = TRUE)
    data <- within(data, nfemale <- ceiling(n_individual / 2))
  }
  
  out <- by(data, data$year, function(x) with(x, 
    {
      index <- sum(concentration * nfemale) / sum(nfemale)
      nfemale <- sum(nfemale)
      data.frame(index = index, nfemale = nfemale)
    }))
  
  out <- do.call("rbind", out)

  out$year <- as.numeric(row.names(out))
  
  row.names(out) <- NULL
  
  out[c("year", "index", "nfemale")]
}        
  

imposex.varmean <- function(cuts) {
  
  delta <- seq(-30, 30, length = 1000)
  
  ilogit <- function(x) exp(x) / (1 + exp(x))
  
  g1 <- ilogit(cuts[1] - delta)
  g2 <- ilogit(cuts[2] - delta)
  g3 <- ilogit(cuts[3] - delta)
  g4 <- ilogit(cuts[4] - delta)
  g5 <- ilogit(cuts[5] - delta)
  
  g6 <- if (length(cuts) == 5) rep(1, length(delta)) else ilogit(cuts[6] - delta)
  
  p <- data.frame(p0 = g1, p1 = g2 - g1, p2 = g3 - g2, p3 = g4 - g3, p4 = g5 - g4, p5 = g6 - g5, p6 = 1 - g6)
  wk.mean <- p$p1 + 2 * p$p2 + 3 * p$p3 + 4 * p$p4 + 5 * p$p5 + 6 * p$p6
  wk.var <- p$p1 + 4 * p$p2 + 9 * p$p3 + 16 * p$p4 + 25 * p$p5 + 36 * p$p6 - wk.mean * wk.mean
  
  data.frame(mean = c(0, wk.mean, length(cuts)), var = c(0, wk.var, 0))
}


cuts6.all <- c(-6.69, -5.98, -4.39, -2.94, 4.80, 7.69)



imposex.link <- list(
  linkfun = function(mu) 
    log(mu/(6 - mu)),
  linkinv = function(eta) {
    junk <- exp(eta)
    6 * junk / (1 + junk)
  },
  mu.eta = function(eta) {
    junk <- exp(eta)
    6 * junk / ((1 + junk) ^ 2)
  },
  valideta = function(eta) TRUE,
  name = "6 point logit: log(mu / (6 - mu))"
)

cuts6.varmean <- imposex.varmean(cuts6.all)

imposex.variance <- list(
  name = "imputed from all non Quasimeme individual data", 
  variance = function(mu) {
    varmean <- cuts6.varmean
    varmean$mean <- varmean$mean 
    varmean$var <- varmean$var
    ifelse(mu %in% c(0, 6), 0, approx(varmean$mean, varmean$var, mu)$y)
  },
  deviance = function(mu, y, w, residuals = F) {
    devi.calc <- function(i, mu, y) {
      mu <- mu[i]
      y <- y[i]
      integrate(function(mu) (y - mu) / imposex.variance$variance(mu), mu, y)$integral
    }
    devi <- 2 * sapply(1:length(mu), devi.calc, mu = mu, y = y)
    if(residuals)
      sign(y - mu) * sqrt(abs(devi) * w)
    else sum(w * devi)
  }
)



imposex.family <- list(
  family = "imposex",
  link = "6 point logit: log(mu / (6 - mu))",
  linkfun = function(mu) 
    log(mu/(6 - mu)),
  linkinv = function(eta) {
    junk <- exp(eta)
    6 * junk / (1 + junk)
  },
  variance = function(mu) {
    # imputed from all non Quasimeme individual data
    varmean <- cuts6.varmean
    varmean$mean <- varmean$mean 
    varmean$var <- varmean$var
    ifelse(mu %in% c(0, 6), 0, approx(varmean$mean, varmean$var, mu)$y)
  },
  dev.resids = function(y, mu, w) {
    devi.calc <- function(i, mu, y) {
      mu <- mu[i]
      y <- y[i]
      out <- integrate(
        function(mu) {(y - mu) / imposex.family$variance(mu)}, 
        lower = mu, 
        upper = y
      )
      out$value
    }
    devi <- 2 * sapply(1:length(mu), devi.calc, mu = mu, y = y)
    #			sign(y - mu) * sqrt(abs(devi) * w)
    sqrt(abs(devi) * w)
  },
  aic = function(y, n, mu, weights, dev) NA,
  mu.eta = function(eta) {
    junk <- exp(eta)
    6 * junk / ((1 + junk) ^ 2)
  },
  initialize = expression({
    mustart<- rep(3, length(y))
  }), 
  validmu = function(mu) TRUE,
  valideta = function(eta) TRUE
)

#make.family("imposex", link = imposex.link, variance = imposex.variance)
#make.family




assess_imposex <- function(
    data, annualIndex, AC, recent.years, determinand, species, 
    station_code, thetaID, max.year, info.imposex, recent.trend = 20) {

  # silence non-standard evaluation warnings
  concentration <- NULL
  
  # main assessment routine for imposex data (imposex functions)
  
  # order data
  
  data <- data[order(data$year), ]
  
  nYearFull <- length(unique(data$year))  
  
  nYearFirst <- min(data$year)
  
  
  # deal with data sets that have crept in by mistake and have no recent data
  
  if (max(data$year) < min(recent.years)) return (NULL)
  
  
  # splits data into e.g. 6 year blocks and checks at least one observation in each block - ensures there 
  # aren't humungous gaps in the time series
  
  reporting.window <- length(recent.years)
  
  year.group <- cut(data$year, seq(max.year, min(data$year) - reporting.window, by = - reporting.window))
  
  in.group <- table(year.group) > 0
  
  if (!all(in.group)) {
    ok.group <- cumprod(rev(in.group)) == 1
    ok.group <- names(which(ok.group))
    
    id <- year.group %in% ok.group
    data <- data[id, ]
  }
  

  # all individual data, a mixture, or just indices
  
  indiID <- with(data, tapply(n_individual, year, function(x) all(x == 1L)))
  
  
  # if a mixture and the most recent three years of data are based on individuals
  # then just model the recent individual data
  
  if (any(indiID) & !all(indiID)) {
    indiRecent <- rev(cumprod(rev(indiID)) == 1)
    if (sum(indiRecent) >= 3) {
      okYears <- as.numeric(names(indiRecent)[indiRecent])
      data <- data[data$year %in% okYears, ]
    }
  }
 
  annualIndex <- annualIndex[annualIndex$year %in% unique(data$year), ]
    
  
  # initialise output

  output <- list(data = data)
  
  nYear <- length(unique(data$year))
  
  summary <- data.frame(
    nyall = nYearFull, nyfit = nYear, nypos = nYear, 
    firstYearAll = nYearFirst, firstYearFit = min(data$year), lastyear = max(data$year), 
    p_nonlinear = NA, p_linear = NA, p_overall = NA, 
    pltrend = NA, ltrend = NA, prtrend = NA, rtrend = NA, 
    meanLY = NA, clLY = NA, class = NA)
  

  # all individual data, a mixture, or just indices
  
  indiID <- with(data, tapply(n_individual, year, function(x) all(x == 1L)))
  
  # NB the following assumes that only a single imposex measure is assessed
  # for each station / species combination - need to make this bullet proof
  
  if (all(indiID) & thetaID %in% names(biota.VDS.estimates)) {
    
    # get theta estimates for modelling
    
    theta <- biota.VDS.estimates[[thetaID]]$par
    theta <- list(est = theta[names(theta) %in% as.character(0:5)])
    
    ntheta <- length(theta$est)
    theta$vcov <- biota.VDS.estimates[[thetaID]]$vcov[1:ntheta, 1:ntheta, drop = FALSE]
    
    # restrict categories where too few data to estimate cut points
    
    data <- within(data, {
      VDS <- concentration
      VDS[VDS > ntheta] <- ntheta
      VDS <- factor(VDS, levels = 0:ntheta)
    })
  
    assessment <- imposex_assess_clm(
      data, theta, annualIndex, species, recent.trend, max.year
    )
  }
  else {
  
    assessment <- imposex.assess.index(annualIndex, species, determinand, info.imposex)

    # ad-hoc adjustment if individuals reported in last monitoring year 
    # needed to deal with steep declines in imposex following the ban
    # take reported index and confidence limit if cl exists (not all species and 
    #   imposex measures) and cl is lower than that from reported trend 

    if ("upper" %in% names(annualIndex) && !is.na(tail(annualIndex$upper, 1))) {
      infoLY <- c(tail(annualIndex, 1))
      if (!("clLY" %in% names(assessment$summary)) || 
          infoLY$upper < assessment$summary$clLY) {
        assessment$summary <- within(assessment$summary, {
          meanLY <- infoLY$index
          clLY <- infoLY$upper
          class <- imposex_class(species, clLY)
        })
      }
    }
  }  
            
  
  # now add on the comparison to ACs - common to both methods  
  
  summary[names(assessment$summary)] <- assessment$summary
  assessment$summary <- NULL
    
  ACsummary <- lapply(names(AC), function(i) {
    upperLimit <- with(summary, if (!is.na(clLY)) clLY else meanLY) 
    diff <- upperLimit - AC[i]
    setNames(data.frame(AC[i], diff), paste0(i, c("", "diff")))
  })
  
  ACsummary <- do.call(cbind, ACsummary)
  
  summary <- cbind(summary, ACsummary) 

  output <- c(output, list(summary = summary), assessment)

  return(output)
}


imposex_class <- function(species, x) {
  
  # calculates imposex class (imposex functions)

  switch(
    species, 
    "Nucella lapillus" = dplyr::case_when(
      x < 0.3  ~ "A", 
      x < 2.0  ~ "B", 
      x < 4.0  ~ "C", 
      x <= 5.0 ~ "D",
      TRUE     ~ "E"
    ), 
    "Tritia nitida / reticulata" = dplyr::case_when(
      x < 0.3  ~ "B",
      x < 2.0  ~ "C",
      x <= 3.5 ~ "D",
      TRUE     ~ "F"
    ), 
    "Buccinum undatum" = dplyr::case_when(
      x < 0.3  ~ "B", 
      x < 2.0  ~ "C", 
      x <= 3.5 ~ "D",
      TRUE     ~ "F"
    ),
    "Neptunea antiqua" = dplyr::case_when(
      x < 0.3  ~ "A",
      x < 2.0  ~ "B",
      x <= 4   ~ "C",
      TRUE     ~ "F"
    ),
    "Littorina littorea" = dplyr::case_when(
      x < 0.3  ~ "C",
      x < 0.5  ~ "D",
      x <= 1.2 ~ "E", 
      TRUE     ~ "F"
    ),  
    NA
  )
}



imposex.assess.index <- function(annualIndex, species, determinand, info.imposex) {
  
  # silence non-standard evaluation warnings
  se <- NULL

  year <- annualIndex$year
  value <- annualIndex$index
  weights <- annualIndex$nfemale

  output <- list()
  summary <- list()
  
  nYear <- length(year)
  

  if (nYear <= 2) {
    summary$meanLY <- max(value)
    summary$class = imposex_class(species, max(value))
    return(list(summary = data.frame(summary)))
  }


  # catch series in which all values are equal - very ad-hoc

  if (diff(range(value)) == 0) {
    if (nYear > 3) {
      summary$p_linear <- summary$p_overall <- summary$pltrend <- summary$prtrend <- 1
      summary$ltrend <- summary$rtrend <- 0
    }
    summary$meanLY <- value[1]
    summary$clLY <- value[1]
    summary$class <- imposex_class(species, value[1])

    output$summary <- data.frame(summary)
        
    if (max(value) == 0) 
      output$pred <- data.frame(
        year = seq(min(year), max(year)),	
        fit = 0, 
        ci.lower = 0,
        ci.upper = 0.05
      )

    return(output)
  } 


  # catch series in which the first value is positive, but all others are zero
  # just add a small value to second index - very ad-hoc and need to resolve
  
  if (nYear > 3 & sum(value[-1]) == 0) 
    value[2:3] <- c(0.02, 0.01)
  

  # catch series in which the last value is positive, but all others are zero
  # just add a small value to second but last index - very ad-hoc and need to resolve
  
  if (nYear > 3 & sum(value[1:(nYear - 1)]) == 0) 
    value[c(nYear-2, nYear - 1)] <- c(0.01, 0.02)

  
  # now do assessment
  
  max.value <- info.imposex[
    info.imposex$species %in% species & info.imposex$determinand %in% determinand, 
    "max_value"
  ]
    
  # scale to unity for everything but Nucella, where this is incoporated in the family object
    
  if (species == "Nucella lapillus") 
    family.choice <- imposex.family
  else {
    value <- value / max.value		
    family.choice <- quasi(link = logit, variance = "mu(1-mu)")
  }
  
  if (nYear == 3) {
    fit <- glm(value ~ 1, weights = weights, family = family.choice, control = glm.control(maxit = 50))
    dfResid <- nYear - 1
  }  
  else {
    fit <- glm(value ~ year, weights = weights, family = family.choice, control = glm.control(maxit = 50))
    output$coefficients <- summary(fit)$coefficients
    dfResid <- nYear - 2
  }
  
  new.year <- seq(min(year), max(year))	
  pred <- predict(fit, newdata = data.frame(year = new.year), se.fit = TRUE)
  pred <- data.frame(fit = pred$fit, se = pred$se.fit)
  pred <- within(pred, {
    ci.lower <- fit + qt(0.05, dfResid) * se
    ci.upper <- fit + qt(0.95, dfResid) * se
  })
  pred <- data.frame(max.value * exp(pred) / (1 + exp(pred)))
  output$pred <- data.frame(year = new.year, pred[c("fit", "ci.lower", "ci.upper")])
  
  if (nYear > 3) {
    summary$p_linear <- summary$p_overall <- summary$pltrend <- summary$prtrend <- 
      round(output$coefficients["year", "Pr(>|t|)"], 4)
    summary$ltrend <- summary$rtrend <- round(output$coefficients["year", "Estimate"], 4)
  }
  summary$meanLY <- round(tail(pred$fit, 1), 3)
  summary$clLY <- round(tail(pred$ci.upper, 1), 3)
  summary$class <- imposex_class(species, tail(pred$ci.upper, 1))

  output$summary <- data.frame(summary)
  
  output
}