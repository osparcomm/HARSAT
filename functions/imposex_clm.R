# cumulative proportional odds model for imposex data

imposex.VDS.p.calc <- function(theta, cumulate = FALSE) {
  if (cumulate) theta <- cumsum(theta)
  cumProb <- c(plogis(theta), 1)
  n <- length(cumProb)
  names(cumProb)[n] <- as.character(n-1)
  c(cumProb[1], diff(cumProb))
}


imposex.clm.fit <- function(data, model, theta, model.control = list(), hessian = FALSE) {

  require(mgcv)
  
  out <- list(model = model, model.control = model.control)
  
  year <- unique(data$year)
  
  if (is.character(model))
    X <- imposex.clm.X(year, model, model.control)
  else
    X <- imposex.clm.X.change(year, model, model.control)
  
  ntheta <- length(theta$est)
  
  par <- rep(0, ncol(X))
  names(par) <- colnames(X)
  
  optimFUN <- function(par, data, theta, X) {
    eta <- setNames(as.vector(X %*% par), row.names(X))
    imposex.clm.loglik.calc(theta, data, eta, minus.twice = TRUE, cumulate = TRUE)
  }
  
  out <- c(out, optim(par, optimFUN, data = data, theta = theta$est, X = X, method = "L-BFGS-B", 
                      control = list(maxit = 500), hessian = hessian))
  
  out$X <- X
  out$nYear <- nrow(X)
  out$pFixed <- ncol(X)
  out$dfResid <- with(out, nYear - pFixed)
  out$twiceLogLik <- - out$value
  out$AIC <- with(out, value + 2 * pFixed)
  #    AICc <- - twiceLogLik + 2 * p * n / (n - p - 1)
  out$AICc <- with(out, if (pFixed < nYear - 1) value + 2 * pFixed * nYear / (nYear - pFixed - 1) else NA)
  
  out
}


imposex.clm.X <- function(
  year, model = c("mean", "linear", "smooth", "full"), model.control = list()) {

  data <- data.frame(year = year, row.names = year)

  # internal data are the years used to construct the smooth
  # the remaining data are those that are either interpolated or extrapolated
  # extrapolation is always linear

  control <- list(internal = year, adj = 2010, dfSmooth = NA)
  control[names(model.control)] <- model.control

  stopifnot(control$internal %in% year)
  if (model == "full" & !all(year %in% control$internal)) stop("can't construct X matrix for unsampled years")
  
  data <- within(data, {
    internal <- year %in% control$internal
    year <- year - control$adj
    response <- rnorm(length(year), year)
  })

  dummyFormula <- switch(
    model, 
    mean = as.formula(response ~ 1), 
    linear = as.formula(response ~ year), 
    smooth = as.formula(paste0("response ~ s(year, k = ", control$dfSmooth + 1, ")")),
    full = as.formula(response ~ as.factor(year))
  )
  
  dummyGam <- gam(dummyFormula, data = data[data$internal, ])
  
  X <- predict(dummyGam, newdata = data, type = "lpmatrix")
  
  np <- ncol(X)
  dimnames(X)[[2]] <- c("intercept", paste0("year", seq(1, np - 1)))[1:np]
  
  return(X)
}


imposex.clm.X.change <- function(year, model, model.control = list()) {

  stopifnot(c("before", "after", "change") %in% names(model)) 

  model <- within(model, {
    before <- match.arg(before, c("mean", "linear"))
    after <- match.arg(after, c("linear", "smooth"))
  })
  
  # model$change is the year when the change happened
  # dfSmooth is the df of the smoother after the break
  
  control <- list(internal = year, adj = 2010, dfSmooth = NA)
  control[names(model.control)] <- model.control
  
  
  # if change year not an element of year, then need to add in to deal with model$before = "mean"; need to 
  # have design matrix in that year to extrapolate back; or indeed if there is model$after = "linear" and there
  # is only one after point
  
  addChangeYear <- ! model$change %in% year
  if (addChangeYear) year <- sort(c(year, model$change))
  
  # ok to make internal == year when model$after == "linear" and deals with the special case 
  # when there is only a single year after the change point and no data in the change year
  
  if (model$after == "linear") 
    internal <- year  
  else
    internal <- year %in% control$internal

  dfBefore <- switch(model$before, mean = 0, linear = 1)
  dfAfter <- switch(model$after, mean = 1, linear = 2, smooth = control$dfSmooth + 1)
  
  X <- matrix(0, nrow = length(year), ncol = dfBefore + dfAfter)

  dimnames(X) <- list(
    as.character(year), 
    c("intercept", paste0("yearAfter", seq(1, dfAfter - 1)), if (dfBefore == 1) "yearBefore")
  )
  
  
  # fill in columns of X based on post-change years 
  # extrapolate back linearly, but then adjust below if necessary
  
  id <- year >= model$change
  idInternal <- id & internal

  X[, 1:dfAfter] <- imposex.clm.X(
    year, model$after, model.control = within(control, internal <- year[idInternal]))


  # if model$before = "mean", then replicate change point in design matrix
  # if model$before = "linear", then add in an additional linear term (relative to extrapolated linear term)
  
  if (model$before == "mean") 
    X[!id, ] <- X[rep(as.character(model$change), sum(!id)), ]

  if (model$before == "linear") 
    X[!id, dfAfter + 1] <- pmin(year[!id] - model$change, 0)
  
  # now remove any additional points needed in construction of design matrix
  
  if (addChangeYear)
    X <- X[-match(as.character(model$change), row.names(X)), ]
  
  X
}

  

imposex.clm.predict <- function(clmFit, theta, data) {

  year <- seq(min(data$year), max(data$year))
  
  args <- list(year = year, model = clmFit$model, model.control = clmFit$model.control)
  
  if (is.character(clmFit$model))
    Xpred <- do.call(imposex.clm.X, args)
  else
    Xpred <- do.call(imposex.clm.X.change, args)
  
  clmFit <- within(clmFit, {
    vcov.unscaled <- 2 * solve(hessian)
    vcov <- vcov.unscaled * disp
    
    summary <- data.frame(est = par, se = sqrt(diag(vcov)))
    summary <- within(summary, {
      t <- est / se
      p <- round(2 * pt(abs(t), dfResid, lower.tail = FALSE), 4)
    })
    summary <- summary[c("est", "se", "t", "p")]
  })
    
  clmFit$pred <- imposex.clm.cl(theta$est, Xpred, clmFit, theta$vcov)

  index <- with(data, tapply(as.numeric(as.character(VDS)), year, mean))
  clmFit$pred$index <- NA
  clmFit$pred[names(index), "index"] <- index
  
  # take the confidence bands that include the uncertainty in theta
  
  clmFit$pred <- clmFit$pred[c("year", "index", "fit", "lower2", "upper2")]
  names(clmFit$pred)[4:5] <- c("ci.lower", "ci.upper")
  
  clmFit$Xpred <- Xpred
  
  clmFit
}



imposex_assess_clm <- function(
    data, theta, annualIndex, species, recent.trend = 20, max.year) {

  output <- list()
  summary <- list()
  
  # decide whether there are sufficient years to model data
  # appropriate type of fit depends on total number of years and 
  # number of years with intermediate values (i.e between 0 and max(VDS))
  
  nYear <- length(unique(data$year))

  if (nYear <= 2) {
    summary$meanLY <- tail(annualIndex$index, 1)
    summary$clLY <- tail(annualIndex$upper, 1)
    summary$class = imposex_class(species, summary$clLY)
    return(list(summary = data.frame(summary)))
  }
  
  
  # do model fitting and comparison using correct data, as likelihood can still be calculated when there are
  # some years in which all observations are zeros or maximums
  # however, confidence intervals can often break down in this case (singularities, infinite estimates) so
  # adjust data to get 'sensible' confidence intervals
  # add in a single category 1 observation for each year that is all zeros
  # and a single category n - 1 observation for each year that is all maxed out

  ntheta <- length(theta$est)
  
  data <- within(data, {
    indexID <- factor(year)
    VDS2 <- VDS
  })
  
  all.min <- with(data, tapply(VDS, indexID, function(x) all(x %in% 0)))
  id <- with(data, indexID %in% names(all.min)[all.min] & !duplicated(indexID))
  data[id, "VDS2"] <- "1"
    
  all.max <- with(data, tapply(VDS, indexID, function(x) all(x %in% as.character(ntheta))))
  id <- with(data, indexID %in% names(all.max)[all.max] & !duplicated(indexID))
  data[id, "VDS2"] <- as.character(ntheta - 1)


  # get 'internal' years which are used to define a smoother - effectively this is to reduce the 
  # influence on the fit when there are many years with all maximum values at the start of a time series
  # or many years with all zero values at the end of a time series
  # if the data begins with a series of maximum values, then use the last one of these
  # if the data ends with a series of zeros, use the first one of these
  
  id1 <- pmax(1, sum(cumprod(all.max)))
  id2 <- nYear + 1 - pmax(1, sum(cumprod(rev(all.min))))

  year <- as.numeric(levels(data$indexID))
  internal <- year[id1:id2]
  nInternal <- length(internal)
  
  control <- list(internal = internal)


  # nYear <= 2 none
  # nYearPos >= 3 linear or smooth (pretty desperate with 3 or 4 years, 
  # but needed with year skipping in recent years)
  # nYearInt >= 5, 8, 11, try smooths on 2, 3, 4 df
  # fewer years required in change point model

  fits <- list(full = imposex.clm.fit(data, "full", theta))

  # now do the remaining fits 

  if (nYear >= 3) {
    fits[[2]] <- imposex.clm.fit(data, "mean", theta)
    fits[[3]] <- imposex.clm.fit(data, "linear", theta)
  }
  
  if (nInternal >= 5) 
    fits[[4]] <- imposex.clm.fit(data, "smooth", theta, model.control = c(control, dfSmooth = 2))
  
  if (nInternal >= 8) 
    fits[[5]] <- imposex.clm.fit(data, "smooth", theta, model.control = c(control, dfSmooth = 3))

  if (nInternal >= 11) 
    fits[[6]] <- imposex.clm.fit(data, "smooth", theta, model.control = c(control, dfSmooth = 4))

  anova <- do.call(rbind, lapply(fits, function(x) 
    as.data.frame(x[c("pFixed", "twiceLogLik", "AIC", "AICc")])))
  row.names(anova) <- 
    c("full", "mean", "linear", "smooth (df = 2)", "smooth (df = 3)", "smooth (df = 4)")[1:nrow(anova)]


  anova.change <- sapply(2004:2008, USE.NAMES = TRUE, simplify = FALSE, FUN = function(iy) {
  
    if (iy <= min(year) | iy >= max(year)) return(NULL)

    # must be at least one year after (not including) the change year

    nBefore <- sum(year < iy)
    nAfter <- sum(year > iy)
    nInternal <- sum(internal >= iy)
    
    beforeOptions <- if (nBefore > 10000) c("mean", "linear") else "mean"
    
    anova <- lapply(beforeOptions, function(iOpt) {
      
      model <- list(before = iOpt, after = "linear", change = iy)

      out <- list(linear = imposex.clm.fit(data, model, theta))

      if (nInternal >= 5) {
        model$after = "smooth"
        out$df2 <- imposex.clm.fit(data, model, theta, model.control = c(control, dfSmooth = 2))
      }
      
      if (nInternal >= 8) 
        out$df3 <- imposex.clm.fit(data, model, theta, model.control = c(control, dfSmooth = 3))
      
      if (nInternal >= 11) 
        out$df4 <- imposex.clm.fit(data, model, theta, model.control = c(control, dfSmooth = 4))
      
      out <- do.call(rbind, lapply(out, function(x) as.data.frame(x[c("pFixed", "twiceLogLik", "AIC", "AICc")])))

      id <- c("linear", "smooth (df = 2)", "smooth (df = 3)", "smooth (df = 4)")
      row.names(out) <- paste(iy, iOpt, id)[1:nrow(out)]
      
      out
    })
      
    do.call(rbind, anova)  
    
  })

  anova <- rbind(anova, do.call(rbind, anova.change))



  # estimate over dispersion by smallest estimate of sigma!

  dispOriginal <- with(anova[-1, ], {
    df <- anova["full", "pFixed"] - pFixed
    sigma <- (anova["full", "twiceLogLik"] - twiceLogLik) / df
    max(1, min(sigma))
  })
 
  # recalculate AIC and AICc (pFixed augmented by 1 for dispersion parameter)
   
  anova <- within(anova, {
    AICc.adj <- ifelse(
      pFixed < nYear - 2, 
      - twiceLogLik / dispOriginal + 2 * nYear * (pFixed + 1) / (nYear - (pFixed + 1) - 1), 
      NA)
    AIC.adj <- - twiceLogLik / dispOriginal + 2 * pFixed
  })
  
  
  # choose best model, but not including full or mean model 
  # round up AICc to 3rd dp to brush over rounding errors due to convergence and to 
  # encourage models that don't involve a step change, or an early step change
  # if few years, then AICc might not be defined - so run with AIC
  
  wk.anova <- anova[-c(1:2), ]
  wk.anova <- within(wk.anova, {
    AIC.adj <- ceiling(AIC.adj * 1000) / 1000
    AICc.adj <- ceiling(AICc.adj * 1000) / 1000
  })

  if (any(!is.na(wk.anova$AICc.adj)))
    bestFit <- row.names(wk.anova)[which.min(wk.anova$AICc.adj)]
  else
    bestFit <- row.names(wk.anova)[which.min(wk.anova$AIC.adj)]
  
  # refit best model - but don't get predictions because will probably fall over if lots of 
  # maxes or zeros
  
  best.id <- unlist(strsplit(bestFit, " "))
  
  fit <- switch(
    best.id[1], 
    linear = imposex.clm.fit(data, "linear", theta), 
    smooth = {
      dfSmooth <- as.numeric(substring(best.id[4], 1, 1))
      imposex.clm.fit(data, "smooth", theta, model.control = c(control, dfSmooth = dfSmooth))
    },
    {
      model <- list(before = best.id[2], after = best.id[3], change = as.numeric(best.id[1]))
      if (model$after == "smooth") 
        dfSmooth <- as.numeric(substring(best.id[6], 1, 1))
      else
        dfSmooth <- NA
      imposex.clm.fit(data, model, theta = theta, model.control = c(control, dfSmooth = dfSmooth))
    }
  )
  

  fit <- within(fit, {
    disp <- (anova["full", "twiceLogLik"] - twiceLogLik) / dfResid
    disp <- max(1, disp)
  })
  
  
  # get Hessian, but using data with no all zero or all maximum indices (since otherwise parameters are
  # unidentifiable / undefined) 
  # NB Not sure this is strictly true now that we have banned linear terms before the change point, but 
  # safer to leave it like this for now
  
  hessianInfo <- imposex.clm.fit(within(data, VDS <- VDS2), fit$model, theta, fit$model.control, 
                                 hessian = TRUE)

  fit[c("par2", "hessian")] <- hessianInfo[c("par", "hessian")]

  fit <- imposex.clm.predict(fit, theta, data)

  output$anova <- anova  
  output$method <- bestFit
  output <- c(output, fit[c("dfResid", "disp", "pFixed", "pred", "summary")])

  
  # get estimated change in concentration over whole time series and in the most recent 
  # e.g. twenty years of monitoring (truncate when data missing)
  # but if have fitted change point model, make recent trend from the change poing
  # NB p value from contrast is NOT the same as from likelihood ratio test even if method = "linear"
  # NB contrasts are actually the negative of the fitted values on the logistic scale, so that a 
  # negative value corresponds to a decrease in VDS

  contrast.whole <- imposex.clm.contrast(fit, start = min(data$year), end = max(data$year))
  row.names(contrast.whole) <- "whole"

  if (best.id[1] %in% c("linear", "smooth")) 
    start.year <- max(max.year - recent.trend + 1, min(data$year))
  else
    start.year <- as.numeric(best.id[1])
  contrast.recent <- imposex.clm.contrast(fit, start = start.year, end = max(data$year))
  row.names(contrast.recent) <- "recent"

  contrast.whole <- rbind(contrast.whole, contrast.recent)
  output$contrasts <- contrast.whole

  # compare fitted concentration in final year against assessment concentrations
  # ACs may be log-transformed to get onto the same scale as the index
  # some determinands have good status when high, some when low
  
  # NB Need to do quite a bit of work to implement this!
  
#   if (output$method %in% c("mean", "linear", "smooth")) 
#     output$reference.values <- do.call(rbind, lapply(AC, function(i) 
#       ctsm.lmm.refvalue(fit, year = max(data$year), 
#                         refvalue = switch(distribution, lognormal = log(i), normal = i), 
#                         lower.tail = switch(good.status, low = TRUE, high = FALSE))))
 
  
  
  # construct summary output - linID is the (piecewise) linear model corresponding to the 
  # best model (might be the best model itself, or the alternative)

  linID <- switch(best.id[1], linear = "linear", smooth = "linear", paste(best.id[1], "mean linear"))

  summary$p_linear <- summary$p_overall <- with(fit, {
    diff.lik <- anova[linID, "twiceLogLik"] - anova["mean", "twiceLogLik"]
    dfFixed <- 1
    Fstat <- (diff.lik / dfFixed) / disp
    pf(Fstat, dfFixed, dfResid, lower.tail = FALSE)
  })
  
  
  if (grepl("smooth", bestFit)) {
  
    summary$p_nonlinear <- with(fit, {
      diff.lik <- twiceLogLik - anova[linID, "twiceLogLik"]
      dfFixed <- pFixed - 2
      Fstat <- (diff.lik / dfFixed) / disp
      pf(Fstat, dfFixed, dfResid, lower.tail = FALSE)
    })
    
    summary$p_overall <- with(fit, {
      diff.lik <- twiceLogLik - anova["mean", "twiceLogLik"]
      dfFixed <- pFixed - 1
      Fstat <- (diff.lik / dfFixed) / disp
      pf(Fstat, dfFixed, dfResid, lower.tail = FALSE)
    })

  }
      


  summary$ltrend <- with(output$contrasts["whole", ], estimate / (end - start))
       
  # for trends, use pltrend (from likelihood ratio test) if method is linear or changepoint linear, 
  # because a better test - really need to go into profile likelihood territory here!
       
  if (grepl("linear", bestFit)) 
    summary$pltrend <- summary$prtrend <- summary$p_linear    
  else {
    summary$pltrend <- output$contrasts["whole", "p"]
    summary$prtrend <- output$contrasts["recent", "p"]
  }
  
  summary$rtrend <- with(output$contrasts["recent", ], estimate / (end - start))

  summary$meanLY <- tail(output$pred$fit, 1)
  summary$clLY <- tail(output$pred$ci.upper, 1)

  summary$class <- imposex_class(species, summary$clLY)
   
  summary <- within(summary, {
    # round for ease of interpretation
    
    if (grepl("smooth", bestFit))
      p_nonlinear <- round(p_nonlinear, 4)
    p_linear <- round(p_linear, 4)
    p_overall <- round(p_overall, 4)

    pltrend <- round(pltrend, 4)
    prtrend <- round(prtrend, 4)
    
    ltrend <- round(ltrend, 4)
    rtrend <- round(rtrend, 4)
  })
  

  output$summary <- data.frame(summary)  
  rownames(output$summary) <- NULL
  
  output
}

  

imposex.clm.loglik.calc <- function(theta, data, eta, minus.twice = FALSE, cumulate = FALSE) {
  
  # eta is mean VDS (relative to theta[1]) on logistic scale
  
  vds <- with(data, table(year, VDS))
  
  out <- sapply(row.names(vds), function(x) {
    theta[1] <- theta[1] + eta[x]
    dmultinom(c(vds[x,]), prob = imposex.VDS.p.calc(theta, cumulate), log = TRUE)
  })
  if (all(is.finite(out))) 
    out <- sum(out)
  else 
    out <- -1e7
  
  if (minus.twice) out <- - 2 * out
  
  out
}


imposex.clm.cl <- function(theta, X, fitOb, theta.vcov, nsim = 1000) {
  
  cutID <- names(theta)
  nCuts <- length(theta)
  categories <- c(0:6)[1:(nCuts+1)]
  
  yearID <- row.names(X)
  eta <- setNames(as.vector(X %*% fitOb$par), yearID)
  
  fit <- sapply(eta, FUN = function(x) {
    theta[1] <- theta[1] + x
    sum(imposex.VDS.p.calc(theta, cumulate = TRUE) * categories)
  })
  
  simEta <- MASS::mvrnorm(nsim, fitOb$par, fitOb$vcov)
  simEta <- t(X %*% t(simEta))
  simEta <- as.data.frame(simEta)
 
  data.cuts <- matrix(theta, nrow = nsim, ncol = nCuts, dimnames = list(NULL, cutID), byrow = TRUE)
  if (nCuts > 1) data.cuts <- t(apply(data.cuts, 1, cumsum))
  
  cl <- sapply(simEta, FUN = function(i) {
    out <- data.cuts + i
    out <- sort(apply(out, 1, function(x) sum(imposex.VDS.p.calc(x) * categories)))
    out <- out[round(nsim * c(0.05, 0.95))]
  }, simplify = FALSE)

  cl <- data.frame(do.call("rbind", cl))
  names(cl) <- c("lower", "upper")
  
  cl$year <- as.numeric(yearID)
  cl$fit <- fit
  

  # now add in uncertainty due to cuts - omit first element becuase it is effectively
  # aliased with the intercept (look at full vcov in biota.vds.estimates)

  theta.vcov[1, ] <- 0
  theta.vcov[, 1] <- 0
  
  data.cuts <- MASS::mvrnorm(nsim, theta, theta.vcov)
  if (nCuts > 1) data.cuts <- t(apply(data.cuts, 1, cumsum))
  
  cl2 <- sapply(simEta, FUN = function(i) {
    out <- data.cuts + i
    out <- sort(apply(out, 1, function(x) sum(imposex.VDS.p.calc(x) * categories)))
    out <- out[round(nsim * c(0.05, 0.95))]
  }, simplify = FALSE)
  
  cl2 <- data.frame(do.call("rbind", cl2))
  names(cl2) <- c("lower2", "upper2")
  
  cbind(cl[c("year", "fit", "lower", "upper")], cl2[c("lower2", "upper2")])
}



imposex.clm.contrast <- function(fitOb, start, end) {
  
  # almost the same function as for contaminants, but can't use predicted values of index
  # must use predicted values on logistic scale
  # could make them consistent by adjusting the lmm version
  # also have made contrasts negative, because +ve diff means decrease in VDS
  
  # error trapping

  if (length(start) > 1 | length(end) > 1) 
    stop('only a single contrast is allowed: start and end must both be scalars')  
  
  pos <- match(c(start, end), fitOb$pred$year)
  if (any(is.na(pos))) stop('start or end year not found in predicted data')
  
  wk <- t(c(-1, 1)) %*% fitOb$Xpred[pos, ]
  contrast <- - wk %*% fitOb$par

  se.contrast <- sqrt(wk %*% fitOb$vcov %*% t(wk))
  
  t.stat <- contrast / se.contrast
  p.contrast <- 1 - pf(t.stat^2, 1, fitOb$dfResid)
  data.frame(start, end, estimate = contrast, se = se.contrast, p = p.contrast)
}