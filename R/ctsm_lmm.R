# edit history

# 30/01/16 allow for SFG (normally distributed data)
# 31/01/16 deal with pathological data where all concentrations are equal
# 04/03/16 allow control parameters to be passed to optim
# 30/10/17 ctsm.lmm.fit replace bobyqa with L-BFGS-B in initial fit
# 08/12/17 ctsm.lmm.fit control fixed effects bounds
# 13/02/18 ctsm.lmm.fit control random effects bounds
# 14/02/18 ctsm.lmm.fit deal with pathological data where all within-year concentrations are equal 


ctsm.lmm.fit <- function(
  data, dfSmooth, varComp.start, fixed_bound = 5, random_bound = 3, controlLmmFit = NULL, ...) {

  library(lme4)
  library(mgcv)
  library(numDeriv)
  library(optimx)

  # get total number of years

  nYear <- length(unique(data$year))  
  

  # get design matrix, based on a single observation each year (to make the smoothers less 
  # sensitive to variation in sample size between years)
  # also ensure fitted values are flat in recent years that are all less-thans

  # create explanatory variable that is constant once there are no more 'real observations'

  maxYearPos <- with(data, max(year[censoring %in% ""]))
  data <- within(data, yearwk <- pmin(year, maxYearPos))
  

  # now get a central year to centre data
  
  dup <- duplicated(data$yearwk)
  yearAdj <- round(median(unique(data$yearwk)))
  data <- within(data, yearwk <- yearwk - yearAdj)

  
  # now get design matrix
  
  dummyFormula <- 
    if (dfSmooth == 0) as.formula(response ~ 1) else
      if (dfSmooth == 1) as.formula(response ~ yearwk) else
        as.formula(paste0("response ~ s(yearwk, k = ", dfSmooth + 1, ")"))

  dummyGam <- gam(dummyFormula, data = data[!dup, ])
  
  X <- as.data.frame(predict(dummyGam, type = "lpmatrix", newdata = data))
  names(X) <- if (dfSmooth == 0) "intercept" else
    c("intercept", paste0("year", seq(1, dfSmooth)))
  
  
  # get good starting values by doing the fit in lmer (if there are multiple samples in some years)  
  #   or lm (otherwise) but ignoring sdAnalytic and censoring
  # when variation in data less than smallest sdAnalytic (e.g. when all less-thans, but only some 
  #   recorded as such) can get pathological behaviour - so add in some 'random' noise to induce a  
  #   small amount of analytical uncertainty - NB don't want to make it truly random for 
  #   reproducibility
  # can get pathological behaviour when there are multiple samples in some year but all the values are the same - 
  #   again add in some noise

  lmerData <- cbind(data, X)
  lmerData <- within(lmerData, {
    year <- factor(year)
    
    if (sd(response) < min(sdAnalytic))
      response <- response + rep(c(-1, 1) * min(sdAnalytic), length = length(response))
    
    if (nYear < nrow(data) && max(tapply(response, year, sd), na.rm = TRUE) < 0.0001)
      response <- response + 0.1 * rep(c(-1, 1) * min(sdAnalytic), length = length(response))
  })

  

  if (nYear < nrow(data)) {

    lmerFormula <- as.formula(
      paste("response ~ 0 +", paste(names(X), collapse = "+"), "+ (1 | year)"))
    
    lmerFit <- lmer(lmerFormula, data = lmerData, control = lmerControl(
      optimizer = "optimx", optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)
    ))

    parInfo <- list(
      fixef = as.data.frame(summary(lmerFit)$coefficients)[1:2],
      varComp = data.frame(start = as.data.frame(VarCorr(lmerFit))$sdcor))
    
    parInfo <- within(parInfo, {
      names(fixef) <- c("start", "se")
      row.names(varComp) <- c("sdYear", "sdSample")
    })
    
  } else {
    
    lmFormula <- as.formula(paste("response ~ 0 +", paste(names(X), collapse = "+")))
    
    lmFit <- lm(lmFormula, data = lmerData)
    
    parInfo <- list(
      fixef = as.data.frame(summary(lmFit)$coefficients)[1:2],
      varComp = data.frame(start = summary(lmFit)$sigma))
    
    parInfo <- within(parInfo, {
      names(fixef) <- c("start", "se")
      row.names(varComp) <- c("sdYear")
    })
    
  }
  
  
  # need to give generous upper bound on standard deviations, in case some of them are estimated to 
  #   be zero
  # allow intercept to go below start - 5 * se as far as minimum response - usually works, but there 
  #   can be strange behaviour when all responses are equal
  
  parInfo <- within(parInfo, {
    fixef <- within(fixef, {
      lower <- start - fixed_bound * se
      upper <- start + fixed_bound * se
      lower[1] <- pmin(lower[1], min(data$response))
    })
    varComp <- within(varComp, {
      lower <- 0
      upper <- random_bound * max(start)
    })
  })
  
  # but use initial values of varComp (e.g. from a previous fit) if provided
  
  if (!missing(varComp.start)) {
    parInfo$varComp[c("start", "lower", "upper")] <- 
      varComp.start[c("est", "bound.lower", "bound.upper")]
  }
  
  # and shrink upper bound if necessary (weird quirk in optim)
  
  if (!is.null(controlLmmFit) && exists("varCompUpper", controlLmmFit, inherits = FALSE))
    parInfo$varComp$upper <- controlLmmFit$varCompUpper
    

  start <- with(parInfo, c(fixef$start, varComp$start))
  lower <- with(parInfo, c(fixef$lower, varComp$lower))
  upper <- with(parInfo, c(fixef$upper, varComp$upper))
  scale <- with(parInfo, c(fixef$se, varComp$upper / 5))
 
  names(start) <- with(parInfo, c(row.names(fixef), row.names(varComp)))
  
  X <- as.matrix(X)
  
  
  # initialise output 
  
  out <- list(
    nYear = nYear, dfSmooth = dfSmooth, pFixed = dfSmooth + 1, pRandom = nrow(parInfo$varComp))
  

  # do optimisation, but calculate hessian in a second step conditional on estimated variance 
  #   components - appears more likely to give a sensible answer!
  # unction hessian from numDeriv is far superior to the hessian from optim
  
  idFixed <- with(out, 1:pFixed)
  idRandom <- with(out, pFixed + (1:pRandom))

  out[c("idFixed", "idRandom")] <- list(idFixed, idRandom)


  # allow for some external control

  control <- list(parscale = scale)
  
  if (!is.null(controlLmmFit) && exists("trace", controlLmmFit, inherits = FALSE)) 
      control <- c(control, list(trace = controlLmmFit$trace, REPORT = 1))

  
  
  fit <- optim(start, 
               fn = function(par, data, X, id) {
                 beta <- par[id$fixed]
                 varComp <- par[id$random]
                 mu <- c(X %*% beta)
                 negTwiceLogLik(data, mu, varComp)
               }, 
               method = "L-BFGS-B", lower = lower, upper = upper, control = control, 
               data = data, X = X, id = list(fixed = idFixed, random = idRandom))
  
  
  out$twiceLogLik <- - fit$value
  
  out$convergence <- list(optim = fit[c("convergence", "message")])
  
  
  # check for fixed effects on bounds or variance components on upper bounds 
  # (usually due to too many less-thans)
  
  coefficients <- data.frame(est = fit$par, scale = scale, bound.lower = lower, bound.upper = upper)
  coefficients <- within(coefficients, {
    onBound <- abs(est - bound.lower) < 0.0001 | abs(est - bound.upper) < 0.0001
  })
  
  out$convergence$bound <- list(
    fixed = !any(coefficients$onBound[idFixed]),
    random = with(coefficients[idRandom, ], all(abs(est - bound.upper) > 0.0001))
  )
  
  
  # AIC: total p = pFixed + pRandom
  # AICc approximated by assuming n is nYearPos (# years with +ve observations) 
  # and that sdYear is dominant variance component so
  # number of parameters used in correction term is pFixed + 1
  
  out$nYearPos <- with(data, sum(tapply(censoring, year, function(x) any(x == ""))))
  out$dfResid <- with(out, nYearPos - pFixed)

  out <- within(out, {
    AIC <- - twiceLogLik + 2 * (pFixed + pRandom)
    #    AICc <- twiceLogLik + 2 * p * n / (n - p - 1)
    AICc <- - twiceLogLik + 2 * (pFixed + pRandom) * nYearPos / (nYearPos - (pFixed + 1) - 1)
  })
  

  # get predicted values
  
  pred <- data.frame(year = with(data, seq(min(year), max(year))))
  pred$yearwk <- pmin(pred$year, maxYearPos) - yearAdj

  Xpred <- predict(dummyGam, type = "lpmatrix", newdata = pred)
  
  pred$fit <- c(Xpred %*% coefficients$est[idFixed])
  pred <- pred[c("year", "fit")]
  

  out$data <- data
  out$X <- X
  out$Xpred <- Xpred
  out$coefficients <- coefficients
  out$pred <- pred

  out
}



ctsm.lmm.hess <- function(ctsm.ob, hess.d = 0.001, hess.r = 6, ...) {

  # now get standard errors on fixed effects only (since harder to get full vcov matrix with 
  # variance estimates if data are not well behaved)

  estimates <- with(ctsm.ob, coefficients$est[idFixed])
  names(estimates) <- with(ctsm.ob, row.names(coefficients)[idFixed])
  
  varComp <- with(ctsm.ob, coefficients$est[idRandom])
  names(varComp) <- with(ctsm.ob, row.names(coefficients)[idRandom])

  # default values of d and r are 0.1 and 4 respectively (contratry to help file)
  # d controls initial value of delta, r controls accuracy
  # lack of convergence can often be solved by increasing r, but at expense of speed
  
  ctsm.ob$hessian <- hessian(
    func = function(par, varComp, data, X) {
      beta <- par
      mu <- c(X %*% beta)
      negTwiceLogLik(data, mu, varComp) 
    }, 
    x = estimates, method.args = list(d = hess.d, r = hess.r),
    data = ctsm.ob$data, X = ctsm.ob$X, varComp = varComp)
  

  # and now use to get appropriate standard errors on parameter estimates and confidence intervals on 
  # fitted values - use 90% pointwise (two-sided) confidence bands on fitted values, so that these also 
  # give an upper (one-sided) 95% confidence limit for comparing with reference values

  within(ctsm.ob, {
    
    vcov <- 2 * solve(hessian)
    
    coefficients <- within(coefficients, se <- t <- p <- NA)
    coefficients[idFixed, ] <- within(coefficients[idFixed, ], {
      se <- sqrt(diag(vcov))
      t <- est / se
      p <- round(2 * pt(abs(t), dfResid, lower.tail = FALSE), 4)
    })
      
    pred <- within(pred, {
      se <- sqrt(diag(Xpred %*% vcov %*% t(Xpred)))
      ci.lower <- fit - qt(0.95, dfResid) * se
      ci.upper <- fit + qt(0.95, dfResid) * se
    })
    pred <- pred[c("year", "fit", "se", "ci.lower", "ci.upper")]
    
  })  
}


# utility function for likelihood calculation

ctsm.lmm.dcalc <- function(x, censoring, mean, sd, log = FALSE) {
  out <- ifelse(censoring == "", dnorm(x, mean, sd), pnorm(x, mean, sd)) 
  if (log) sum(log(out)) else prod(out)
}


#' @export
negTwiceLogLik <- function(data, mu, varComp) {

  library(mvtnorm)

  sdSample <- if ("sdSample" %in% names(varComp)) varComp["sdSample"] else 0
  sdYear <- varComp["sdYear"]
  
  data <- within(data, {
#    response <- log(concentration)
    mu <- mu
#    sdAnalytic <- aweight / concentration
    sdIndep <- sqrt(sdAnalytic ** 2 + sdSample ** 2)
    sdTotal <- sqrt(sdIndep ** 2 + sdYear ** 2)
  })


  logLik <- by(data, as.factor(data$year), function(data) {
    
    n <- nrow(data)

    if (n == 1 | sdYear <= 0.0001) {

      # univariate, or effectively independent observations, so can use standard functions
      
      with(data, ctsm.lmm.dcalc(response, censoring, mu, sdTotal, log = TRUE))

    } else if (all(data$censoring %in% "")) {

      # can use dmvnorm with impunity
        
      sigma <- (sdYear ** 2) * matrix(1, n, n) + (data$sdIndep ** 2) * diag(n)
      dmvnorm(data$response, data$mu, sigma = sigma, log = TRUE)

    } else {
        
      # correlated data, with a mixture of reals and less-thans (possibly all less-thans), so integrate
      # out the year effect, since conditional on this, everything is independent
      # need to ensure a scalar mean for z (the year effect) - mu is constant within years
      # NB don't use pmvnorm for optimization, without asking for great accuracy as it is an MC procedure
      
      zmu <- unique(data$mu)
      if (length(zmu) != 1) stop("different mean response within years")
       
      out <- try(integrate(function(z) {
        dYear <- dnorm(z, zmu, sdYear)
        dIndep <- sapply(z, function(x) ctsm.lmm.dcalc(data$response, data$censoring, x, data$sdIndep))
        dYear * dIndep
      }, zmu - 10 * sdYear, zmu + 10 * sdYear)$value)
      if (class(out) %in% "try-error")
        class(out)
      else 
        log(out)
    }      
  })

  if ("try-error" %in% logLik) return(100000)
  
  out <- - 2 * sum(logLik)
  if (out == "Inf") return(100000)

  out
}

