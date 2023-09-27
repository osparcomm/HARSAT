# 2_66 (OSPAR 2022)

# ctsm.VDS.cl - tidy up code dealing with Tritia (more needs to be done) and 
#  change factors to characters in output


#' ctsm.VDS.varlist
#' 
#' A list of names for functions and values necessary to export to
#' cluster prcesses for parallel computation
#' 
#' @export
ctsm.VDS.varlist <- paste("ctsm.VDS", c("p.calc", "loglik.calc", "index.opt", "cl"), sep = ".")

#' Detects the environment from the call
#' 
#' This is a utility function that detects the package environment.
#' If you have imported `harsat` as a package, it returns the package
#' environment. Otherwise, it returns the global environment. You can safely
#' export functions from the result of this, for example, when you are
#' setting up a cluster of child processes for parallel computation.
#' 
#' @export 
cstm.VDS.environment <- function() {
  environment(sys.function(sys.nframe()))
}

#' ctsm.VDS.p.calc
#' 
#' @param theta
#' @param cumulate a boolean, whether to use cumulative probabilities
#' @export
ctsm.VDS.p.calc <- function(theta, cumulate = FALSE) {
  if (cumulate) theta <- cumsum(theta)
  cumProb <- c(plogis(theta), 1)
  n <- length(cumProb)
  names(cumProb)[n] <- as.character(n-1)
  c(cumProb[1], diff(cumProb))
}

#' ctsm.VDS.loglik.calc
#' 
#' @param theta
#' @param data
#' @param index.theta
#' @param minus.twice
#' @param cumulate a boolean, whether to use cumulative probabilities
#' @export
ctsm.VDS.loglik.calc <- function(
  theta, data, index.theta, minus.twice = FALSE, cumulate = FALSE) {

  vds <- with(data, table(indexID, VDS))
  
  if (missing(index.theta)) {
    index.theta <- rep(0, nrow(vds))
    names(index.theta) <- row.names(vds)
  }
  
  out <- sapply(row.names(vds), function(x) {
    theta[1] <- theta[1] + index.theta[x]
    dmultinom(c(vds[x,]), prob = ctsm.VDS.p.calc(theta, cumulate), log = TRUE)
  })
  if (all(is.finite(out))) 
    out <- sum(out)
  else 
    out <- -1e7
  
  if (minus.twice) out <- - 2 * out
  
  out
}

#' ctsm.VDS.index.opt
#' 
#' @param data
#' @param theta
#' @param refLevel
#' @param calc.vcov
#' @export
ctsm.VDS.index.opt <- function(data, theta, refLevel, calc.vcov = FALSE) {

  data <- droplevels(data)
  
  if (missing(theta)) {
    theta <- with(data, table(VDS))
    theta <- cumsum(theta) / sum(theta)
    theta <- qlogis(theta[as.character(0:(max(data$VDS) - 1))])
    theta <- c(theta[1], diff(theta))
  }
  
  ntheta <- length(theta)
  nindex <- nlevels(data$indexID)
  
  index.theta <- rep(0, nindex)
  names(index.theta) <- levels(data$indexID)
  
  if (missing(refLevel)) {
    ref.choose <- with(data, tapply(VDS, indexID, mean))
    refLevel <- names(sort(ref.choose))[nindex %/% 2]
  }
  
  # add in a single category 1 observation for each index that is all zeros
  # and a single category n - 1 observation for each index that is all maxed out
  
  wk.zero <- with(data, tapply(VDS, indexID, function(x) all(x == 0)))
  wk.id <- with(data, indexID %in% names(wk.zero)[wk.zero] & !duplicated(indexID))
  data[wk.id, "VDS"] <- 1
  
  wk.max <- with(data, tapply(VDS, indexID, function(x) all(x == ntheta)))
  wk.id <- with(data, indexID %in% names(wk.max)[wk.max] & !duplicated(indexID))
  data[wk.id, "VDS"] <- ntheta - 1
  
  in.par <- c(theta, index.theta[-match(refLevel, names(index.theta))])
  
  wk.optim <- function(par, data, ntheta) {
    
    cutID <- as.character(0:(ntheta - 1))
    
    theta <- par[cutID]
    
    par <- par[setdiff(names(par), cutID)]
    
    index.theta <- rep(0, nlevels(data$indexID))
    names(index.theta) <- levels(data$indexID)
    index.theta[names(par)] <- par
    
    ctsm.VDS.loglik.calc(theta, data, index.theta, minus.twice = TRUE, cumulate = TRUE)
  }
  
  out <- optim(in.par, wk.optim, data = data, ntheta = ntheta, method = "L-BFGS-B", 
               lower = c(-Inf, rep(0, ntheta - 1), rep(-Inf, length(in.par) - ntheta)), 
               control = list(trace = 1, maxit = 500, REPORT = 10), hessian = calc.vcov)
  
  
  out$refLevel <- refLevel
  
  if (calc.vcov) {

    out$vcov <- 2 * solve(out$hessian)
    out$summary <- data.frame(est = out$par, se = sqrt(diag(out$vcov)))
    out$summary <- within(out$summary, {
      t <- est / se
      p <- round(2 * pnorm(abs(t), lower.tail = FALSE), 4)
    })
    out$summary <- subset(out$summary, select = c(est, se, t, p))
  }
  
  out$K <- ntheta
  
  out
}


#' ctsm.VDS.cl
#' 
#' @export
ctsm.VDS.cl <- function(fit, nsim = 1000) {
  
  nCuts <- fit$K
  cutsID <- as.character(0:(nCuts-1))
  categories <- 0:nCuts
  
  indexID <- setdiff(names(fit$par), cutsID)

  data <- MASS::mvrnorm(nsim, fit$par, fit$vcov)

  data.cuts <- data[, cutsID, drop = FALSE]
  if (nCuts > 1) data.cuts <- t(apply(data.cuts, 1, cumsum))

  data.index <- matrix(0, nrow = nsim, ncol = length(indexID) + 1, 
                       dimnames = list(NULL, c(indexID, fit$refLevel)))
  
  data.index[, indexID] <- data[, indexID]
  data.index <- as.data.frame(data.index)

  cl <- sapply(data.index, FUN = function(i) {
    out <- data.cuts + i
    out <- sort(apply(out, 1, function(x) sum(ctsm.VDS.p.calc(x) * categories)))
    out <- out[round(nsim * c(0.05, 0.95))]
  }, simplify = FALSE)

  cl <- data.frame(do.call("rbind", cl))
  names(cl) <- c("lower", "upper")

  n_tail = 2L
  if (any(grepl("Tritia nitida (reticulata)", row.names(cl), fixed = TRUE))) {
    # warning("ad-hoc fix for Tritia nitida (reticulata)")
    n_tail = 3L
  }
  
  namesID <- strsplit(row.names(cl), " ", fixed = TRUE)
  cl$species <- sapply(namesID, function(x) paste(tail(x, n_tail), collapse = " "))
  cl$year <- as.numeric(sapply(namesID, function(x) x[length(x) - n_tail]))
  cl$station_code <- sapply(
    namesID, 
    function(x) paste(head(x, length(x) - n_tail - 1), collapse = " ")
  )

  cl
}


# finds station year species combinations for which both individual and pooled 
# data have been submitted and checks for consistency

# ctsm.VDS.check <- function(ctsmOb) {
#   
#   data <- droplevels(subset(ctsmOb$data, determinand %in% determinands$Biota$imposex))
#   data <- data[c("seriesID", "station", "year", "species", "determinand", "concentration", "n_individual", 
#                  "%FEMALEPOP")]
#   
#   data[c("country", "region")] <- 
#     ctsmOb$stations[as.character(data$station), ][c("country", "region")]
#   data <- droplevels(data)
#   
#   data <- within(data, {
#     indexID <- factor(paste(station, determinand, species, year))
#     stopifnot(round(n_individual) == n_individual)
#     n_individual <- as.integer(n_individual)
#   })
#   
#   
#     # identify points with both individual and pooled data
#   
#   mixedData <- with(data, tapply(n_individual, indexID, function(x) any(x == 1) & any(x > 1)))
#   mixedData <- names(mixedData)[mixedData]
#   data <- droplevels(subset(data, indexID %in% mixedData))
#   
#   
#   # calculate VDSI based on individuals
#   
#   splitID <- factor(data$n_individual == 1, levels = c(FALSE, TRUE), labels = c("index", "stage"))
#   data <- split(data, splitID)
#   
#   stage <- with(data$stage, aggregate(concentration, by = list(indexID = indexID), function(x) 
#     c("SIndex" = mean(x), "SSum" = sum(x), "SN" = length(x))))
#   stage <- with(stage, cbind(indexID, as.data.frame(x)))
#   
#   data <- merge(data$index, stage, all = TRUE)
#   
#   names(data)[match("concentration", names(data))] <- "Index"
#   
#   data <- data[c("country", "station", "year", "species", "determinand", "n_individual", "%FEMALEPOP", 
#                  "Index", "SIndex", "SN", "SSum")]
#   names(data)[match("%FEMALEPOP", names(data))] <- "fprop"
#   
#   data <- within(data, {
#     okValue <- abs(Index - SIndex) <= 0.01  
#     okN <- round(n_individual * fprop / 100) == SN
#   })
#   
#   return(data)
# }
