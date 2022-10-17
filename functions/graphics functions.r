# edit history

# 17/12/14 label.units allow j/h/g
# 17/12/14 various ensure SFG not log transformed
# 05/01/15 various work out which variables are log transformed from info.determinand
#   simplify code for plotting AC labels
# 06/01/15 ctsm.webstats ensure revised summary column names are picked up correctly
#   do appropriate summaries for normally distributed variables (e.g. SFG)
#   ensure Pr(?t) has the correct < or > sign
# 07/01/15 ctsm.web.getKey allow GI and SB as multiple matrices
# 01/03/15 ctsm.webstats status assessment for imposex picks up from clLY rather than 
#   nyall (now that based on individual data)
# 19/10/15 ctsm.webstats use file.path to get appropriate path
# 28/01/16 change names of assessment objects (e.g. fullData, annualIndex)
# 06/06/16 sort out paths, so that consistent with file.path
# 13/06/16 ctsm.web.setup make filenames unique for e.g. biological effects with multiple timeseries
#   for the same determinand
# 16/06/16 plot.data pick up correct names for imposex assessments
# 15/07/16 ctsm.web.getKey fix typo in creation of start.text 
# 15/07/16 throughout info$recent.years -> info$recentYears
# 31/10/16 ctsm.webplot and ctsm.multiplot change path so that default correctly points to current directory
# 04/11/16 plot.auxiliary allow for auxiliary variables that don't have an associated qflag
# 05/11/16 label.units update since get.basis has changed to give NA for bioeffects and imposex
# 20/11/16 plot.auxiliary fix bug when an aux variable has no qflag and all values are missing
# 19/01/17 add in water
# 18/06/18 label.units remove 'wet weight' from units
# 22/10/19 pick up basis from timeSeries structure (not get.basis)
# 08/11/19 qflag can now be <, D and Q

# 2_61 OSPAR 2020
# ctsm.web.setup 
# - comment out filename construction (hopefully no longer needed)
# - convert some code to tidyverse
# plot.data.xlim 
# - ensure maximum year with tick mark is last possible data year (previously could exceed it
#   if a very long timeseries); not sorted for auxiliary data, which is a bit more fiddly; would be 
#   better to go direct to the axis function and control drawing of tick marks direct
# plot.panel
# - adjusted cex values so that multipanel assessments look better
# ctsm.web.getKey 
# - ad hoc fix to match matrix names to markdown names - needs to be revisited

# 2_64 OSPAR 2021
# plot.data 
# - plot closest AC outside range of data, even if there are some AC within the data
# panel.data
# - finetune cex for multipanel assessments
# - type argument extended to allow ratio plots and make multi_assessment plots explicit
# plot.multiassessment
# - finetune xykey.cex so year labels don't overlap (so much)
# - now only shows assessment indices rather than data
# plot.multidata
# - adjust height of varnames so that they fit
# ctsm.web.getKey
# - allow multiple units for sediment in multiplots (dioxins)
# plot.ratio, plot.ratio.data, plot.ratio.pred
# - introduce ratio plotting functions - still rather hardwired
# plot.data, plot.auxiliary, ctsm.web.getKey
# - replace 'concentration' where not appropriate - could be streamlined / made external

# 2_66 OSPAR 2022
# plot.ratio 
# - deal with situation when there are within-year replicates, but insufficient 
#   to fit a meaningful random effects model
# various - remove manual fix for survival distribution

# 2_67 CSSEG 2020
# ctsm.web.setup moved to reporting functions

ctsm.webplot <- function(assessmentObject, determinands, path = ".", filetype = c("jpg", "ps"), 
  xykey.cex = list(data = 1.4, auxiliary = 1.2), select = c("data", "assessment", "auxiliary"),  
  add.one = TRUE, ...) {

  require(lattice)
  require(grid)

  
  # check arguments
  
  filetype <- match.arg(filetype)
  
  select <- intersect(c("data", "assessment", "auxiliary"), select)
  if (length(select) == 0) stop('no valid graph types selected')

  
  # restrict to working determinands

  assessmentObject <- ctsm.subset.assessment(assessmentObject, determinand %in% determinands)

  
  # construct working data structures
  
  info <- assessmentObject$info
  data <- assessmentObject$data
  timeSeries <- assessmentObject$timeSeries
  
  
  outpath <- file.path(path, "graphics")

  timeSeries <- within(timeSeries, fileName <- paste(filePrefix, fileSuffix, sep = "_"))

  outnames <- expand.grid(select, timeSeries$fileName, stringsAsFactors = FALSE)
  outnames <- do.call("paste", rev(outnames))
  for (i in select)
    outnames <- gsub(paste(" ", i, sep = ""), substring(i, 1, 2), outnames)
  outnames <- paste(outnames, ".", filetype, sep = "")

  outwidth <- ceiling(log10(length(outnames)+1))            # need to add one to deal with dummy 
  # opening file (newPage should be FALSE, but messy to code)
  outfile <- file.path(outpath, paste0("tmpwebplot%0", outwidth, "d.", filetype))

  switch(filetype, 
    jpg = jpeg(outfile, quality = 100, width = 580, height = 480), 
    ps = postscript(outfile))


  data <- split(data, data$seriesID, drop = TRUE)
  

  
  # do the plots, and store whether graph exists

  outok <- lapply(rownames(timeSeries), function(id)
  {
    print(id)
    info <- c(info, timeSeries[id, ])

    data <- data[[id]]
    if (nrow(data) == 0 || all(is.na(data$concentration))) return(rep(FALSE, length(select)))

    assessment <- assessmentObject$assessment[[id]]
        
    info$matrix <- as.character(data$matrix[1])
    info$group <- as.character(data$group[1])

    for (i in 1:length(select))
    {
      switch(select[i], 
        data = plot.data(data, assessment, info, type = "data", xykey.cex = xykey.cex$data, ...),
        assessment = plot.data(data, assessment, info, type = "assessment", xykey.cex = xykey.cex$data, ...),
        auxiliary = plot.auxiliary(data, info, xykey.cex = xykey.cex$auxiliary, ...))
    }

    rep(TRUE, length(select))
  })

  dev.off()


  # rename files

  outnames <- outnames[unlist(outok)]

  if (add.one) 
  {
    outfiles <- file.path(
      outpath, 
      paste0("tmpwebplot", gsub(" ", "0", format(1:(length(outnames)+1), width = outwidth)), ".", filetype))
    file.remove(outfiles[1])    # dummy file caused by newPage = TRUE in first plot
  }
  else
    outfiles <- file.path(
      outpath, 
      paste0("tmpwebplot", gsub(" ", "0", format(1:(length(outnames)), width = outwidth)), ".", filetype))

  if (length(outnames) == 0) return(invisible())

  outnames <- gsub(" ", "_", outnames, fixed = TRUE)

  if (info$purpose == "CSEMP") 
  {
    if (max(nchar(outnames)) > 50) warning('some file names have more than 50 characters')
  }
  outnames <- file.path(outpath, outnames)
  
  if (add.one) outfiles <- outfiles[-1]
  lapply(1:length(outfiles), function(i) file.rename(outfiles[i], outnames[i]))

  invisible()
}






ctsm.webmultiplot <- function(assessmentWebObject, determinandGroups, path = ".",
  filetype = c("jpg", "ps"), select = c("data", "assessment"), ...) {

  require(lattice)
  require(grid)

  # check arguments for validity

  filetype <- match.arg(filetype)
  
  select <- intersect(c("data", "assessment"), select)
  if (length(select) == 0) stop('no valid graph types selected')
  

  # restrict to working determinandGroups - first need to get abbreviated names 
  # need to do this here to ensure compatability with parser file
  
  assessmentObject <- assessmentWebObject$assessment
  determinands <- assessmentWebObject$determinands
  

  # fileName identifies all timeSeries in same station / species / determinand group
  # convert detGroup to character just to be on safe side!
  
  assessmentObject$timeSeries <- within(assessmentObject$timeSeries, {
    detGroup <- as.character(detGroup)
    fileName <- factor(paste(filePrefix, fileGroup, sep = "_"))
  })
  
  assessmentObject <- ctsm.subset.assessment(assessmentObject, detGroup %in% determinandGroups)
  
  
  # construct working data structures
  
  info <- assessmentObject$info
  data <- assessmentObject$data
  timeSeries <- assessmentObject$timeSeries


  # sort, using determinands, to ensure everything is in the 'correct' presentational order
  
  timeSeries <- timeSeries[order(match(timeSeries$determinand, determinands)), ]
  
  
  # get determinands present in each group (could have done this neater and earlier perhaps)

  determinandGroups <- with(timeSeries, tapply(
    determinand, detGroup, function(x) unique(as.character(x)), simplify = FALSE))


  outpath <- file.path(path, "graphics")
  
  outnames <- expand.grid(select, levels(timeSeries$fileName), stringsAsFactors = FALSE)
  
  outnames <- do.call("paste", rev(outnames))
  for (i in select)
    outnames <- gsub(paste(" ", i, sep = ""), substring(i, 1, 2), outnames)
  outnames <- paste(outnames, filetype, sep = ".")

  # need to add one to deal with dummy opening file (newPage should be FALSE, but messy to code)                              
  outwidth <- ceiling(log10(length(outnames)+1))            
  outfile <- file.path(outpath, paste0("tmpwebplot%0", outwidth, "d.", filetype))   

  switch(filetype, 
    jpg = jpeg(outfile, quality = 100, width = 580, height = 480), 
    ps = postscript(outfile))

  
  timeSeries <- split(timeSeries, timeSeries$fileName, drop = FALSE)
  
  
  outok <- lapply(timeSeries, function(fSeries) {
  
    print(do.call("paste", fSeries[1, c("station", "detGroup")]))

    info$seriesID <- row.names(fSeries)
    
    data <- droplevels(data[data$seriesID %in% info$seriesID, ])
    assessment <- assessmentObject$assessment[info$seriesID]
    
    # to check - not sure which bits of fSeries are needed, presumably for the key, could simplify

    info <- c(info, fSeries[1, ])

    info$determinand <- as.character(fSeries$determinand)
    names(info$determinand) <- info$seriesID

    
    # get series names for labelling plots - usually just determinand, but could be more complex if e.g. 
    # measured in multiple tissues, or when dealing with biological effects
    
    info$plotNames <- list(data = info$determinand, assessment = info$determinand)
    if (any(duplicated(info$plotNames))) {
      dups <- duplicated(info$determinand) | duplicated(info$determinand, fromLast = TRUE)
      info$plotNames$data[dups] <- paste(info$determinand, fSeries$level6name, sep = "\n")[dups]
      info$plotNames$assessment[dups] <- paste(info$determinand, fSeries$level6name)[dups]
      if (any(duplicated(info$plotNames$data)))
        stop("duplicate plotting names - need to extend coding")
    }
    
    info$matrix <- with(data, tapply(matrix, seriesID, function(x) as.character(x[1])))
    info$group <- as.character(get.info("determinand", info$determinand, "group", info$compartment))
    
    if ("data" %in% select) 
      plot.multidata(data, info, ...)
    
    if ("assessment" %in% select) 
      plot.multiassessment(data, assessment, info, type = "assessment", ...)

    rep(TRUE, length(select))
  })

  dev.off()


  
  # rename files

  outnames <- outnames[unlist(outok)]

  outfiles <- file.path(
    outpath, 
    paste0("tmpwebplot", gsub(" ", "0", format(1:(length(outnames)), width = outwidth)), ".", filetype)
  )
#  outfiles <- paste(outpath, "tmpwebplot", gsub(" ", "0", format(1:(length(outnames)+1), width = outwidth)), 
#    ".", filetype, sep = "")
#  file.remove(outfiles[1])    # dummy file caused by newPage = TRUE in first plot  - doesn't seem to be 
#  needed here need to investigate why

  if (length(outnames) == 0) return(invisible())

  outnames <- gsub(" ", "_", outnames, fixed = TRUE)

  if (max(nchar(outnames)) > 50) warning('some file names have more than 50 characters')

  outnames <- file.path(outpath, outnames)
  
#  outfiles <- outfiles[-1]
  lapply(1:length(outfiles), function(i) file.rename(outfiles[i], outnames[i]))

  invisible()
}


ctsm.format <- function(x, y = x, nsig = 3) {

  # get arguments to format y to nsig significant figures
  # and use these to format x
  
  wk <- floor(log10(abs(y))) + 1 # digits to left of dp
  digits <- max(wk, nsig)
  nsmall <- max(0, nsig - wk)
  format(x, digits = digits, nsmall = nsmall, scientific = nsig)
}




ctsm.webstats <- function(assessmentObject, determinands, path = ".", printLevel = 10) {

  require(lattice)
  require(grid)
  require(hwriter)


  # restrict to working determinands
  
  assessmentObject <- ctsm.subset.assessment(assessmentObject, determinand %in% determinands)
  
  
  # construct working data structures
  
  info <- assessmentObject$info
  timeSeries <- assessmentObject$timeSeries

  timeSeries <- within(timeSeries, fileName <- paste0(filePrefix, "_", fileSuffix, ".html"))
  if (max(nchar(timeSeries$fileName)) > 50) warning('some file names have more than 50 characters')
  
  
  # get matrix and group for each time series for use in key function
  
  matrixID <- with(assessmentObject$data, 
    tapply(matrix, seriesID, function(x) as.character(x[1])))
  
  groupID <- with(assessmentObject$data, 
    tapply(group, seriesID, function(x) as.character(x[1])))
  
  
  stylename <- paste(path, 'ctsm_style.css', sep = "")

  outpath <- file.path(path, "graphics")
  
  hwrite.gap <- function(page, nLines = 2) hwrite(rep("", nLines), page, br = TRUE, table = FALSE)

  lapply(1:nrow(timeSeries), function(i) {

    id <- rownames(timeSeries)[i]
    if (i %/% printLevel == i / printLevel) print(id)
    
    assessment <- assessmentObject$assessment[[id]]
    
    if (is.null(assessment$summary)) return()

    info <- c(info, timeSeries[id, ])
  
    info$matrix <- matrixID[id]
    info$group <- groupID[id]


    # open html page

    outPage <- openPage(info$fileName, dirname = outpath, link.css = stylename)


    # key

    key <- ctsm.web.getKey(info, html = TRUE)
    hwrite(unlist(key), outPage, br = TRUE, table = FALSE)
      

    # useful formatting functions
    
    hformat <- function(x, k) 
      ifelse(is.na(x), "", format(round(x, k), nsmall = k, scientific = FALSE))

    hstyle <- function(x) 
      matrix(c('text-align:left', rep('text-align:right', ncol(x))), 
             nrow = nrow(x) + 1, ncol = ncol(x) + 1, byrow = TRUE)
    
    hpformat <- function(x, k = 4) 
      ifelse(x < 10^(-k), 
             paste0("<", hformat(10^(-k), k)), 
             hformat(x, k))
    
    hzap <- function(x, digits = 4)
      format(zapsmall(x, digits = digits), scientific = FALSE)
    
    
    # trend assessment: analysis of variance and change in last decade
                   
    hwrite.gap(outPage)
    hwrite("Trend assessment", outPage, br = TRUE, div = TRUE, class = "title")
    
    switch(
      as.character(info$detGroup), 
      Imposex = {
        
        if(is.null(assessment$coefficients) & is.null(assessment$anova))
          
          hwrite("Insufficient data", outPage, br = TRUE)
        
        else if (!is.null(assessment$coefficients)) {

          # estimates from index based assessment - hopefully to be phased out
          
          hwrite("Estimate of trend", outPage, br = TRUE, div = TRUE, class = "subtitle")
          
          coeff <- as.data.frame(assessment$coefficients)
          coeff <- coeff["year", ]
          row.names(coeff) <- "trend"
          
          names(coeff) <- c("est", "se", "t", "p")

          id <- c("est", "se")
          coeff[id] <- hzap(unlist(coeff[id]), digits = 4)
          
          id <- "t"
          coeff[id] <- lapply(coeff[id], hformat, k = 3)
          
          id <- "p"
          coeff[id] <- lapply(coeff[id], hpformat, k = 4)
          
          names(coeff) <- c("Estimate", "Std error", "t", "Pr(>|t|)")
          
          hwrite(coeff, outPage, border = 0, cellpadding = 2, style = hstyle(coeff))
        }
        else {
          
          hwrite("Analysis of variance", outPage, br = TRUE, div = TRUE, class = "subtitle")
          
          anova <- within(assessment$anova, {
            deviance <- - twiceLogLik
            logLik <- twiceLogLik / 2
            AIC <- AIC.adj
            AICc <- AICc.adj
            dfResid <- max(pFixed) - pFixed
            disp <- ifelse(dfResid > 0, (max(twiceLogLik) - twiceLogLik) / dfResid, NA)
            disp <- pmax(1, disp)
            rm(twiceLogLik, AIC.adj, AICc.adj)
          })
          
          if (nrow(anova) <= 2) stop("no models fitted - need to recode")
          
          
          # get relevant rows of anova table (don't want to show all change point model fits)
          # if change point model optimal, also include linear and smooth fit
          
          anovaID <- sapply(strsplit(row.names(anova), " "), "[[", 1)
          
          modelID <- unlist(strsplit(assessment$method, " "))[1]
          isChange <- !(modelID %in% c("linear", "smooth"))
          
          ok <- anovaID %in% c("mean", "linear", "smooth", modelID[1])
          anova <- anova[ok, ]
          anovaID <- anovaID[ok]
          
          if (isChange) {
            
            # order anova table so that change point models follow mean model, with 
            # linear and smooth later on - also relabel row names to make more interpretable
            
            row.names(anova) <- sapply(row.names(anova), function(x) {
              xsplit <- unlist(strsplit(x, " "))
              if (xsplit[1] != modelID) return(x)
              out <- paste(xsplit[-c(1:2)], collapse = " ")
              paste(out, "from", xsplit[1])
            })
            
            ok <- anovaID %in% c("mean", modelID)
            anova <- rbind(anova[ok, ], anova[!ok, ])
            ok <- c(ok[ok], ok[!ok])
            
            anova["other models:", ] <- NA
            anova <- anova[order(c(1:length(ok), sum(ok)+0.5)), ]  
            ok <- append(ok, FALSE, sum(ok))
          }
          
          
          # now in correct order get F values and p values
          
          anova <- within(anova, {
            Fstat <- c(NA, - diff(deviance)) / disp
            Fdf <- c(NA, diff(pFixed))
            if (isChange) {
              Fstat[!ok] <- NA
              Fdf[!ok] <- NA
              dfResid[!ok] <- NA
            }
            p <- pf(Fstat, Fdf, dfResid, lower.tail = FALSE)
          })
          
          anova["mean", "dfResid"] <- NA
          
          # format columns
          
          id <- c("pFixed", "dfResid", "Fdf")
          anova[id] <- lapply(anova[id], hformat, k = 0)
          
          id <- c("AIC", "AICc", "logLik", "deviance", "disp")
          anova[id] <- lapply(anova[id], hformat, k = 2)
          
          id <- "Fstat"
          anova[id] <- lapply(anova[id], hformat, k = 3)
          
          id <- "p"
          anova[id] <- lapply(anova[id], hpformat, k = 4)

          anova <- anova[c("pFixed", "AIC", "AICc", "logLik", "deviance", "disp", "Fstat", "Fdf", 
                           "dfResid", "p")]
          
          names(anova) <- c("Df", "AIC", "AICc", "Log lik", "Deviance", "Dispersion", 
                            "F", "F_df1", "F_df2", "Pr(>F)")
          
          hwrite(anova, outPage, border = 0, cellpadding = 2, style = hstyle(anova))
          
          
          # change in e.g. last decade
          
          hwrite.gap(outPage)
          
          hwrite("Change in logit cumulative probability", outPage, class = "subtitle", br = TRUE)
          hwrite("Fit on stage scale, change on logit scale", outPage, div = TRUE, class = "subtitle", br = TRUE)
          
          contrasts <- assessment$contrasts 
          
          row.names(contrasts)[1] <- "overall"
          if (nrow(contrasts) == 2)
            if (isChange)
              row.names(contrasts)[2] <- paste("since", modelID)
            else
              row.names(contrasts)[2] <- paste("last", assessmentObject$info$recent.trend, "years")
          
          pred <- assessment$pred
          row.names(pred) <- pred$year
          
          contrasts$fit1 <- pred[as.character(contrasts$start),"fit"]
          contrasts$fit2 <- pred[as.character(contrasts$end),"fit"]
          
          contrasts <- within(contrasts, t <- estimate / se)
          
          contrasts <- contrasts[c("start", "end", "fit1", "fit2", "estimate", "se", "t", "p")]
          
          id <- c("start", "end")
          contrasts[id] <- lapply(contrasts[id], hformat, k = 0)
          
          id <- c("fit1", "fit2")
          contrasts[id] <- lapply(contrasts[id], hformat, k = 3)
          
          id <- c("estimate", "se")
          contrasts[id] <- hzap(unlist(contrasts[id]), digits = 4)
          
          id <- "t"
          contrasts[id] <- lapply(contrasts[id], hformat, k = 3)
          
          id <- "p"
          contrasts[id] <- lapply(contrasts[id], hpformat, k = 4)
          
          names(contrasts) <- c("Year start", "Year end", "Fit start", "Fit end", "Change", 
                                "Std error", "t", "Pr(>|t|)")
          
          hwrite(contrasts, outPage, border = 0, cellpadding = 2, style = hstyle(contrasts))
        }
      },
      { 
        anova <- assessment$anova

        if(is.null(anova) || nrow(anova) == 1)
          hwrite("Insufficient data", outPage, br = TRUE)
        else {      

          # anova table
          
          hwrite("Analysis of variance", outPage, br = TRUE, div = TRUE, class = "subtitle")
         
          anova <- within(anova, {
            df <- 1:nrow(anova)
            deviance <- - twiceLogLik
            logLik <- twiceLogLik / 2
            Chisq <- c(NA, - diff(deviance))
            Chidf <- c(NA, diff(df))
            p <- pchisq(Chisq, Chidf, lower.tail = FALSE)
            rm(twiceLogLik)
          })

          id <- c("df", "Chidf")
          anova[id] <- lapply(anova[id], hformat, k = 0)
          
          id <- c("AIC", "AICc", "logLik", "deviance")
          anova[id] <- lapply(anova[id], hformat, k = 2)

          id <- "Chisq"
          anova[id] <- lapply(anova[id], hformat, k = 3)

          id <- "p"
          anova[id] <- lapply(anova[id], hpformat, k = 4)

          anova["mean", c("Chisq", "Chidf", "p")] <- ""
          
          anova <- anova[c("df", "AIC", "AICc", "logLik", "deviance", "Chisq", "Chidf", "p")]
          names(anova)[c(1, 4, 5, 7, 8)] <- c("Df", "Log lik", "Deviance", "Chi df", "Pr(>Chisq)")
          
          hwrite(anova, outPage, border = 0, cellpadding = 2, style = hstyle(anova))


          # change in e.g. last decade
          
          hwrite.gap(outPage)
          
          if (get.info("determinand", info$determinand, "distribution") %in% "lognormal") 
            txt <- "Change in log concentration"
          else 
            txt <- "Change in concentration"
      
          hwrite(txt, outPage, div = TRUE, class = "subtitle", br = TRUE)

          contrasts <- assessment$contrasts 
          row.names(contrasts) <- c(
            "overall", 
            paste("last", assessmentObject$info$recent.trend, "years")
          )[1:nrow(contrasts)]
  
          pred <- assessment$pred
          row.names(pred) <- pred$year
          
          contrasts$fit1 <- pred[as.character(contrasts$start),"fit"]
          contrasts$fit2 <- pred[as.character(contrasts$end),"fit"]
          
          contrasts <- within(contrasts, t <- estimate / se)

          contrasts <- contrasts[c("start", "end", "fit1", "fit2", "estimate", "se", "t", "p")]

          id <- c("start", "end")
          contrasts[id] <- lapply(contrasts[id], hformat, k = 0)
          
          id <- c("fit1", "fit2")
          contrasts[id] <- hzap(unlist(contrasts[id]), digits = 4)

          id <- c("estimate", "se")
          contrasts[id] <- hzap(unlist(contrasts[id]), digits = 4)

          id <- "t"
          contrasts[id] <- lapply(contrasts[id], hformat, k = 3)

          id <- "p"
          contrasts[id] <- lapply(contrasts[id], hpformat, k = 4)
          
          names(contrasts) <- c(
            "Year start", "Year end", "Fit start", "Fit end", "Change", "Std error", "t", "Pr(>|t|)")

          hwrite(contrasts, outPage, border = 0, cellpadding = 2, style = hstyle(contrasts))
        }
      }
    )
  
    
    # status assessment     
           
    hwrite.gap(outPage)
    hwrite("Status assessment", outPage,  br = TRUE, div = TRUE, class = "title")
    
    switch(as.character(info$detGroup),
      Imposex = 
      {
        if(is.na(assessment$summary$clLY))
          hwrite("Insufficient data", outPage, br = TRUE)               
        
        else if (all(is.na(assessment$AC))) 
          hwrite("No assessment criteria", outPage, br = TRUE)
        
        else {

          rv <- with(assessment, {
            AC <- AC[!is.na(AC)]
            sapply(names(AC), USE.NAMES = TRUE, simplify = FALSE, FUN = function(id) {
              out <- summary[c("meanLY", "clLY", id)]
              names(out)[3] <- "AC"
              out
            }) 
          })

          rv <- do.call("rbind", rv)
          
          rv <- within(rv,  p <- ifelse(clLY < AC, "<0.05", ">0.05"))

          id <- c("meanLY", "clLY")
          rv[id] <- lapply(rv[id], hpformat, k = 3)

          id <- "AC"
          rv[id] <- lapply(rv[id], hzap, digits = 4)

          names(rv) <- c("Fitted value", "Upper CL", "Ref value", "Pr(&gt;Ref)")
          
          hwrite(rv, outPage, border = 0, cellpadding = 2, style = hstyle(rv))
        }
        
      },       
      {  
        nyear <- assessment$summary$nyfit
        
        if(is.null(nyear) || nyear < 3)
          hwrite("Insufficient data", outPage, br = TRUE)               
        
        else if (all(is.na(assessment$AC))) 
          hwrite("No assessment criteria", outPage, br = TRUE)
        
        else {

          good.status <- as.character(get.info("determinand", info$determinand, "good.status"))
          ptxt <- switch(good.status, low = "Pr(&gt;t)", high = "Pr(&lt;t)")
          
          is.lognormal <- get.info("determinand", info$determinand, "distribution") %in% "lognormal"
          
          rv <- within(assessment$reference.values, {
            FittedConc <- if (is.lognormal) exp(fit) else fit
            RefConc <- assessment$AC
            tvalue <- difference / se
          })
          
          rv <- rv[!is.na(rv$RefConc), c("FittedConc", "RefConc",  "difference", "se", "tvalue", "p")]

          id <- c("FittedConc", "RefConc")
          rv[id] <- hzap(unlist(rv[id]), digits = 4)

          id <- c("difference", "se")
          rv[id] <- hzap(unlist(rv[id]), digits = 4)

          id <- "tvalue"
          rv[id] <- lapply(rv[id], hformat, k = 3)
          
          id <- "p"
          rv[id] <- lapply(rv[id], hpformat, k = 4)
          
          names(rv) <- c(
            "Conc fitted", "Conc ref", if (is.lognormal) "Log ratio" else "Difference", "Std error", "t", ptxt)
          
          hwrite(rv, outPage, border = 0, cellpadding = 2, style = hstyle(rv))
        }
      })  
        
    closePage(outPage)    
                                              
  })
  invisible()
}




ctsm.web.getKey <- function(info, auxiliary.plot = FALSE, html = FALSE) {

  compartment <- switch(info$compartment, biota = "Biota", sediment = "Sediment", water = "Water")
  
  txt <- paste("Compartment:", compartment)
  
  matrixID <- unique(na.omit(info$matrix))
  matrixID <- sort(matrixID)
  
  if (length(matrixID) == 0) stop('no valid matrix information')
  
  matrixNames <- get.info("matrix", matrixID, "name")
  
  # ad-hoc fix to make names consistent with markdown
  
  matrixNames[matrixNames %in% "erythrocytes (red blood cells in vertebrates)"] <- "red blood cells"
  matrixNames[matrixNames %in% "egg homogenate of yolk and albumin"] <- "egg yolk and albumin"
  matrixNames[matrixNames %in% "hair/fur"] <- "hair"


  txt <- switch(compartment,
    Biota = {
      txt <- paste0(txt, " (", get.info("species", info$species, "common.name"), " ")
      if (length(matrixID) == 1) 
        paste0(txt, matrixNames, ")") 
      else {
        out <- sapply(matrixID, function(i) {
          seriesID <- names(info$matrix)[info$matrix == i]
          detID <- unique(info$determinand[seriesID])
          paste0("(", paste0(detID, collapse = ", "), ")")
        })
        out <- paste(matrixNames, out, sep = " ")
        paste0(txt, paste(out, collapse = "; "), ")")
      }
    },
    Sediment = {
      if (length(matrixID) > 1)
        warning('multiple matrices not supported for sediment in ctsm.web.getKey')
      paste0(txt, " (", matrixNames[1], ")")
    },
    Water = {
      if (length(matrixID) > 1) 
        warning('multiple matrices not supported for water in ctsm.web.getKey')
      paste0(txt, " (", matrixNames[1], ")")
    }
  )

  out <- list(media = txt)

  txt <- paste("Station: ", info$station, sep = "")
  if (!is.na(info$stationName)) txt <- paste(txt, " (", info$stationName, ")", sep = "")

  out$station <- if (html) convert.html.characters(txt) else txt

  out$determinand <- paste("Determinand:", get.info("determinand", info$determinand, "common.name"))


  unitID <- as.character(get.info("determinand", info$determinand, "unit", info$compartment))

  groupID = unique(info$group)
  if (length(groupID) > 1)
    stop('multiple determinand groups not supported in ctsm.web.getKey')

  basis <- as.character(info$basis)
  
  sep.html <- if (html) " " else "~"
  
  if (auxiliary.plot) { 
    start.text <- paste(
      switch(
        info$group, 
        Effects = "Effect", 
        Imposex = "Imposex", 
        "Concentration"
      ),
      "units:", 
      sep = sep.html
    )
  } else {
    start.text <- "Units:"
  }


  # extra.text deals with normalised sediments

  is_extra <- info$compartment == "sediment" && (is.null(info$country) || info$country != "Spain")
  
  if (is_extra) {
    if (html) {
      extra.text <- paste(
        "normalised to", 
        switch(groupID, Metals = "5% aluminium", "2.5% organic carbon")
      )
    } else {
    extra.text <- paste(
        '"normalised to"', 
        switch(groupID, Metals = '"5% aluminium"', '"2.5% organic carbon"'), 
        sep = "~"
      )
    }    
    
    if (length(unitID) > 1L) {
      extra.text <- paste("all", extra.text, sep = sep.html)
    }
    
  }  


  unitText <- mapply(
    label.units, 
    units = unitID, 
    basis = basis, 
    MoreArgs = list(html = html, compartment = info$compartment)
  )
  
  names(unitText) <- info$seriesID
  
  unitID <- unique(unitText)
  
  # if (info$compartment == "sediment" & length(unitID) > 1)
  #   stop("unsupported multiple units for sediment in ctsm.web.getKey")
  
  # out$units <- paste(start.text, do.call("label.units", args), sep = ifelse(html, " ", "~"))

  if (length(unitID) == 1) 
    out$units <- paste(start.text, unitID, sep = sep.html) 
  else {
    wk <- sapply(unitID, function(i) {
      seriesID <- names(unitText)[unitText == i]
      detID <- unique(info$determinand[seriesID])
      paste0('"', "(", paste0(detID, collapse = ", "), ")", '"')
    })
    wk <- paste(unitID, wk, sep = sep.html)
    wk <- paste(wk, collapse = if (html) "; " else paste0('~', '"; "', "~"))
    out$units <- paste(start.text, wk, sep = sep.html)
  }
  
  if (is_extra) {
    out$units <- paste(out$units, extra.text, sep = sep.html)
  }
  
  
  out$extraction <- paste("Data extraction: ", info$extraction, sep = "")

  out
}



plot.data.ylim <- function(...) {

  ylim <- range(..., na.rm = T)
  if (diff(ylim) < 0.00001) 
    ylim + c(-0.1, 0.1) 
  else 
    extendrange(ylim, f = 0.04)
}

plot.data.xlim <- function(...) {

  xrange <- range(..., na.rm = T)
  xlim <- extendrange(xrange, f = 0.04)
  xlim[2] <- min(xlim[2], xrange[2] + 0.99)  # ensures maximum tick mark is last possible data year
  xlim
}

plot.axis <- function(side, ntick.x = 4, ntick.y = 5, xykey.cex = 1, plot.type = c("data", "auxiliary"), 
                      is.data, add.xlab = TRUE, useLogs = TRUE, ...) {

  plot.type <- match.arg(plot.type)

  switch(side,
    bottom = {
      grid.lines(x = c(0, 1), y = c(0, 0), default.units = "npc", gp = gpar(col = "black")) 
      xlim <- current.panel.limits()$xlim
      at <- seq(ceiling(xlim[1]), floor(xlim[2]))
      panel.axis(side = side, outside = TRUE, at = at, labels = FALSE, tck = 0.5, line.col = "black")
      at <- plot.scales(xlim, n = ntick.x)
      add.labels <- add.xlab & current.row() == 1
      panel.axis(
        side = side, outside = TRUE, at = at, labels = add.labels, tck = 1, line.col = "black", 
        text.cex = xykey.cex
      )
    },
    left = {
      grid.lines(x = c(0, 0), y = c(0, 1), default.units = "npc")
      if (!missing(is.data))
        if (!is.data[which.packet()]) return()
      ylim <- current.panel.limits()$ylim
      if (plot.type == "data" | 
          (plot.type == "auxiliary" && names(is.data)[which.packet()] %in% c("value", "concentration"))
      ) {
        tmp <- plot.scales(ylim, n = ntick.y, logData = useLogs)
        at <- if (useLogs) log(tmp) else tmp
      } else {
        tmp <- plot.scales(ylim, n = ntick.y)
        at <- tmp
      }
      panel.axis(
        side = side, outside = TRUE, at = at, labels = format(tmp), tck = 1, line.col = "black", 
        text.cex = xykey.cex
      )
    }
  )
}



plot.AC <- function(AC, ylim, useLogs = TRUE) {

  AC <- AC[!is.na(AC)]
  AC <- sort(if (useLogs) log(AC) else AC)
  AC <- data.frame(id = names(AC), value = AC, ok = AC >= ylim[1] & AC <= ylim[2], stringsAsFactors = FALSE)
  within(AC, pos <- convertY(unit(value, "native"), "npc", valueOnly = TRUE))
}



plot.data <- function(data, assessment, info, type = c("data", "assessment"), 
                      xykey.cex = 1.0, ntick.x = 4, ntick.y = 3, newPage = TRUE, ...) {

  type <- match.arg(type) 

  is.pred <- "pred" %in% names(assessment)
  #   if (is.pred & info$determinand %in% c("VDS", "IMPS", "INTS")) 
  #     assessment$pred <- swap.names(assessment$pred, c("lower", "upper"), 
  #                                   c("ci.lower", "ci.upper"))
  
  is.AC <- !all(is.na(assessment$AC))

  # make data types compatible - i.e. raw data or assessment indices
 
  distribution <- get.info("determinand", info$determinand, "distribution")
  
  if (info$determinand %in% c("MNC")) {
    warning("remember to fix distribution changes")
    distribution <- "normal"
  }
  
  
  useLogs <- distribution %in% "lognormal"

  data <- switch(type, 
    data = 
    {
      out <- subset(data, select = c(year, concentration, qflag))
      if (useLogs) out <- within(out, concentration <- log(concentration))
      out
    },
    assessment = 
    {
      out <- assessment$annualIndex
      names(out)[2] <- "concentration"
      if (info$determinand %in% c("VDS", "IMPS", "INTS")) out$qflag <- rep("", nrow(out))
      out
    }
  )    

  
  if (info$distribution == "beta" && is.pred) {
    assessment$pred <- dplyr::mutate(
      assessment$pred, 
      fit = 100 * plogis(.data$fit),
      ci.lower = 100 * plogis(.data$ci.lower),
      ci.upper = 100 * plogis(.data$ci.upper)
    )
  }
  

  # set up graphical structures
    
  args.list <- list(data$concentration)         # NB have taken logs above, so this is log concentration
  
  if (is.pred) {
    args.list <- c(
      args.list, 
      assessment$pred$ci.lower, 
      assessment$pred$ci.upper
    )
  }
  
  if (info$determinand %in% c("VDS", "IMPS", "INTS")) {
    wk <- info.imposex[
      info.imposex$species %in% info$species & info.imposex$determinand %in% info$determinand, 
      "max_value"]
    ylim <- extendrange(c(0, wk), f = 0.04)
  }            
  else {
    ylim <- c(do.call("plot.data.ylim", args.list))
  }
    
  xlim <- plot.data.xlim(data$year, info$recentYears)

  plot.formula <- data$concentration ~ data$year

  
  # ensures ylabels are formatted correctly and fit into viewport

  ykey <- format(plot.scales(ylim, n = ntick.y, logData = useLogs)) 
  key.ylab.padding <- max(nchar(ykey))


  # sets up viewport so that assessment concentrations and ylabel fit correctly
  
  AC.width <- unit(0, "npc")
  if (is.AC) {

    AC <- plot.AC(assessment$AC, ylim, useLogs)

    # if (any(AC$ok))
    #   AC <- AC[AC$ok,]
    # else if (all(AC$value < ylim[1]))
    #   AC <- tail(AC, 1)
    # else if (all(AC$value > ylim[2]))
    #   AC <- head(AC, 1)
    # else
    # {
    #   wk <- max(which(AC$value < ylim[1]))
    #   AC <- AC[c(wk, wk+1), ]
    # }

    id <- which(AC$ok)
    
    # expand to catch the closest AC below and above the range of the data
    
    if (any(AC$value < ylim[1])) {
      id <- c(id, max(which(AC$value < ylim[1])))
    }

    if (any(AC$value > ylim[2])) {
      id <- c(id, min(which(AC$value > ylim[2])))
    }
    
    AC <- AC[id, ]

    AC.width <- max(unit(rep(xykey.cex, nrow(AC)), "strwidth", as.list(AC$id))) + unit(xykey.cex, "char")
    if (!all(AC$ok)) AC.width <- AC.width + unit(xykey.cex * (0.8 + 0.35), "char")
  }

  wk.viewport <- viewport(
    x = unit(xykey.cex, "char"), y = unit(2 * xykey.cex, "char"), just = c("left", "bottom"),
    width = unit(1, "npc") - AC.width - unit(2 * xykey.cex, "char"), 
    height = unit(1, "npc") - unit(5.5 * xykey.cex, "char"))

  data.plot <- xyplot(
    plot.formula, 
    ylim = ylim, 
    xlim = xlim, 
    xlab = "", 
    ylab = "", 
   	par.settings = list(
   	  axis.line = list(col = "transparent"), 
      layout.widths = list(
        left.padding = 2, axis.left = 0, ylab.axis.padding = 0, ylab = 0, 
        key.ylab.padding = key.ylab.padding * xykey.cex, right.padding = 0, key.right = 0, 
        axis.key.padding = 0, axis.right = 0
      ),
      layout.heights = list(
        axis.bottom = 0, bottom.padding = 2, axis.xlab.padding = 0, xlab = 0, 
        xlab.key.padding = xykey.cex, key.sub.padding = 0, 
        axis.top = 0, top.padding = 0, main = 0, main.key.padding = 0, key.top = 0, key.axis.padding = 0
      )
   	), 
    axis = function(side, ...) plot.axis(side, ntick.x = ntick.x, ntick.y = ntick.y, 
                                         xykey.cex = xykey.cex, useLogs = useLogs, ...), 
    panel = function(x, y)
    {
      plot.panel(x, y, data$qflag, type, AC = assessment$AC, 
                 pred = if (is.pred) assessment$pred else NULL, ylim, useLogs = useLogs,
                 indiCL = assessment$annualIndex)

      if (is.AC) {
        AC <- plot.AC(assessment$AC, ylim, useLogs)        # needs to be before pushViewport
        pushViewport(viewport(clip = "off"))

        if (any(AC$ok)) {
          with(AC, {
            grid.text(
              id[ok], 
              x = unit(1, "npc") + unit(1, "char"), 
              y = pos[ok], 
              just = c("left", "centre"), 
              gp = gpar(cex = xykey.cex)
            )
          })
        } 

        if (any(AC$value < ylim[1])) {
          wk <- max(which(AC$value < ylim[1]))
          id <- AC$id[wk]
          grid.text(
            id, 
            x = unit(1, "npc") + unit(1, "char"), 
            y = 0, 
            just = c("left", "centre"), 
            gp = gpar(cex = xykey.cex)
          )
          grid.lines(
            x = unit(1, "npc") + unit(1.8, "char") + unit(1, "strwidth", id), 
            y = unit.c(unit(0.8, "char"), unit(-0.8, "char")), 
            arrow = arrow(length = unit(0.7, "char")), 
            gp = gpar(cex = xykey.cex)
          )
        }
          
        if (any(AC$value > ylim[2])) {
          wk <- min(which(AC$value > ylim[2]))
          id <- AC$id[wk]
          grid.text(
            id, 
            x = unit(1, "npc") + unit(1, "char"), 
            y = 1, 
            just = c("left", "centre"), 
            gp = gpar(cex = xykey.cex)
          )
          grid.lines(
            x = unit(1, "npc") + unit(1.8, "char") + unit(1, "strwidth", id), 
            y = unit(1, "npc") + unit.c(unit(-0.8, "char"), unit(0.8, "char")), 
            arrow = arrow(length = unit(0.7, "char")), 
            gp = gpar(cex = xykey.cex)
          )
        }
        upViewport()
      }

      pushViewport(viewport(clip = "off"))
      extra <- dplyr::case_when(
        info$group %in% "Effects"               ~ "",
        info$determinand %in% "VDS"             ~ "stage",
        info$determinand %in% c("IMPS", "INTS") ~ "",
        TRUE                                    ~ "concentration"
      )
      ylabel <- paste(
        get.info("determinand", info$determinand, "common.name"), 
        extra
      )
      grid.text(
        ylabel, 0, unit(1, "npc") + unit(1.5, "char"), just = c("left", "bottom"), gp = gpar(cex = xykey.cex))
      upViewport()
    })

  plot.setup(newPage)
  pushViewport(viewport(layout.pos.row = 1))
  pushViewport(wk.viewport)
  print(data.plot, newpage = FALSE)
  upViewport()
  upViewport()
  plot.info(info, ...)
}



plot.setup <- function(newPage) {

  require("lattice")
  require("grid")

  if (newPage) grid.newpage()

  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 4), c("null", "lines")))))
}



label.units <- function(
  units = c("ug/kg", "mg/kg", "ng/ml", "pmol/min/mg protein", "ug/ml", "ug/l", 
            "nmol/min/mg protein", "ng/min/mg protein", "stg", "j/h/g", "mins",
            "d", "%", "nr/1000 cells", "TEQ ug/kg", "ng/l"),
  basis, html = FALSE, compartment, extra.text = NA) {

  units <- match.arg(units)
  
  ok <- basis %in% c("W", "D", "L") |
    (is.na(basis) & units %in% c(
      "ng/ml", "ug/ml", "pmol/min/mg protein", "nmol/min/mg protein", "ng/min/mg protein", "stg", 
      "j/h/g", "mins", "d", "%", "nr/1000 cells"))
  
  if (!ok)
    stop("basis not recognised")
  
    
  if (html)
    units <- switch(
      units, 
      "ug/kg" = "&mu;g kg<sup>-1</sup>", 
      "mg/kg" = "mg kg <sup>-1</sup>",
      "ng/ml" = "ng ml <sup>-1</sup>",
      "ug/ml" = "&mu;g ml<sup>-1</sup>", 
      "ug/l" = "&mu;g l<sup>-1</sup>",
      "ng/l" = "ng l<sup>-1</sup>",
      "TEQ ug/kg" = "TEQ &mu;g kg<sup>-1</sup>", 
      "stg" = "stage",
      "j/h/g" = "J h <sup>-1</sup> g <sup>-1</sup>",
      "pmol/min/mg protein" = "pmol min <sup>-1</sup> mg protein <sup>-1</sup>",
      "nmol/min/mg protein" = "nmol min <sup>-1</sup> mg protein <sup>-1</sup>",
      "ng/min/mg protein" = "ng min <sup>-1</sup> mg protein <sup>-1</sup>", 
      "mins" = "min",
      "d" = "d",
      "%" = "%", 
      "nr/1000 cells" = "nr/1000 cells")
  else
    units <- switch(
      units, 
      "ug/kg" = 'paste(mu, "g") ~ "kg"^-1', 
      "mg/kg" = '"mg kg"^-1',
      "ng/ml" = '"ng ml"^-1',
      "ug/ml" = 'paste(mu, "g") ~ "ml"^-1', 
      "ug/l" = 'paste(mu, "g") ~ "l"^-1', 
      "ng/l" = '"ng l"^-1', 
      "TEQ ug/kg" = 'paste("TEQ", mu, "g") ~ "kg"^-1', 
      "stg" = '"stage"',
      "j/h/g" = '"J h"^-1 ~ "g"^-1',
      "pmol/min/mg protein" = '"pmol min"^-1 ~ "mg protein"^-1',
      "nmol/min/mg protein" = '"nmol min"^-1 ~ "mg protein"^-1',
      "ng/min/mg protein" = '"ng min"^-1 ~ "mg protein"^-1', 
      "mins" = '"min"', 
      "d" = '"d"',
      "%" = '"%"',
      "nr/1000 cells" = '"nr/1000 cells"')
  
  args <- list(units, sep = if (html) " " else "~")
  
  if(!is.na(basis) & compartment != "water") {
    if (html)
      basis <- switch(basis, D = "dry weight", W = "wet weight", L = "lipid weight")
    else
      basis <- switch(basis, D = '"dry weight"', W = '"wet weight"', L = '"lipid weight"')

    args <- c(args, basis)
  }
  
  if (!is.na(extra.text)) 
    args <- c(args, as.list(extra.text))
  do.call("paste", args)
}



plot.panel <- function(
  x, y, qflag, type, AC = NA, pred = NULL, ylim, layout.row, useLogs = TRUE, 
  indiCL) {

  # type
  # data is standard (single panel) data plot
  # assessment is standard (single panel) assessment plot
  # ratio is multipanel ratio plot (under development)
  # multi_assessment is multipanel assessment plot of related compounds
  
  type <- match.arg(type, c("data", "assessment", "ratio", "multi_assessment"))
  
  if (!is.null(pred)) 
    lpolygon(c(pred$year, rev(pred$year)), c(pred$ci.lower, rev(pred$ci.upper)), 
             border = FALSE, col = grey(0.8))
  
  if (!all(is.na(AC))) {
    AC <- if (useLogs) log(na.omit(AC)) else na.omit(AC)
    ok <- AC >= ylim[1] &  AC <= ylim[2]
    if (any(ok)) panel.abline(h = AC[ok], lty = 8, lwd = 0.5)
  }

  if (!is.null(pred)) llines(pred$year, pred$fit, lwd = 2, col = "black")

  if (any(duplicated(x))) x <- jitter(x, amount = 0.1)

  wk.cex <- switch(
    type, 
    data = 2.5, 
    assessment = 2,
    ratio = switch(layout.row, 2.0, 1.4, 0.9, 0.7, 0.6), 
    multi_assessment = switch(layout.row, 2.0, 1.4, 0.9, 0.7, 0.6) 
  )
  
  wk.pch <- switch(
    type, 
    data = "+", 
    assessment = 16, 
    ratio = "+", 
    multi_assessment = 16
  )
  
  wk.cex.qflag = if (wk.pch == "+") wk.cex * 0.8 else wk.cex
  
  # recognise following symbols for qflag: ">" and "?" can arise in ratio plots

  qflag <- as.character(qflag)
  
  is_qflag <- qflag %in% c("<", "D", "Q", ">", "?")

  qflag[is_qflag] <- tolower(qflag[is_qflag])
  
  if (any(is_qflag)) {
    lpoints(
      x[is_qflag], 
      y[is_qflag], 
      pch = qflag[is_qflag], 
      cex = wk.cex.qflag, 
      col = "black"
    )
  }
    
  if (any(!is_qflag)) {
    lpoints(
      x[!is_qflag], 
      y[!is_qflag], 
      pch = wk.pch, 
      cex = wk.cex, 
      col = "black"
    )
  }

  if (!missing(type) && type == "assessment" && "lower" %in% names(indiCL)) {
    with(indiCL, lsegments(year, lower, year, upper, lwd = 2, col = "black"))
  }
}


plot.auxiliary <- function(data, info, auxiliary_id = "default", xykey.cex = 1.0, ntick.x = 3, ntick.y = 3, 
                           newPage = TRUE, ...) {

  # auxiliary_id specifies the choice of 'auxiliary' variables to plot: 
  # default:
  #   sediment = value, concentration, AL, CORG
  #   biota = concentration, LNMEA, DRYWT%, LIPIDWT%
  #   water = not specified yet
  # otherwise must contain four relevant variables
  
  distribution <- get.info("determinand", info$determinand, "distribution")
  
  if (info$determinand %in% "MNC") {
    warning("remember to fix distribution changes")
    distribution <- "normal"
  }
  
  useLogs <- distribution %in% "lognormal"
  

  data <- within(data, {
    if (useLogs) concentration <- log(concentration)
    concentration.qflag <- qflag

    if (exists("concOriginal")) {
      if (useLogs) value <- log(concOriginal)
      value.qflag <- qflagOriginal
    }
  })

  
  stopifnot(
    is.character(auxiliary_id),
    length(auxiliary_id) %in% c(1, 4)
  )
  
  if (length(auxiliary_id) == 1 && auxiliary_id == "default") {
    auxiliary <- switch(
      info$compartment, 
      sediment = c("value", "concentration", "AL", "CORG"), 
      biota = c("concentration", "LNMEA", "DRYWT%", "LIPIDWT%")
    )
  } else {
    auxiliary <- auxiliary_id
  }
    
  auxiliary.qflag <- paste(auxiliary, "qflag", sep = ".")
  
  # not all auxiliary variables have qflags, so create dummy columns
  
  ok <- auxiliary.qflag %in% names(data)
  if (!all(ok))
    data[auxiliary.qflag[!ok]] <- lapply(data[auxiliary[!ok]], function(x) ifelse(is.na(x), NA, ""))
    
  data <- data[c("year", auxiliary, auxiliary.qflag)]


  data <- reshape(
    data, varying = list(auxiliary, auxiliary.qflag), v.names = c("value", "qflag"), direction = "long", 
    timevar = "type", times = auxiliary)
  
  data <- within(data, type <- ordered(type, levels = auxiliary))

  xlim <- range(data$year, info$recentYears)

  is.data <- unlist(with(data, tapply(value, type, function(i) !all(is.na(i)))))
  data <- within(data, value[type %in% names(is.data[!is.data])] <- 0)

  ylim <- with(data, tapply(value, type, extendrange, f = 0.07))        # this is what R will use in xyplot
  if (info$compartment == "sediment")
  {
    ylim$value <- with(subset(data, type %in% c("value", "concentration")), extendrange(value, f = 0.07))
    ylim$concentration <- ylim$value
  }
  
  if (info$determinand %in% c("VDS", "IMPS", "INTS")) 
  {
    wk <- info.imposex[
      info.imposex$species %in% info$species & info.imposex$determinand %in% info$determinand, "max_value"]
    ylim$concentration <- extendrange(c(0, wk), f = 0.07)
  }            
  
  
  ykey <- sapply(levels(data$type), function(i)
  {
    if (is.data[i])
    {  
      wk <- plot.scales(
        ylim[[i]], n = ntick.y, logData = (useLogs & i %in% c("value", "concentration")))
      max(nchar(format(wk)))
    }
    else 0
  })
  key.ylab.padding <- max(ykey[c(1, 3)])
  
  ylim <- with(data, tapply(value, type, range, na.rm = TRUE))          # but this is what must get passed in!
  if (info$compartment == "sediment")
  {
    ylim$value <- with(subset(data, type %in% c("value", "concentration")), range(value, na.rm = TRUE))
    ylim$concentration <- ylim$value
  }

  if (info$determinand %in% c("VDS", "IMPS", "INTS")) 
  {
    wk <- info.imposex[
      info.imposex$species %in% info$species & info.imposex$determinand %in% info$determinand, "max_value"]
    ylim$concentration <- c(0, wk)
  }            
  
  # not perfect, but it does the job without plotting everything in their own viewport
  # first element is tick mark plus padding, then width of key, then a gap for prettiness
  
  between.x <- 1 + max(ykey[c(2, 4)], na.rm = TRUE) * xykey.cex / 2 + 1.5               
  between.y <- 1 + xykey.cex * 2 + 1
  
  wk.viewport <- viewport(
    x = unit(xykey.cex, "char"), 
    y = unit(1 * xykey.cex, "char"), 
    just = c("left", "bottom"),
    width = unit(1, "npc") - unit(2 * xykey.cex, "char"), 
    height = unit(1, "npc") - unit(4 * xykey.cex, "char")
  )

  data.plot <- with(data, {
    xyplot(
      value ~ year | type, 
      xlim = xlim, 
      ylim = ylim,
      xlab = "", 
      ylab = "", 
      scales = list(relation = "free"), 
      strip = FALSE, 
      between = c(list(x = between.x, y = between.y)), 
      par.settings = list(
        axis.line = list(col = "transparent"),
        layout.widths = list(
          left.padding = 2, axis.left = 0, ylab.axis.padding = 0, ylab = 0, 
          key.ylab.padding = key.ylab.padding * xykey.cex, 
          right.padding = 0, key.right = 0, axis.key.padding = 0, axis.right = 0, strip.left = 0, 
          key.left = 0, axis.panel = 0
        ),
        layout.heights = list(
          axis.bottom = 0, bottom.padding = 2, axis.xlab.padding = 0, xlab = 0, 
          xlab.key.padding = xykey.cex, key.sub.padding = 0, 
          axis.top = 0, top.padding = 0, main = 0, main.key.padding = 0, key.top = 0, 
          key.axis.padding = 0, axis.panel = 0)
      ), 
      axis = function(side, ...) 
        plot.axis(side, ntick.x = ntick.x, ntick.y = ntick.y, xykey.cex = xykey.cex, 
                  plot.type = "auxiliary", is.data = is.data, useLogs = useLogs, ...),
      panel = function(x, y, subscripts) {
        type.id <- levels(type)[which.packet()]
        
        if (info$compartment == "sediment" && "country" %in% names(info) && 
            info$country == "Spain" && type.id == "concentration")
          grid.text("data not-normalised", 0.5, 0.5, gp = gpar(cex = xykey.cex))
        else {
          if (any(duplicated(x))) x <- jitter(x, amount = 0.1)
          qflag <- tolower(as.character(qflag[subscripts]))
          qflag <- ifelse(qflag %in% "", "+", qflag)
          lpoints(x, y, pch = qflag, cex = 2.5, col = "black")
        }

        pushViewport(viewport(clip = "off"))
        ylabel <- switch(
          type.id, 
          concentration = switch(
            info$compartment, 
            biota = paste(
              info$determinand, 
              dplyr::case_when(
                info$group %in% "Effects"               ~ "",
                info$determinand %in% "VDS"             ~ "stage",
                info$determinand %in% c("IMPS", "INTS") ~ "",
                TRUE                                    ~ "concentration"
              )
            ),
            sediment = paste(info$determinand, "normalised")
          ),
          value = paste(info$determinand, "non-normalised"),
          LNMEA = {
            family <- as.character(get.info("species", info$species, "family"))
            unit <- get.info("determinand", type.id, "unit", info$compartment)
            paste(
              "Mean ", 
              switch(
                family, 
                Fish = "fish", 
                Bivalvia = "shell", 
                "Gastropoda" = "shell", 
                ""
              ), 
              " length (", unit, ")", 
              sep = ""
            )
          },
          {
            unit <- get.info("determinand", type.id, "unit", info$compartment)
            paste(get.info("determinand", type.id, "common.name"), " (", unit, ")", sep = "")
          }
        )
        grid.text(ylabel, 0, unit(1, "npc") + unit(1, "char"), just = c("left", "bottom"), 
                  gp = gpar(cex = xykey.cex))
        upViewport()
      }
    )
  })

  plot.setup(newPage)
  pushViewport(viewport(layout.pos.row = 1))
  pushViewport(wk.viewport)
  print(data.plot, newpage = FALSE)
  upViewport()
  upViewport()

  plot.info(info, plot.type = "auxiliary", ...)
}


plot.scales <- function(x, n = 5, min.n = 3, logData = FALSE, f = 0.05) {

  # x gives data on log (base e) scale as used in e.g. cstm
  # n is desired number of ticks
  # adapted from axTicks

  x <- x[!is.na(x)]
  if (logData) x <- exp(x)
	rng <- range(x)
	
	small <- .Machine$double.eps^0.5
  if (diff(rng) < small) 
  {
    rng <- 
      if (abs(rng[1]) < small) c(-f, f)
      else rng + c(-f, f) * abs(rng[1])
  }
   

  lin.calc <- !logData | (rng[2] / rng[1] < 10)     # linear scale
  
  if (lin.calc)
  {
    scales <- mean(rng)
    wk.n <- n
    while (length(scales) < min.n)
    {
      scales <- pretty(rng, wk.n, min.n)
      scales <- scales[scales >= rng[1] & scales <= rng[2]]
      wk.n <- wk.n + 1
    }    
    if (n == 3 & length(scales) %in% c(5, 6)) scales <- scales[c(1,3,5)]
    if (n == 3 & length(scales) == 7) scales <- scales[c(1,4,7)]
  	scales.lin <- scales
 	}

  if (!logData) return(scales.lin)

  ii <- c(ceiling(log10(rng[1])), floor(log10(rng[2])))
  scales.log <- lapply(1:3, function(j)
    {
      x10 <- 10^((ii[1] - (j >= 2)):ii[2])
      scales <- switch(j, x10, c(outer(c(1, 3), x10))[-1], c(outer(c(1, 2, 5), x10))[-1])
      scales[scales >= rng[1] & scales <= rng[2]]  
    })

  n.choice <- which.min(abs(sapply(scales.log, length) - n))
  if (length(scales.log[[n.choice]]) < min.n & n.choice < 3) n.choice <- n.choice + 1
  scales.log <- scales.log[[n.choice]]
  if (n == 3 & length(scales.log) %in% c(5, 6)) scales.log <- scales.log[c(1,3,5)]
  
  if (lin.calc && (length(scales.lin) < length(scales.log) | length(scales.log) < n)) scales.lin else scales.log
}



plot.multiassessment <- function(data, assessment, info, ...) {

  is.data <- sapply(assessment, function(i) !is.null(i))
  
  is.pred <- sapply(assessment, function(i) !is.null(i) && !is.null(i$pred))
  is.AC <- sapply(assessment, function(i) !is.null(i) && !all(is.na(i$AC)))
  
  series_distribution <- get.info("determinand", info$determinand, "distribution")
  
  if (any(info$determinand %in% "MNC")) {
    warning("remember to fix distribution changes")
    series_distribution[info$determinand %in% "MNC"] <- "normal"
  }
  
  useLogs <- series_distribution %in% "lognormal"
  
  names(useLogs) <- info$seriesID
  
  
  # make data types compatible - i.e. raw data or assessment indices

  data <- sapply(info$seriesID, simplify = FALSE, function(i) {
    if (is.data[i]) {
      out <- assessment[[i]]$annualIndex
      names(out)[2] <- "concentration"
      out
    }
    else data.frame(year = info$max.year, concentration = 0, qflag = factor("", levels = c("", "<")))
  })  
    

  # transform fitted values for beta distributed data
  
  is.beta <- series_distribution %in% "beta" & is.pred
  
  if (any(is.beta)) {
    assessment[is.beta] <- sapply(assessment[is.beta], simplify = FALSE, function(x) {
      x$pred <- dplyr::mutate(
        x$pred, 
        fit = 100 * plogis(.data$fit),
        ci.lower = 100 * plogis(.data$ci.lower),
        ci.upper = 100 * plogis(.data$ci.upper)
      )
      x
    })
  }
  

  # set up graphical structures

  ylim <- sapply(info$seriesID, simplify = FALSE, function(i) {
    args.list <- list(data[[i]]$concentration)         # NB have taken logs above, so this is log concentration
    if (is.pred[[i]]) args.list <- c(args.list, assessment[[i]]$pred$ci.lower, assessment[[i]]$pred$ci.upper)
    do.call("plot.data.ylim", args.list)
  })

  args.list <- sapply(info$seriesID[is.data], simplify = FALSE, function(i) data[[i]]$year)
  args.list <- c(args.list, info$recentYears)
  xlim <- do.call("plot.data.xlim", args.list)

  plot.formula <- data$concentration ~ data$year


  # ensures ylabels are formatted correctly and fit into viewport

  ntick.y <- ntick.x <- 3
  ykey <- sapply(info$seriesID, simplify = FALSE, function(i) 
    format(plot.scales(ylim[[i]], n = ntick.y, logData = useLogs[i])))
  key.ylab.padding <- max(sapply(ykey, function(i) max(nchar(i))))
  

  ndet <- length(info$seriesID)
  layout.row <- ceiling(sqrt(ndet))
  layout.col <- ceiling(ndet / layout.row)

  add.xlab = 1:ndet <= layout.col
  names(add.xlab) <- info$seriesID

  xykey.cex <- switch(layout.row, 1.4, 1.1, 0.9, 0.7, 0.6)


  # sets up viewport so that assessment concentrations and ylabel fit correctly
  
  xAC <- 0.5

  AC.width <- sapply(info$seriesID[is.data], simplify = FALSE, function(i) {
    
    if (!is.AC[[i]]) return(unit(0, "npc"))

    ylim <- ylim[[i]]
    AC <- plot.AC(assessment[[i]]$AC, ylim, useLogs[i])

    # if (any(AC$ok))
    #   AC <- AC[AC$ok,]
    # else if (all(AC$value < ylim[1]))
    #   AC <- tail(AC, 1)
    # else if (all(AC$value > ylim[2]))
    #   AC <- head(AC, 1)
    # else {
    #   wk <- max(which(AC$value < ylim[1]))
    #   AC <- AC[c(wk, wk+1), ]
    # }
  
    id <- which(AC$ok)
    
    # expand to catch the closest AC below and above the range of the data
    
    if (any(AC$value < ylim[1])) {
      id <- c(id, max(which(AC$value < ylim[1])))
    }
    
    if (any(AC$value > ylim[2])) {
      id <- c(id, min(which(AC$value > ylim[2])))
    }
    
    AC <- AC[id, ]
    
    out <- max(unit(rep(xykey.cex, nrow(AC)), "strwidth", as.list(AC$id))) + 
      unit(xykey.cex, "char")
    if (!all(AC$ok)) out <- out + unit(xykey.cex * (0.8 + 0.35), "char")
    out

  })
  
  AC.width <- do.call("max", AC.width)


  data.plot <- sapply(info$seriesID, simplify = FALSE, function(i) {
    xyplot(
      data[[i]]$concentration ~ data[[i]]$year, 
      ylim = ylim[[i]], 
      xlim = xlim, 
      xlab = "", 
      ylab = "", 
      aspect = 0.7,
     	par.settings = list(
     	  axis.line = list(col = "transparent"), 
     	  layout.widths = list(
     	    left.padding = 2, axis.left = 0, ylab.axis.padding = 0, ylab = 0, 
     	    key.ylab.padding = key.ylab.padding * xykey.cex, right.padding = 0, key.right = 0, 
     	    axis.key.padding = 0, axis.right = 0),
     	  layout.heights = list(
     	    axis.bottom = 0, bottom.padding = 2, axis.xlab.padding = 0, xlab = 0, xlab.key.padding = xykey.cex, 
     	    key.sub.padding = 0, axis.top = 0, top.padding = 0, main = 0, main.key.padding = 0, key.top = 0, 
     	    key.axis.padding = 0)), 
      axis = function(side, ...) {
        plot.axis(
          side, ntick.x = ntick.x, ntick.y = ntick.y, xykey.cex = xykey.cex, 
          is.data = is.data[i], add.xlab = add.xlab[i], useLogs = useLogs[i], ...
        )
      }, 
      panel = function(x, y) {
        if (is.data[i]) {
          plot.panel(
            x, y, data[[i]]$qflag, 
            type = "multi_assessment",
            layout.row = layout.row, 
            AC = assessment[[i]]$AC, 
            pred = if (is.pred[[i]]) assessment[[i]]$pred else NULL, 
            ylim = ylim[[i]], 
            useLogs = useLogs[i], 
            indiCL = assessment[[i]]$data
          )
        }

        if (is.data[i] && is.AC[[i]]) {
          # needs to be before pushViewport (not sure any more following correction to plot.AC)
          AC <- plot.AC(assessment[[i]]$AC, ylim[[i]], useLogs[i])        
          pushViewport(viewport(clip = "off"))

          if (any(AC$ok)) {
            with(AC, grid.text(id[ok], x = unit(1, "npc") + unit(xAC, "char"), y = pos[ok], 
                               just = c("left", "centre"), gp = gpar(cex = xykey.cex)))
          }
          
          if (any(AC$value < ylim[[i]][1])) {
            wk <- max(which(AC$value < ylim[[i]][1]))
            id <- AC$id[wk]
            grid.text(id, x = unit(1, "npc") + unit(xAC, "char"), y = 0, 
                      just = c("left", "centre"), gp = gpar(cex = xykey.cex))
            grid.lines(x = unit(1, "npc") + unit(xAC + 0.8, "char") + unit(1, "strwidth", id), 
                       y = unit.c(unit(0.8, "char"), unit(-0.8, "char")), 
                       arrow = arrow(length = unit(0.7, "char")), gp = gpar(cex = xykey.cex))
          }
            
          if (any(AC$value > ylim[[i]][2])) {
            wk <- min(which(AC$value > ylim[[i]][2]))
            id <- AC$id[wk]
            grid.text(id, x = unit(1, "npc") + unit(xAC, "char"), y = 1, 
                      just = c("left", "centre"), gp = gpar(cex = xykey.cex))
            grid.lines(x = unit(1, "npc") + unit(xAC + 0.8, "char") + unit(1, "strwidth", id), 
                       y = unit(1, "npc") + unit.c(unit(-0.8, "char"), unit(0.8, "char")), 
                       arrow = arrow(length = unit(0.7, "char")), gp = gpar(cex = xykey.cex))
          }
        
          upViewport()
        }

        pushViewport(viewport(clip = "off"))
        grid.text(info$plotNames$assessment[i], 0, unit(1, "npc") + unit(1, "char"), 
                  just = c("left", "bottom"), gp = gpar(cex = xykey.cex))
        upViewport()
      })
    })
#  pushViewport(wk.viewport)


  plot.setup(newPage = TRUE)
  pushViewport(viewport(layout.pos.row = 1))
  
  pushViewport(
    viewport(y = unit(xykey.cex, "char"), just = "bottom", height = unit(1, "npc") - unit(xykey.cex, "char")))
  pushViewport(viewport(layout = grid.layout(layout.row, layout.col)))

  lapply(1:ndet, function(idet) {
    icol <- 1 + (idet - 1) %% layout.col
    irow <- layout.row - (idet - 1) %/% layout.col
    pushViewport(viewport(layout.pos.row = irow, layout.pos.col = icol))
    pushViewport(
      viewport(x = unit(xykey.cex, "char"), y = 0, just = c("left", "bottom"), 
               width = unit(1, "npc") - AC.width - unit((1 + xAC) * xykey.cex, "char"), 
               height = unit(1, "npc") - unit(3 * xykey.cex, "char")))
    print(data.plot[[idet]], newpage = FALSE)
    upViewport()
    upViewport()
   })
   upViewport()
   upViewport()
   upViewport()
   
   plot.info(info, ...)
}



plot.multidata <- function(data, info,  ...) {

  data <- subset(data, !is.na(concentration))

  series_distribution <- get.info("determinand", data$determinand, "distribution")
  
  if (any(data$determinand %in% "MNC")) {
    warning("remember to fix distribution changes")
    series_distribution[data$determinand %in% "MNC"] <- "normal"
  }
  
  useLogs <- series_distribution %in% "lognormal"
  
  data <- within(data, concentration[useLogs] <- log(concentration[useLogs]))
  data <- data[c("year", "sampleID", "seriesID", "qflag", "concentration")]


  data <- reshape(data, direction = "wide", idvar = c("sampleID", "year"), timevar = "seriesID")

  # add in extra rows for recent years to ensure there is a sensible range of years and to cover 
  # situation where a determinand has only a single value and a range can't be calculated - could be done
  # much more elegantly

#  if (max(data$year) < max(info$recentYears) | min(data$year) > min(info$recentYears)) 
#  {
#    last.id <- nrow(data) + c(1, 2)
#    data[last.id,] <- NA
#    data[last.id, "year"] <- range(info$recentYears)
#  }

  # varnames doesn't appear to work in splom, so make sure we give the columns the names
  # we want printed
  
  plot.data <- data[c("year", paste("concentration", info$seriesID, sep = "."))]
  names(plot.data) <- varnames <- c("year", info$plotNames$data)
  
  colID <- c("year", info$seriesID)


  pscales <- lapply(names(plot.data), function(i) {
    if (i == "year") 
      limits <- plot.data.xlim(plot.data[i], info$recentYears) 
    else 
      limits <- plot.data.ylim(plot.data[i])
    
    list(at = NULL, labels = NULL, limits = limits)
  })

  # reduce size of varname text 
  # reduction increases with maximum length of string and number of strings
  
  # n_adj <- length(varnames) %/% 5
  # l_adj <- 0.031 * max(nchar(varnames))
  # varname.cex <- 1 - n_adj * l_adj

  adj <- length(varnames) * max(nchar(varnames))
  varname.cex <- 1 - 0.007 * (adj - 30)
  
    
  data.plot <- splom(~ plot.data, 
    xlab = "", pscales = pscales, 
    varnames = varnames, varname.cex = varname.cex,  
    superpanel = function(z, panel, ...) ctsm.panel.pairs(z, panel = panel, ...),
    par.settings = list(
      layout.heights = list(bottom.padding = 0, axis.bottom = 0, axis.xlab.padding = 0, xlab = 0)),
    panel = function(x, y, i, j) {
      qflag <- if (i > 1) data[paste("qflag", colID[i], sep = ".")] else rep("", length(x))
      lpoints(x, y, col = "black", pch = ifelse(qflag == "", "+", "<"), cex = 1.1)
    }  
  )

  plot.setup(newPage = TRUE)
  pushViewport(viewport(layout.pos.row = 1))
  print(data.plot, newpage = FALSE)
  upViewport()
  plot.info(info, ...)
}


ctsm.panel.pairs <- function (z, panel = lattice.getOption("panel.splom"), lower.panel = panel, 
    upper.panel = panel, diag.panel = "diag.panel.splom", as.matrix = FALSE, 
    groups = NULL, panel.subscripts, subscripts, pscales = 5, 
    prepanel.limits = function(x) if (is.factor(x)) levels(x) else extendrange(range(as.numeric(x), 
        finite = TRUE)), varname.col = add.text$col, varname.cex = add.text$cex, 
    varname.font = add.text$font, varname.fontfamily = add.text$fontfamily, 
    varname.fontface = add.text$fontface, axis.text.col = axis.text$col, 
    axis.text.cex = axis.text$cex, axis.text.font = axis.text$font, 
    axis.text.fontfamily = axis.text$fontfamily, axis.text.fontface = axis.text$fontface, 
    axis.line.col = axis.line$col, axis.line.lty = axis.line$lty, 
    axis.line.lwd = axis.line$lwd, axis.line.alpha = axis.line$alpha, 
    axis.line.tck = 1, ...) 
{
    lower.panel <- if (is.function(lower.panel)) 
        lower.panel
    else if (is.character(lower.panel)) 
        get(lower.panel)
    else eval(lower.panel)
    upper.panel <- if (is.function(upper.panel)) 
        upper.panel
    else if (is.character(upper.panel)) 
        get(upper.panel)
    else eval(upper.panel)
    diag.panel <- if (is.function(diag.panel)) 
        diag.panel
    else if (is.character(diag.panel)) 
        get(diag.panel)
    else eval(diag.panel)
    add.text <- trellis.par.get("add.text")
    axis.line <- trellis.par.get("axis.line")
    axis.text <- trellis.par.get("axis.text")
    n.var <- ncol(z)
    if (n.var == 0) 
        return()
    lim <- vector("list", length = n.var)
    for (i in seq_len(n.var)) lim[[i]] <- if (is.list(pscales) && 
        !is.null(pscales[[i]]$lim)) 
        pscales[[i]]$lim
    else prepanel.limits(z[, i])
    if (length(subscripts)) {
        draw <- is.list(pscales) || (is.numeric(pscales) && pscales != 
            0)
        splom.layout <- grid.layout(nrow = n.var, ncol = n.var)
        pushViewport(viewport(layout = splom.layout, name = "pairs"))
        for (i in 1:n.var) for (j in 1:n.var) {
            if (as.matrix) 
                pushViewport(viewport(layout.pos.row = i, layout.pos.col = j, 
                  name = paste("subpanel", j, i, sep = "."), 
                  clip = trellis.par.get("clip")$panel, xscale = if (is.character(lim[[j]])) 
                    c(0, length(lim[[j]]) + 1)
                  else lim[[j]], yscale = if (is.character(lim[[i]])) 
                    c(0, length(lim[[i]]) + 1)
                  else lim[[i]]))
            else pushViewport(viewport(layout.pos.row = n.var - 
                i + 1, layout.pos.col = j, name = paste("subpanel", 
                j, i, sep = "."), clip = trellis.par.get("clip")$panel, 
                xscale = if (is.character(lim[[j]])) 
                  c(0, length(lim[[j]]) + 1)
                else lim[[j]], yscale = if (is.character(lim[[i]])) 
                  c(0, length(lim[[i]]) + 1)
                else lim[[i]]))
            if (i == j) {
                axls <- if (is.list(pscales) && !is.null(pscales[[i]]$at)) 
                  pscales[[i]]$at
                else if (is.character(lim[[i]])) 
                  seq_along(lim[[i]])
                else pretty(lim[[i]], n = if (is.numeric(pscales)) 
                  pscales
                else 5)
                labels <- if (is.list(pscales) && !is.null(pscales[[i]]$lab)) 
                  pscales[[i]]$lab
                else if (is.character(lim[[i]])) 
                  lim[[i]]
                else rep("", length(axls))
                if (is.numeric(lim[[i]])) {
                  axlims <- range(lim[[i]])
                  axid <- axls > axlims[1] & axls < axlims[2]
                  axls <- axls[axid]
                  labels <- labels[axid]
                }
                diag.panel(x = z[subscripts, j], varname = colnames(z)[i], 
                  limits = lim[[i]], at = axls, lab = labels, 
                  draw = draw, varname.col = varname.col, varname.cex = varname.cex, 
                  varname.font = varname.font, varname.fontfamily = varname.fontfamily, 
                  varname.fontface = varname.fontface, axis.text.col = axis.text.col, 
                  axis.text.cex = axis.text.cex, axis.text.font = axis.text.font, 
                  axis.text.fontfamily = axis.text.fontfamily, 
                  axis.text.fontface = axis.text.fontface, axis.line.col = axis.line.col, 
                  axis.line.lty = axis.line.lty, axis.line.lwd = axis.line.lwd, 
                  axis.line.alpha = axis.line.alpha, axis.line.tck = 0, 
                  ...)
                grid.rect(gp = gpar(col = axis.line.col, lty = axis.line.lty, 
                  lwd = axis.line.lwd, fill = "transparent"))
            }
            else {
                pargs <- if (!panel.subscripts) 
                  c(list(x = z[subscripts, j], y = z[subscripts, 
                    i]), i = i, j = j, list(...))
                else c(list(x = z[subscripts, j], y = z[subscripts, 
                  i], i = i, j = j, groups = groups, subscripts = subscripts), 
                  list(...))
                if (!("..." %in% names(formals(panel)))) 
                  pargs <- pargs[intersect(names(pargs), names(formals(panel)))]
                if (as.matrix) 
                  do.call(if (i > j) 
                    "lower.panel"
                  else "upper.panel", pargs)
                else do.call(if (i < j) 
                  "lower.panel"
                else "upper.panel", pargs)
                grid.rect(gp = gpar(col = axis.line.col, lty = axis.line.lty, 
                  lwd = axis.line.lwd, fill = "transparent"))
            }
            upViewport()
        }
        upViewport()
    }
}



plot.info <- function(info, plot.type = c("data", "auxiliary"), ...) {

  plot.type <- match.arg(plot.type)
  key <- ctsm.web.getKey(info, auxiliary.plot = plot.type == "auxiliary")
  
  pushViewport(viewport(layout.pos.row = 2))
  
  grid.text(key$media, x = unit(1, "char"), y = unit(4, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  grid.text(key$station, x = unit(1, "char"), y = unit(3, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  grid.text(parse(text = key$units), x = unit(1, "char"), y = unit(2, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  grid.text(key$extraction, x = unit(1, "char"), y = unit(1, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  
  invisible()
}



plot.ratio.data <- function(data, numerator, denominator, type = c("logistic", "log")) {

  type <- match.arg(type)
  
  id <- c(
    "year", 
    paste(c("concentration", "qflag"), numerator, sep = "."),
    paste(c("concentration", "qflag"), denominator, sep = ".")
  )

  
  # calculate ratios
    
  data <- data[id]
  names(data) <- c("year", "n_conc", "n_qflag", "d_conc", "d_qflag")
  
  data <- dplyr::mutate(
    data, 
    ratio = switch(
      type,
      logistic = .data$n_conc / (.data$n_conc + .data$d_conc),
      log = .data$n_conc / .data$d_conc
    ), 
    qflag = dplyr::case_when(
      .data$n_qflag %in% "" & .data$d_qflag %in% ""               ~ "+",
      .data$n_qflag %in% "" & .data$d_qflag %in% c("<", "D", "Q") ~ ">",
      .data$n_qflag %in% c("<", "D", "Q") & .data$d_qflag %in% "" ~ "<",
      !is.na(.data$n_qflag) & !is.na(.data$d_qflag)               ~ "?",
      TRUE                                                        ~ NA_character_
    )
  )
  
  data <- data[c("year", "ratio", "qflag")]
  
  data
}


plot.ratio.pred <- function(
  data, 
  type = c("logistic", "log"), 
  control = list(nyear = 5, prop_qflag = 0.1)
) {

  type <- match.arg(type)
  
  data <- na.omit(data)
  
  # only fit smoother if number of years >= 5 and number of 'less-thans' < 10%
  
  nyear <- length(unique(data$year))
  
  prop_qflag <- sum(!data$qflag %in% "+") / length(data$qflag) 
  
  if (nyear < control$nyear || prop_qflag >= control$prop_qflag ) {
    return(NULL)
  }
  
  
  # fit smoother with (optional) random year effect
  # random effect in mgcv requires k <= number of replicated years
  # only fit random effect if 5 or more replicated years

  data$yfac <- factor(data$year)
  
  nyear <- nlevels(data$yfac)

  if (nyear < nrow(data)) {
    rep_year <- data$year[duplicated(data$year)]
    nrep <- length(unique(rep_year))
  } else {
    nrep <- 0
  }
  
  if (nrep < 5) {
    k_choice <- min(10, nyear)
    formula <- ratio ~ s(year, k = k_choice)
  } else {
    k_choice <- min(10, nrep)
    formula <- ratio ~ s(year, k = k_choice) + s(yfac, bs = "re")
  }

  fit <- mgcv::gam(
    formula, 
    data = data, 
    family = switch(type, logistic = "betar", log = "gaussian"), 
    method = "REML"
  )
  
  new_data <- data.frame(
    year = seq(min(data$year), max(data$year)), 
    yfac = factor(min(data$year))
  )
  
  pred <- mgcv::predict.gam(fit, new_data, type = "iterms", se.fit = TRUE)
  
  new_data$fit <- pred$fit[, "s(year)"] + unname(attr(pred, "constant"))
  new_data$se <- pred$se.fit[, "s(year)"]
  new_data <- dplyr::mutate(
    new_data, 
    ci.lower = fit - 2 * se,
    ci.upper = fit + 2 * se
  )
  
  if (type == "logistic") {
    var_id <- c("fit", "ci.lower", "ci.upper")
    new_data[var_id] <- lapply(new_data[var_id], plogis)
  }
  
  new_data
}



plot.ratio <- function(data, info, ...) {
  
  require("lattice")
  require("grid")
  
  # get working data 
  # sediment - use non-normalised concentrations
  # biota - could use values before conversion to target bases, but would need to 
  #   ensure comparability - maybe later
  # don't need to log transform, because looking at ratios - maybe in plots

  if (info$compartment == "sediment") {
    data$concentration <- data$concOriginal
    data$qflag <- data$qflagOriginal
  }
  
  # set up ratios
  
  det_id <- switch(
    info$group, 
    Metals = c("SE", "HG"),
    PAH_parent = c("ANT", "PA", "FLU", "PYR", "ICDP", "BGHIP", "BAA", "CHR"), 
    PBDEs = c("BDE47", "BD153"), 
    Organofluorines = c("PFNA", "PFOA", "PFUNDA", "PFDA", "PFTRDA", "PFDOA"),
    Organochlorines = c("DDEPP", "DDTPP", "DDTOP")
  )
  
  ratio_id <- switch(
    info$group, 
    Metals = "SE / HG",
    PAH_parent = c(
      "ANT / (ANT + PA)", "FLU / (FLU + PYR)", "ICDP / (ICDP + BGHIP)", 
      "BAA / (BAA + CHR)"
    ), 
    PBDEs = "BDE47 / BD153",
    Organofluorines = c("PFNA / PFOA", "PFUNDA / PFDA", "PFTRDA / PFDOA"), 
    Organochlorines = c("DDEPP / DDTPP", "DDTOP / DDTPP")
  )
  
  
  numerator_id <- switch(
    info$group,
    Metals = "SE",
    PAH_parent = c("ANT", "FLU", "ICDP", "BAA"), 
    PBDEs = "BDE47",
    Organofluorines = c("PFNA", "PFUNDA", "PFTRDA"),
    Organochlorines = c("DDEPP", "DDTOP")
  )
  
  denominator_id <- switch(
    info$group,
    Metals = "HG",
    PAH_parent = c("PA", "PYR", "BGHIP", "CHR"),
    PBDEs = "BD153",
    Organofluorines = c("PFOA", "PFDA", "PFDOA"),
    Organochlorines = c("DDTPP", "DDTPP")
  )
  
  ratio_type <- switch(
    info$group,
    PAH_parent = "logistic",
    "log"
  )
  
  use_logs <- ratio_type %in% "log"
  
  ref_lines <- switch(
    info$group, 
    Metals = list(78.96 / 200.59),
    PAH_parent = list(0.1, c(0.4, 0.5), c(0.2, 0.5), c(0.2, 0.35)),
    PBDEs = list(NA), 
    Organofluorines = list(NA, NA, NA),
    Organochlorines = list(1, c(0.3, 1))
  )
  
  ref_txt <- switch(
    info$group, 
    Metals = list(c("no mediation", "possible mediation")),
    PAH_parent = list(
      c("petrogenic", "pyrolitic"), 
      c("petrogenic", "oil combustion", "coal combustion"),
      c("petrogenic", "oil combustion", "coal combustion"),
      c("petrogenic", "coal combustion", "combustion")
    ), 
    PBDEs = list(NA),
    Organofluorines = list(NA, NA, NA),
    Organochlorines = list(
      c("new contamination", "old contamination"), 
      c("technical DDT", "", "technical dicofol")
    )
  )
  
  names(ref_lines) <- names(ref_txt) <- ratio_id  
  
  
  # restrict to relevant determinands   
  
  data <- data[data$determinand %in% det_id, ]
  
  
  # identify determinands in det_id which are not reported with the data
  
  missing_det <- setdiff(det_id, unique(as.character(data$determinand)))
  
  
  # widen data 
  
  data <- data[c("year", "sampleID", "determinand", "qflag", "concentration")]
  
  data <- reshape(
    data, 
    direction = "wide", 
    idvar = c("sampleID", "year"), 
    timevar = "determinand"
  )

  
  # add in dummy columns to deal with variables that are not reported
  
  if (length(missing_det) > 0) {
    new_id <- paste("qflag", missing_det, sep = ".")
    data[new_id] <- rep(NA_character_, nrow(data))
    
    new_id <- paste("concentration", missing_det, sep = ".")
    data[new_id] <- rep(NA_real_, nrow(data))
  }
  

  # calculate ratios 
  # returns NULL if any ratio has no data (either because a variable was not
  #   reported, or because no samples have both variables reported)
  
  data <- mapply(
    numerator = numerator_id, 
    denominator = denominator_id, 
    FUN = plot.ratio.data,
    MoreArgs = list(data = data, type = ratio_type), 
    SIMPLIFY = FALSE
  )
  
  is_data <- sapply(data, function(x) !is.null(x) && any(!is.na(x$ratio)))
  
  
  # plot logistic ratios on raw scale, log ratios on log scale
  
  if (ratio_type == "log") {
    data[is_data] <- sapply(data[is_data], simplify = FALSE, FUN = function(x) {
      x$ratio <- log(x$ratio)
      x
    })
  }
  
  
  pred <- lapply(data, plot.ratio.pred, type = ratio_type)
  
  names(data) <- names(pred) <- names(is_data) <- ratio_id  
  
  
  # set up graphical structures
  
  ylim <- sapply(ratio_id, simplify = FALSE, function(i) {
    if (!is_data[i]) {
      out <- switch(ratio_type, logistic = c(0, 1), log = c(0.5, 2))
      return(out)
    }
    args.list <- list(data[[i]]$ratio)         
    if (!is.null(pred[[i]])) 
      args.list <- c(args.list, pred[[i]]$ci.lower, pred[[i]]$ci.upper)
    do.call("plot.data.ylim", args.list)
  })
  
  args.list <- sapply(ratio_id, simplify = FALSE, function(i) data[[i]]$year)
  args.list <- c(args.list, info$recentYears)
  xlim <- do.call("plot.data.xlim", args.list)
  
  plot.formula <- data$ratio ~ data$year
  
  
  # ensures ylabels are formatted correctly and fit into viewport
  
  ntick.y <- ntick.x <- 3
  ykey <- sapply(ratio_id, simplify = FALSE, function(i) 
    format(plot.scales(ylim[[i]], n = ntick.y, logData = FALSE)))
  key.ylab.padding <- max(sapply(ykey, function(i) max(nchar(i))))
  
  
  ndet <- length(ratio_id)
  layout.row <- ceiling(sqrt(ndet))
  layout.col <- ceiling(ndet / layout.row)
  
  add.xlab = 1:ndet <= layout.col
  names(add.xlab) <- ratio_id
  
  xykey.cex <- switch(layout.row, 1.4, 1.1, 0.9, 0.7, 0.6)
  ref.cex <- switch(layout.row, 1.2, 1.0)
  
  # get positions for plotting reference text  

  ref_txt <- mapply(
    ref_txt, ref_lines, ylim, 
    FUN = function(txt, values, ylim) {

      if (any(is.na(values))) {
        return(NA)
      }
      
      if (use_logs) {
        values = log(values)
      }
      
      lower <- c(-Inf, values)
      upper <- c(values, Inf)
      
      ok <- upper > ylim[1] & lower < ylim[2]
      
      lower <- pmax(lower, ylim[1])
      upper <- pmin(upper, ylim[2])
      
      out <- (lower + upper) / 2
      out[!ok] <- NA
      
      names(out) <- txt
      out
    },
    SIMPLIFY = FALSE
  )
 
  
  # sets up viewport so that assessment concentrations and ylabel fit correctly
  
  xAC <- 0.5
  
  # AC.width <- unit(0, "npc")
  
  AC.width <- sapply(ratio_id, simplify = FALSE, function(i) {
    
    if (!is_data[i] || all(is.na(ref_txt[[i]]))) {
      return(unit(0, "npc"))
    }
    
    txt <- ref_txt[[i]]
    
    max(unit(rep(ref.cex, length(txt)), "strwidth", as.list(names(txt)))) + 
      unit(ref.cex, "char")
  })
  
  AC.width <- do.call("max", AC.width)

  
  data.plot <- sapply(ratio_id, simplify = FALSE, function(i) {
    xyplot(
      ratio ~ year, 
      data = data[[i]],
      ylim = ylim[[i]], 
      xlim = xlim, 
      xlab = "", 
      ylab = "", 
      aspect = 0.7,
      par.settings = list(
        axis.line = list(col = "transparent"), 
        layout.widths = list(
          left.padding = 2, axis.left = 0, ylab.axis.padding = 0, ylab = 0, 
          key.ylab.padding = key.ylab.padding * xykey.cex, right.padding = 0, 
          key.right = 0, axis.key.padding = 0, axis.right = 0
        ),
        layout.heights = list(
          axis.bottom = 0, bottom.padding = 2, axis.xlab.padding = 0, xlab = 0, 
          xlab.key.padding = xykey.cex, key.sub.padding = 0, axis.top = 0, 
          top.padding = 0, main = 0, main.key.padding = 0, key.top = 0, 
          key.axis.padding = 0
        )
      ), 
      axis = function(side, ...) {
        plot.axis(
          side, ntick.x = ntick.x, ntick.y = ntick.y, xykey.cex = xykey.cex, 
          is.data = is_data[i], add.xlab = add.xlab[i], useLogs = use_logs, ...
        )
      }, 
      panel = function(x, y) {
        if (is_data[i]) {
          plot.panel(
            x, y, data[[i]]$qflag,
            type = "ratio",
            AC = ref_lines[[i]],
            pred = pred[[i]],
            layout.row = layout.row,
            ylim = ylim[[i]],
            useLogs = use_logs
          )
        }
        
        if (is_data[i] && !all(is.na(ref_txt[[i]]))) {
          
          AC <- plot.AC(ref_txt[[i]], ylim[[i]], useLogs = FALSE)        
          pushViewport(viewport(clip = "off"))
          
          with(
            AC, 
            grid.text(
              id, x = unit(1, "npc") + unit(xAC, "char"), y = pos, 
              just = c("left", "centre"), gp = gpar(cex = ref.cex)
            ) 
          )
          
          upViewport()
        }
        
        pushViewport(viewport(clip = "off"))
        grid.text(i, 0, unit(1, "npc") + unit(1, "char"), 
                  just = c("left", "bottom"), gp = gpar(cex = xykey.cex))
        upViewport()
      })
  })
  #  pushViewport(wk.viewport)
  
  plot.setup(newPage = TRUE)
  pushViewport(viewport(layout.pos.row = 1))
  
  pushViewport(
    viewport(
      y = unit(xykey.cex, "char"), 
      just = "bottom", 
      height = unit(1, "npc") - unit(xykey.cex, "char")
    )
  )
  
  pushViewport(viewport(layout = grid.layout(layout.row, layout.col)))
  
  lapply(1:ndet, function(idet) {
    icol <- 1 + (idet - 1) %% layout.col
    irow <- layout.row - (idet - 1) %/% layout.col
    pushViewport(viewport(layout.pos.row = irow, layout.pos.col = icol))
    pushViewport(
      viewport(x = unit(xykey.cex, "char"), y = 0, just = c("left", "bottom"), 
               width = unit(1, "npc") - AC.width - unit((1 + xAC) * xykey.cex, "char"), 
               height = unit(1, "npc") - unit(3 * xykey.cex, "char")))
    print(data.plot[[idet]], newpage = FALSE)
    upViewport()
    upViewport()
  })
  upViewport()
  upViewport()
  upViewport()
  
  
  # hack of plot.info to allow for dimensionless units
  
  key <- ctsm.web.getKey(info, auxiliary.plot = FALSE)
  
  plot_key <- function(txt, y_lines) {
    grid.text(
      txt, 
      x = unit(1, "char"), 
      y = unit(y_lines, "lines"), 
      gp = gpar(cex = 0.8), 
      just = c("left", "centre")
    )
  }
  
  pushViewport(viewport(layout.pos.row = 2))
  
  plot_key(key$media, 4) 
  plot_key(key$station, 3) 
  plot_key("Units: dimensionless", 2) 
  plot_key(key$extraction, 1)
  
  invisible()
}
