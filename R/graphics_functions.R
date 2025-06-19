# Graphical routines ----

#' Graphical summaries of an assessment
#'
#' Generates a series of assessment plots for each time series. The plots are 
#' exported as either png or pdf files.
#'
#' @param assessment_obj An assessment object resulting from a call to
#'   run_assessment
#' @param subset An optional vector specifying which timeseries are to be
#'   plotted. An expression will be evaluated in the timeSeries component of
#'   assessment_obj; use `series` to identify individual timeseries.
#' @param output_dir The output directory for the assessment plots (possibly
#'   supplied using `file.path`). The default is the working directory. The
#'   output directory must already exist.
#' @param file_type A character vector specifying the types of assessment plot. 
#'   The default `c("data", "index", "auxiliary")` produces three plots for 
#'   each time series. See details 
#' @param file_format A character string specifying Whether the files should be 
#'   png (the default) or pdf.
#' @param auxiliary A character string specifying the auxiliary variables 
#'   plotted if `file_type = "auxiliary"`. See details
#'
#' @returns A series of png or pdf files with graphical summaries of an
#'   assessment. 
#'   
#' @details
#' 
#'   ## Types of assessment plots 
#' 
#'   * `file_type = "data"` shows the raw data with the fitted trend and 
#'   pointwise two-sided 90% confidence bands   
#'   * `file_type = "index"` shows annual indices that summarise the data for 
#'   each year with the fitted trend and pointwise two-sided 90% confidence 
#'   bands   
#'   * `file_type = "auxiliary"` shows the raw data and key auxiliary variables 
#'   see below)
#'   
#'   ## Auxiliary variables
#' 
#'   The default (`auxiliary = "default"`) is to plot the following variables:  
#'   
#'   * biota: determinand concentration, LNMEA (mean length), DRYWT% (dry weight 
#'   content), LIPIDWT% (lipi weight content)
#'   * sediment: non-normalised determinand concentration, normalised 
#'   determinand concentration, AL (aluminium concentration), CORG (organic 
#'   carbon content)
#'   * water: no plots are generated at present
#'   
#'   For biota, the determinand concentration will always be plotted, but it is 
#'   possible to change the three auxiliary variables. For example, to plot 
#'   WTMEA (mean weight) instead of LIPIDWT% you would set `auxiliary = 
#'   c("LNMEA", "WTMEA", "DRYWT%)`. For this to work, WTMEA must previously have
#'   been specified as an auxiliary variable for the determinand in question 
#'   using the `biota_auxliary` column in the determinand reference table. At 
#'   present, there must always be three auxiliary variables for biota.
#'   
#'   For sediment, the non-normalised determinand concentration and the 
#'   normalised determinand concentration will always be plotted, but it is 
#'   possible to change the two auxiliary variables. For example, for metals in 
#'   sediment, you might set `auxiliary = c("AL", "LI")` to plot aluminium
#'   and lithium concentrations instead of aluminium and organic carbon 
#'   concentrations. Again, for this to work, LI must previously have been 
#'   specified as an auxiliary variable for the determinand in question using 
#'   the `sediment_auxliary` column in the determinand reference table. 
#'   At present, there must always be two auxiliary variables for sediment.
#'   
#'   At present, plots for only a limited range of auxiliary variables are 
#'   supported. More flexibility in these plots, such as changing the number of 
#'   auxiliary variables, is desirable and will emerge in due course.
#'
#' @export
plot_assessment <- function(
    assessment_obj, 
    subset = NULL, 
    output_dir = ".",
    file_type = c("data", "index", "auxiliary"),
    file_format = c("png", "pdf"), 
    auxiliary = "default") {

  # silence non-standard evaluation warnings
  .data <- NULL
  
  # graphics_functions.R

  # check file_type, file_format and output_dir are valid
  
  file_format = match.arg(file_format)
  

  if (!all(file_type %in% c("data", "index", "auxiliary"))) {
    stop(
      "\nArgument 'file_type' is invalid: ", 
      "must be 'data','index', 'auxiliary' or some combination of these", 
      call. = FALSE
    )
  }
  
  if ("auxiliary" %in% file_type && assessment_obj$info$compartment == "water") {
    warning(
      "Auxiliary plots are not currently generated for water assessments",
      call. = FALSE
    )
    file_type <- setdiff(file_type, "auxiliary")
    
    if (length(file_type) == 0L) {
      return(invisible())
    }
  }  

  
  if (!dir.exists(output_dir)) {
    stop(
      "\nThe output directory '", output_dir, "' does not exist.\n", 
      "Create it or check the information supplied to argument 'output_dir'",
      " is correct.",
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


  # ad-hoc fix to deal with EROD labelling (which is sex-specific)
  # have raised an issue to deal with this after next release by making the sex
  # classification part of the subseries column
  
  if (info$compartment == "biota") {
    timeSeries <- dplyr::mutate(
      timeSeries, 
      subseries = dplyr::case_when(
        .data$sex == "M" ~ "males",
        .data$sex == "F" ~ "females",
        .default = .data$subseries
      )
    )
  }  

  
  # get relevant subset 
  
  timeSeries <- apply_subset(timeSeries, subset, parent.frame())

  if (nrow(timeSeries) == 0L) {
    warning(
      "no timeseries were selected - nothing has been plotted"
    )
    return(invisible())
  }
    
  series_id <- row.names(timeSeries)


  # identify variables for auxiliary plots - only applies to biota and 
  # sediment at present
  
  if ("auxiliary" %in% file_type) {

    if (length(auxiliary) == 1 && auxiliary == "default") {
      auxiliary <- switch(
        info$compartment, 
        sediment = c("AL", "CORG"), 
        biota = c("LNMEA", "DRYWT%", "LIPIDWT%")
      )
    } 
    
    ok_length <- switch(info$compartment, sediment = 2, biota = 3)

    if (length(auxiliary) != ok_length) {
      stop(
        "\nArgument 'auxiliary' is invalid: ", 
        "must be either 'default' or a length-", ok_length, 
        " character vector\nspecifying the auxiliary variables to be plotted",  
        call. = FALSE
      )
    }

    ok <- auxiliary %in% names(assessment_obj$data)    
    if (!all(ok)) {
      stop(
        "\nArgument 'auxiliary' is invalid: ",
        "the following variables are not in the data:\n", 
        paste(auxiliary[!ok], collapse = ", "),   
        call. = FALSE
      )
    }
    
    auxiliary_id <- switch(
      info$compartment,
      sediment = c("value", "concentration", auxiliary),
      biota = c("concentration", auxiliary)
    ) 

  }

  
  # plot each timeSeries
  
  lapply(series_id, function(id) {
    
    data <- dplyr::filter(assessment_obj$data, .data$seriesID == id)
    
    assessment <- assessment_obj$assessment[[id]]


    # get relevant series info
        
    series <- timeSeries[id, ]
    

    # get file name from id, and add country and station name 
    # for easier identification
    
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
    
    
    # plot assessment with index
    
    if ("index" %in% file_type) {
      
      output_file <- paste0(output_id, " index.", file_format)
      output_file <- file.path(output_dir, output_file)
      
      switch(
        file_format, 
        png = png(output_file, width = 680, height = 480), 
        pdf = pdf(output_file, width = 7, height = 7 * 12 / 17)
      )
      
      plot_data(data, assessment, series, info, type = "index", xykey.cex = 1.4) 
      dev.off()
      
    }    
    
    
    # plot assessment with data
    
    if ("data" %in% file_type) {
      
      output_file <- paste0(output_id, " data.", file_format)
      output_file <- file.path(output_dir, output_file)
      
      switch(
        file_format, 
        png = png(output_file, width = 680, height = 480), 
        pdf = pdf(output_file, width = 7, height = 7 * 12 / 17)
      )
      
      plot_data(data, assessment, series, info, type = "data", xykey.cex = 1.4)
      dev.off()
      
    }  

    
    # plot data with auxiliary variables
    
    if ("auxiliary" %in% file_type) {
      
      output_file <- paste0(output_id, " auxiliary.", file_format)
      output_file <- file.path(output_dir, output_file)
      
      switch(
        file_format, 
        png = png(output_file, width = 680, height = 480), 
        pdf = pdf(output_file, width = 7, height = 7 * 12 / 17)
      )
     
      plot_auxiliary(data, assessment, series, info, auxiliary_id, xykey.cex = 1.4) 
      dev.off()
      
    }    
    
  })
  
  invisible() 
}  


#' Apply subsetting to a time series
#' @description 
#' This is an internal function that applies a subsetting function to a timeseries.
#' It is somewhat complex, due to way it is designed to be called. The subset is designed
#' to be passed either as a value (`NULL` or  logical), as a expression, or as a variable
#' holding an expression. When it is an expression, it is applied in the context of the
#' timeseries data frame to generate a vector of booleans for subsetting. The complexity
#' comes from the way this implements lazy evaluation, so we cannot evaluate the expression
#' in the normal calling context.
#' 
#' The function also removes row names, and converts the series column to row names.
#' 
#' @param timeSeries a time series data frame
#' @param subset (default `NULL`) either `NULL`, which selects all entries, or a logical,
#'   which `TRUE` selects all entries and `FALSE` none of them, or an expression, which
#'   will be evaluated in the context of the dataframe to generate a subsetting vector.
#' @param env the calling environment for variable values.
apply_subset <- function(timeSeries, subset, env = parent.frame()) {
  subsetExpression <- substitute(subset, parent.frame())
  subsetExpressionType <- typeof(subsetExpression)
  if (subsetExpressionType == 'NULL') {
    ## Do nothing
  } else if (subsetExpressionType == 'logical') {
    timeSeries <- timeSeries[subset, ]
    row.names(timeSeries) <- NULL
  } else if (subsetExpressionType == 'symbol' || subsetExpressionType == 'language') {
    ok <- eval(subsetExpression, timeSeries, env)
    if (is.null(ok)) {
      ok <- TRUE
    } else if (is.expression(ok)) {
      ok <- eval(ok, timeSeries, env)
    }
    timeSeries <- timeSeries[ok, ]
    row.names(timeSeries) <- NULL
  }
  
  timeSeries <- tibble::column_to_rownames(timeSeries, "series")
  return(timeSeries)
}


ctsm.format <- function(x, y = x, nsig = 3) {

  # get arguments to format y to nsig significant figures
  # and use these to format x
  
  wk <- floor(log10(abs(y))) + 1 # digits to left of dp
  digits <- max(wk, nsig)
  nsmall <- max(0, nsig - wk)
  format(x, digits = digits, nsmall = nsmall, scientific = nsig)
}


#' Make key for annotating assessment plots
#' 
#' Provides the meta-data that accompanies each plot in `plot_assessment`: 
#' compartment, station, units and data extraction. It also provides extended 
#' information for plots involving multiple timeseries e.g. in `plot_multidata`, 
#' `plot_multiassessment` and `plot_ratio`. 
#'
#' @param series List describing the timeseries in the plot. If there are 
#' multiple timeseries, then some entries of the list require multiple entries, 
#' one for each timeseries. See details.
#' @param info The standard `info` element of the harsat object.
#' @param type Character scalar (default `"data"`) describing the type of plot:
#' other options are `"index"`, `"auxiliary"`, `"data_splom"`, `"index_mp"`, 
#' `"ratio_mp"`. 
#'
#' @details
#' `series` is a list with giving information about a single timeseries 
#' (when called by `plot_assessment`) or multiple related timeseries (when 
#' called by `plot_multidata`, `plot_multiassessment`, and `plot_ratio`). It
#' must contain the following (compartment-specific) elements:
#' * `determinand` (all): a character vector giving the 
#' determinand for each timeseries; if there are multiple timeseries, then the 
#' vector must have names which are the corresponding series identifiers 
#' * `matrix` (biota, sediment): a character vector giving the matrix for 
#' each timeseries; if there are multiple timeseries, then the vector must have
#' names which are the corresponding series identifiers; 
#' the matrices might vary across timeseries, but this is currently likely only 
#' for biological effects 
#' * `filtration` (water): a character vector giving the filtration method for 
#' each timeseries; if there are multiple timeseries, then the vector must have
#' names which are the corresponding series identifiers; 
#' the filtration methods might vary across timeseries, but this is currently 
#' unlikely in practice
#' * `basis` (all): a character vector giving the basis for 
#' each timeseries; if there are multiple timeseries, then the vector must have
#' names which are the corresponding series identifiers; 
#' the bases might vary across timeseries, but this is currently likely only for 
#' biological effects
#' * `species` (biota): a character giving the species (common to all 
#' timeseries)
#' * `subseries` (all): a character giving the subseries (common to all 
#' timeseries); can be `NA_character_`
#' 
#' In addition, if `info$normalise` is `TRUE`, then it must also contain:
#' * `normaliser`, `normaliser_value` and `normaliser_unit`: character, numeric,
#' and character vectors giving the normaliser, it's value and the unit on which
#' it is expressed for each timeseries; missing values are allowed if there is 
#' no normalisation for the timeseries involved
#'
#' @return a list with four text (or expression) elements:
#' * `media`: compartment , species and matrix (biota), matrix (sediment), 
#' filtration (water), subseries
#' * `station`: station_name and station_longname
#' * `units`: an expression giving the units including any normalisation
#' * `extraction`: the extraction date
#' 
#' Both `media` and `units` can be complicated if there are multiple 
#' timeseries, particularly if e.g. there are different normalisers (sediment) 
#' or matrices (biological effects) involved 
#'    `
make_key <- function(series, info, type) {

  # type only affects the bahaviour of the unit element of the key (which 
  # gets very messy)
  
  type <- match.arg(
    type, 
    c("data", "index", "auxiliary", "data_splom", "index_mp", "ratio_mp")
  )

  
  # media
  
  # This combines: 
  # - compartment  
  # - species (biota)
  # - matrix (biota and sediment); filtration (water) 
  # - subseries
  
  compartment <- stringr::str_to_sentence(info$compartment)
  
  txt <- paste0("Compartment: ", compartment, " (") 

  if (info$compartment == "biota") {
    species <- ctsm_get_info(info$species, series$species, "common_name")
    txt <- paste0(txt, species, " ")
  }
  
  # for convenience in code below, think of filtration as a matrix

  matrix_all <- switch(
    compartment,
    Biota = series$matrix,
    Sediment = series$matrix,
    Water = series$filtration
  )
  
  matrix_id <- sort(unique(matrix_all))

  matrix_name <- switch(
    compartment,
    Biota = ctsm_get_info(info$matrix, matrix_id, "name"),
    Sediment = ctsm_get_info(info$matrix, matrix_id, "name"),
    Water = matrix_id
  )

  # ad-hoc fix to make names consistent with markdown
  
  matrix_name[matrix_name %in% "erythrocytes (red blood cells in vertebrates)"] <- 
    "red blood cells"
  matrix_name[matrix_name %in% "egg homogenate of yolk and albumin"] <- 
    "egg yolk and albumin"
  matrix_name[matrix_name %in% "hair/fur"] <- "hair"
  
  if (length(matrix_id) == 1L) { 
    txt <- paste0(txt, matrix_name) 
  } else {
    out <- sapply(matrix_id, function(i) {
      series_id <- names(matrix_all)[matrix_all == i]
      det_id <- unique(series$determinand[series_id])
      det_id <- paste0(det_id, collapse = ", ")
      paste0("(", det_id, ")")
    })
    out <- paste(matrix_name, out, sep = " ")
    txt <- paste0(txt, paste(out, collapse = "; "))
  }

  if (is.na(series$subseries)) {
    txt <- paste0(txt, ")")
  } else {
    txt <- paste0(txt, " - ", series$subseries, ")")
  }
    
  out <- list(media = txt)

  
  # station
  
  txt <- paste("Station: ", series$station_name, sep = "")
  if (!is.na(series$station_longname) && 
      series$station_name != series$station_longname) {
    txt <- paste(txt, " (", series$station_longname, ")", sep = "")
  }
  
  out$station <- txt
  
  
  # unit
  # need an expression for sending to grid graphics to deal with superscripts 
  # and greek characters (e.g. mg/kg)
  # text should be separated by a tilde
  # however need to put % in quotes because it is a special character!
  # this also applies to determinands which might e.g. have a % in their name
  # things get very messy because of this - I'm sure it can be done better!
  
  # there are many special cases!
  # for water, the basis is just assumed to be wet weight

  
  if (type == "ratio_mp") {
    
    out$unit = "units:~dimensionless"
  
  } else {
    
    # set up basic structures used by everything
    
    unit_id <- ctsm_get_info(
      info$determinand, 
      series$determinand, 
      "unit", 
      info$compartment, 
      sep = "_"
    )

    basis <- if (info$compartment == "water") {
      NA_character_ 
    } else {
      series$basis
    }

    unit_txt <- label_unit(unit_id, basis)
    
    names(unit_txt) <- series$seriesID
    
    
    is_normalised <- info$normalise && any(!is.na(series$normaliser))
    
    if (is_normalised) {
      normaliser_txt <- label_normaliser( 
        series$normaliser, 
        series$normaliser_value, 
        series$normaliser_unit
      )
      
      names(normaliser_txt) <- series$seriesID
    }
    

    if (type %in% c("data", "index")) {

      out$unit <- paste("Units:", unit_txt, sep = "~")
      
      if (is_normalised) {
        out$unit <- paste(out$unit, "normalised~to", normaliser_txt, sep = "~")
      }
      
    } else if (type %in% "auxiliary") {
    
      group_id <- ctsm_get_info(
        info$determinand, 
        series$determinand, 
        "group", 
        info$compartment,
        sep = "_"
      )
      
      start_txt <- switch(
        group_id, 
        Effects = "Effect~units:", 
        Imposex = "Imposex~units:", 
        "Concentration~units:"
      )
      
      out$unit <- paste(start_txt, unit_txt, sep = "~")
      
      if (is_normalised) {
        
        if (info$compartment %in% c("biota", "water")) {
          out$unit <- paste(out$unit, "normalised~to", normaliser_txt, sep = "~")
        } else {
          n <- nchar(normaliser_txt)
          normaliser_txt <- paste0(substring(normaliser_txt, 1, n-1), ')"')
          out$unit <- paste(
            out$unit, 
            '"(normalised to"',
            normaliser_txt, 
            sep = "~"
          )
        }
      } 
      
    } else {
      
      # multipanel index plots or splom data plots where all sorts of things can
      # happen
      
      n_unit <- length(unique(unit_txt))
      
      n_normaliser <- 
        if (is_normalised) {
          length(unique(normaliser_txt))
        } else {
          0L
        }
      
      
      if (n_unit == 1L) {
        
        unit_txt <- unique(unit_txt)
        
        out$unit <- paste("Units:", unit_txt, sep = "~")
        
        if (n_normaliser == 1L) {
          out$unit <- paste(out$unit, "normalised~to", normaliser_txt, sep = "~")
        }
        
        if (n_normaliser > 1L) {
          
          normaliser_unique <- sort(unique(normaliser_txt), na.last = TRUE)
          
          # not everything might be normalised - sorting ensures the non-
          # normalised items are listed last
          
          if (any(is.na(normaliser_txt))) {
            normaliser_txt[is.na(normaliser_txt)] <- '"non-normalised"'
            normaliser_unique[is.na(normaliser_unique)] <- '"non-normalised"'
          }  
          
          txt <- sapply(normaliser_unique, function(i) {
            series_id <- names(normaliser_txt)[normaliser_txt == i]
            det_id <- unique(series$determinand[series_id])
            det_id <- paste0(det_id, collapse = ", ")
            paste0("(", det_id, ")")
          })
          
          # can't simplify text below using collapse = ",~" because the , is  
          # then adjacent to a quoted text leading to strange behaviour
          
          txt <- paste0(txt, c(rep(",", n_normaliser - 2), " or", ""))
          txt <- paste0('"', txt, '"')
          txt <- paste(normaliser_unique, txt, sep = "~")
          txt <- paste(txt, collapse = "~")
          
          out$unit <- paste(out$unit, "normalised~to", txt, sep = "~")
        }
        
      } else {
        
        if (is_normalised) {
          # not everything might be normalised
          ok <- !is.na(normaliser_txt)
          unit_txt[ok] <- paste(
            unit_txt[ok], 
            "normalised~to", 
            normaliser_txt[ok], 
            sep = "~"
          )
        }
        
        unit_unique <- sort(unique(unit_txt))
        
        txt <- sapply(unit_unique, function(i) {
          series_id <- names(unit_txt)[unit_txt == i]
          det_id <- unique(series$determinand[series_id])
          det_id <- paste0(det_id, collapse = ", ")
          paste0("(", det_id, ")")
        })
        
        # can't simplify text below by using collapse = ";~" because the ; is 
        # then adjacent to a quoted text leading to strange behaviour
        
        txt <- paste0(txt, c(rep(";", n_unit - 1), ""))
        txt <- paste0('"', txt, '"')
        txt <- paste(unit_unique, txt, sep = "~")
        txt <- paste(txt, collapse = "~")
        
        out$unit <- paste("Units:", txt, sep = "~")
      }
      
    }    

  }    
        
  out$unit <- parse(text = out$unit)
    

  # extraction
  
  out$extraction <- "Data extraction:"
  if (!is.null(info$extraction)) {
    out$extraction <- paste(out$extraction, format(info$extraction, "%d %B %Y"))
  }

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

  # silence non-standard evaluation warnings
  value <- NULL

  AC <- AC[!is.na(AC)]
  AC <- sort(if (useLogs) log(AC) else AC)
  AC <- data.frame(id = names(AC), value = AC, ok = AC >= ylim[1] & AC <= ylim[2], stringsAsFactors = FALSE)
  within(AC, pos <- convertY(unit(value, "native"), "npc", valueOnly = TRUE))
}



plot_data <- function(
    data, assessment, series, info, type = c("data", "index"), 
    xykey.cex = 1.0, ntick.x = 4, ntick.y = 3, ...) {

  # silence non-standard evaluation warnings
  .data <- year <- censoring <- NULL

  type <- match.arg(type) 

  is.pred <- "pred" %in% names(assessment)
  #   if (is.pred & info$determinand %in% c("VDS", "IMPS", "INTS")) 
  #     assessment$pred <- swap.names(assessment$pred, c("lower", "upper"), 
  #                                   c("ci.lower", "ci.upper"))
  
  is.AC <- !all(is.na(assessment$AC))

  # make data types compatible - i.e. raw data or assessment indices
 
  useLogs <- series$distribution %in% "lognormal"

  data <- switch(type, 
    data = {
      out <- subset(data, select = c(year, concentration, censoring))
      if (useLogs) out <- within(out, concentration <- log(concentration))
      out
    },
    index = {
      out <- assessment$annualIndex
      names(out)[2] <- "concentration"
      if (series$determinand %in% c("VDS", "IMPS", "INTS")) {
        out$censoring <- rep("", nrow(out))
      }
      out
    }
  )    

  
  if (series$distribution == "beta" && is.pred) {
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
  
  if (series$determinand %in% c("VDS", "IMPS", "INTS")) {
    wk <- info$imposex[
      info$imposex$species %in% series$species & info$imposex$determinand %in% series$determinand, 
      "max_value"]
    ylim <- extendrange(c(0, wk), f = 0.04)
  }            
  else {
    ylim <- c(do.call("plot.data.ylim", args.list))
  }

  xlim <- plot.data.xlim(data$year, info$recent_years)

  plot.formula <- data$concentration ~ data$year

  
  # ensures ylabels are formatted correctly and fit into viewport

  ykey <- format(plot.scales(ylim, n = ntick.y, logData = useLogs)) 
  key.ylab.padding <- max(nchar(ykey))


  # sets up viewport so that assessment concentrations and ylabel fit correctly
  
  AC.width <- unit(0, "npc")
  if (is.AC) {

    AC <- plot.AC(assessment$AC, ylim, useLogs)  ## This throws a newPage

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
      plot.panel(x, y, data$censoring, type, AC = assessment$AC, 
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
        series$group %in% "Effects"               ~ "",
        series$determinand %in% "VDS"             ~ "stage",
        series$determinand %in% c("IMPS", "INTS") ~ "",
        TRUE                                    ~ "concentration"
      )
      ylabel <- paste(
        ctsm_get_info(info$determinand, series$determinand, "common_name"), 
        extra
      )
      grid.text(
        ylabel, 0, unit(1, "npc") + unit(1.5, "char"), just = c("left", "bottom"), gp = gpar(cex = xykey.cex))
      upViewport()
    })
  
  newPage <- !is.AC
  plot.setup(newPage)
  pushViewport(viewport(layout.pos.row = 1))
  pushViewport(wk.viewport)
  print(data.plot, newpage = FALSE)
  upViewport()
  upViewport()
  plot_key(series, info, type)
}



plot.setup <- function(newPage) {

  if (newPage) grid.newpage()

  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 4), c("null", "lines")))))
}






#' Text representation of units for plots and reports
#' 
#' Utilty function for adding units to plots in e.g. `plot_assessment` or 
#' html reports in e.g. `report_assessment`. Standard (recognised) units are 
#' prettified, with non-standard units returned unchanged.
#'
#' @param unit Character vector giving the units
#' @param basis Character vector giving the bases. Must be one of "D", "L", 
#' "W" or NA_character_
#' @param normaliser Character vector giving the names of the normaliser (if 
#' there are any)
#' @param normaliser_value Numeric vector giving the values of the normaliser 
#' (if there are any)
#' @param normaliser_unit Character vector giving the units of the normaliser 
#' (if there are any)
#' @param html Logical with TRUE returning an html representation for use in 
#' markdown and FALSE (default) returning an text expression for use in lattice 
#' (grid) graphics
#'
#' @return 
#' A character vector that can be used in markdown (`html = TRUE`) or in 
#' grid graphics ('html = FALSE`). In the latter, need to turn the strings into
#' an expression by `parse(text = result)`.  
#'
#' @examples
#' label_unit("mg/kg", "D", "AL", 5, "%")
#' label_unit("ug/kg", "L", "LIPIDWT%", 5, "%", html = TRUE)
#' 
label_unit <- function(
    unit, 
    basis = NA_character_, 
    normaliser = NA_character_, 
    normaliser_value = NA_real_, 
    normaliser_unit = NA_character_,
    html = FALSE) {

  # when html = FALSE:
  #  incorporate spaces by using ~ or by enclosing in quotes
  #  need to ensure all special characters are enclosed in quotes
  
  if (!is.character(unit)) {
    stop("argument 'unit' must be a character")
  }  
  
  if (!(is.character(basis))) {
    stop("argument 'basis' must be a character")
  }  
  
  if (!all(basis %in% c("W", "D", "L", NA_character_))) {
    stop("unrecognised basis")
  }
  

  # make unit pretty
  
  out <- pretty_unit(unit, html)
  

  # append basis

  if(!all(is.na(basis))) {

    b_txt <- dplyr::case_match(
      basis, 
      "D" ~ "dry weight", 
      "W" ~ "wet weight", 
      "L" ~ "lipid weight", 
      .default = NA_character_
    )
    
    if (html) {
      out_txt <- paste(out, b_txt)
    } else {
      b_txt <- gsub(" ", "~", b_txt, fixed = TRUE)
      out_txt <- paste(out, b_txt, sep = "~")
    }
    
    out <- dplyr::if_else(is.na(basis), out, out_txt)
    
  }
  

  # append normaliser 
  
  if(!all(is.na(normaliser))) {

    n_txt <- label_normaliser(
      normaliser, 
      normaliser_value, 
      normaliser_unit, 
      html = html
    ) 
      
    if (html) {
      out_txt <- paste(out, "normalised to", n_txt)
    } else {
      out_txt <- paste(out, "normalised~to", n_txt, sep = "~") 
    }
    
    out <- dplyr::if_else(is.na(normaliser), out, out_txt)
  }
  
  out
}


#' Text representation of normaliser information for plots and reports
#'
#' Utilty function for providing inforamtion about normalistation to plots in 
#' e.g. `plot_assessment` (where it is used with `label_unit`) or 
#' html reports in e.g. `report_assessment`. 
#'
#' @param normaliser Character vector giving the names of the normaliser (ICES 
#' code rather than common name)
#' @param value Numeric vector giving the value of the normaliser 
#' @param unit Character vector giving the unit of the normaliser 
#' @param html Logical with TRUE returning an html representation for use in 
#' markdown and FALSE (default) returning an text expression for use in lattice 
#' (grid) graphics
#'
#' @return
#' A character string that can be used in markdown (`html = TRUE`) or used in 
#' grid graphics ('html = FALSE`). In the latter, need to turn the string into
#' an expression by `parse(text = result)`.  
#' 
#' @examples 
#' label_normaliser("AL", 2.5, "%")
#' label_normaliser("AL", 2.5, "%", html = TRUE)
#' 
label_normaliser <- function(normaliser, value, unit, html = FALSE) {

  # ideally want to get the translations from ICES code to text from 
  # info$determinand, but too many structures called info flying around -  
  # deal with later

  # when html = FALSE:
  # incorporate spaces by using ~ or by enclosing in quotes
  # need to ensure all special characters are enclosed in quotes

  sep <- if (html) " " else "~"
    
  normaliser <- dplyr::case_match(
    normaliser, 
    "AL" ~ "aluminium",
    "LI" ~ "lithium",
    "CORG" ~ "organic carbon",
    "LIPIDWT%" ~ "lipid",
    .default = normaliser
  )
  
  # put normaliser in quotes (when html = FALSE) in case any special characters 
  # are carried through
  
  if (!html) {
    normaliser <- dplyr::if_else(
      is.na(normaliser), 
      NA_character_, 
      paste0('"', normaliser, '"')
    )
  }

  
  # no space between value and unit (used for % and %%)

  out <- dplyr::case_when(
    unit %in% c("%", "%%") ~ 
      if (html) {
        paste0(value, unit)
      } else {
        paste0('"', value, unit, '"')
      },
    !is.na(unit) ~ paste(value, pretty_unit(unit), sep = sep), 
    .default = NA_character_
  )
    
  out <- dplyr::if_else(
    is.na(normaliser), 
    NA_character_, 
    paste(out, normaliser, sep = sep)
  )
  
  out
}



#' Pretty representation of units for plots and reports
#' 
#' Utilty function for adding units to plots in e.g. `plot_assessment` or 
#' html reports in e.g. `report_assessment`. Standard (recognised) units are 
#' prettified, with non-standard units returned unchanged.
#'
#' @param unit 
#' @param html A logical with TRUE returning an html representation for use in 
#' markdown and FALSE (default) returning an text expression for use in lattice 
#' (grid) graphics
#'
#' @return 
#' A vector of character strings that can be used in markdown (`html = TRUE`) or 
#' in grid graphics ('html = FALSE`). In the latter, need to turn the strings 
#' into expressions by `parse(text = result)`.  
#'
#' @examples
#' pretty_unit("ug/kg")
#' pretty_unit("ug/kg", html = TRUE)
#' 
pretty_unit <- function(unit, html = FALSE) {
  
  if (!is.character(unit)) {
    stop("argument 'unit' must be a character")
  }
    
  # when html = FALSE:
  # incorporate spaces by using ~ or by enclosing in quotes
  # need to ensure all special characters are enclosed in quotes
  
  if (html) {
    unit <- dplyr::case_match(
      unit, 
      "mg/kg" ~ "mg kg<sup>-1</sup>",
      "mg/g"  ~ "mg g<sup>-1</sup>",
      "ug/kg" ~ "&mu;g kg<sup>-1</sup>", 
      "ug/g"  ~ "&mu;g g<sup>-1</sup>", 
      "ng/kg" ~ "ng kg<sup>-1</sup>",
      "ng/g"  ~ "ng g<sup>-1</sup>",
      "ug/l"  ~ "&mu;g l<sup>-1</sup>",
      "ug/ml" ~ "&mu;g ml<sup>-1</sup>", 
      "ng/l"  ~ "ng l<sup>-1</sup>",
      "ng/ml" ~ "ng ml<sup>-1</sup>",
      "mins"  ~ "min",
      "stg"   ~ "stage",
      "j/h/g" ~ "J h <sup>-1</sup> g <sup>-1</sup>",
      "pmol/min/mg protein" ~ "pmol min<sup>-1</sup> mg protein<sup>-1</sup>",
      "nmol/min/mg protein" ~ "nmol min<sup>-1</sup> mg protein<sup>-1</sup>",
      "ng/min/mg protein"   ~ "ng min<sup>-1</sup> mg protein<sup>-1</sup>", 
      .default = unit
    )
  } else {
    unit <- dplyr::case_match(
      unit, 
      "mg/kg" ~ 'mg ~ kg^{-1}',
      "mg/g"  ~ 'mg ~ g^{-1}',
      "ug/kg" ~ 'paste(mu, "g") ~ kg^{-1}', 
      "ug/g"  ~ 'paste(mu, "g") ~ g^{-1}', 
      "ng/kg" ~ 'ng ~ kg^{-1}',
      "ng/g"  ~ 'ng ~ g^{-1}',
      "ug/l"  ~ 'paste(mu, "g") ~ l^{-1}', 
      "ug/ml" ~ 'paste(mu, "g") ~ ml^{-1}', 
      "ng/l"  ~ 'ng ~ l^{-1}', 
      "ng/ml" ~ 'ng ~ ml^{-1}',
      "mins"  ~ 'min', 
      "stg"   ~ 'stage',
      "j/h/g" ~ 'J ~ h^{-1} ~ g^{-1}',
      "pmol/min/mg protein" ~ 'pmol ~ min^{-1} ~ "mg protein"^{-1}',
      "nmol/min/mg protein" ~ 'nmol ~ min^{-1} ~ "mg protein"^{-1}',
      "ng/min/mg protein"   ~ 'ng ~ min^{-1} ~ "mg protein"^{-1}', 
      .default = paste0('"', unit, '"')
    )
  }
  
  unit
}



plot.panel <- function(
  x, y, censoring, type, AC = NA, pred = NULL, ylim, layout.row, useLogs = TRUE, 
  indiCL) {

  # type
  # data is standard (single panel) data plot with (optional) assessment
  # index is standard (single panel) index plot with (optional) assessment 
  # ratio_mp is multipanel ratio plot (under development)
  # index_mp is multipanel index plot of related compounds with (optional) 
  #  assessment
  
  type <- match.arg(type, c("data", "index", "ratio_mp", "index_mp"))
  
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
    index = 2,
    ratio_mp = switch(layout.row, 2.0, 1.4, 0.9, 0.7, 0.6), 
    index_mp = switch(layout.row, 2.0, 1.4, 0.9, 0.7, 0.6) 
  )
  
  wk.pch <- switch(
    type, 
    data = "+", 
    index = 16, 
    ratio_mp = "+", 
    index_mp = 16
  )
  
  wk.cex.censoring = if (wk.pch == "+") wk.cex * 0.8 else wk.cex
  
  # recognise following symbols for censoring: ">" and "?" can arise in ratio plots

  censoring <- as.character(censoring)
  
  is_censoring <- censoring %in% c("<", "D", "Q", ">", "?")

  censoring[is_censoring] <- tolower(censoring[is_censoring])
  
  if (any(is_censoring)) {
    lattice::lpoints(
      x[is_censoring], 
      y[is_censoring], 
      pch = censoring[is_censoring], 
      cex = wk.cex.censoring, 
      col = "black"
    )
  }
    
  if (any(!is_censoring)) {
    lattice::lpoints(
      x[!is_censoring], 
      y[!is_censoring], 
      pch = wk.pch, 
      cex = wk.cex, 
      col = "black"
    )
  }

  if (!missing(type) && type == "index" && "lower" %in% names(indiCL)) {
    lattice::lsegments(
      indiCL$year, 
      indiCL$lower, 
      indiCL$year, 
      indiCL$upper, 
      lwd = 2, 
      col = "black"
    )
  }
}


plot_auxiliary <- function(
    data, assessment, series, info, auxiliary,  
    xykey.cex = 1.0, ntick.x = 3, ntick.y = 3, newPage = TRUE, ...) {

  # silence non-standard evaluation warnings
  .data <- info.imposex <- NULL

  # auxiliary specifies the choice of 'auxiliary' variables to plot: 
  # default:
  #   sediment = value, concentration, AL, CORG
  #   biota = concentration, LNMEA, DRYWT%, LIPIDWT%
  #   water = not specified yet
  # otherwise must contain four relevant variables
  
  useLogs <- series$distribution %in% "lognormal"
  
  data <- dplyr::mutate(
    data, 
    concentration = if (useLogs) log(.data$concentration) else .data$concentration,
    concentration.censoring = .data$censoring
  )
  
  if ("concOriginal" %in% names(data)) {
    data <- dplyr::mutate(
      data, 
      value = if (useLogs) log(.data$concOriginal) else .data$concOriginal,
      value.censoring = .data$censoringOriginal
    )
  }


  auxiliary.censoring <- paste(auxiliary, "censoring", sep = ".")
  
  # not all auxiliary variables have censorings, so create dummy columns
  
  ok <- auxiliary.censoring %in% names(data)
  if (!all(ok)) {
    data[auxiliary.censoring[!ok]] <- lapply(
      data[auxiliary[!ok]], 
      function(x) ifelse(is.na(x), NA, "")
    )
  }
  data <- data[c("year", auxiliary, auxiliary.censoring)]


  data <- reshape(
    data, 
    varying = list(auxiliary, auxiliary.censoring), 
    v.names = c("value", "censoring"), 
    direction = "long", 
    timevar = "type", 
    times = auxiliary
  )
  
  data <- within(data, type <- ordered(type, levels = auxiliary))

  xlim <- range(data$year, info$recent_years)

  is.data <- unlist(with(data, tapply(value, type, function(i) !all(is.na(i)))))
  data <- within(data, value[type %in% names(is.data[!is.data])] <- 0)

  # this is what R will use in xyplot
  ylim <- with(data, tapply(value, type, extendrange, f = 0.07))        
  if (info$compartment == "sediment") {
    ylim$value <- with(
      subset(data, type %in% c("value", "concentration")), 
      extendrange(value, f = 0.07)
    )
    ylim$concentration <- ylim$value
  }
  
  if (series$determinand %in% c("VDS", "IMPS", "INTS")) {
    wk <- dplyr::filter(
      info$imposex, 
      .data$species == series$species,
      .data$determinand == series$determinand
    )
    imposex_max_value <- wk$max_value
    ylim$concentration <- extendrange(c(0, imposex_max_value), f = 0.07)
  }            
  
  
  ykey <- sapply(levels(data$type), function(i) {
    if (is.data[i]) {
      wk <- plot.scales(
        ylim[[i]], 
        n = ntick.y, 
        logData = (useLogs & i %in% c("value", "concentration"))
      )
      max(nchar(format(wk)))
    }
    else 0
  })
  key.ylab.padding <- max(ykey[c(1, 3)])
  
  # but this is what must get passed in!
  ylim <- with(data, tapply(value, type, range, na.rm = TRUE))          
  if (info$compartment == "sediment") {
    ylim$value <- with(
      subset(data, type %in% c("value", "concentration")), 
      range(value, na.rm = TRUE)
    )
    ylim$concentration <- ylim$value
  }

  if (series$determinand %in% c("VDS", "IMPS", "INTS")) {
    ylim$concentration <- c(0, imposex_max_value)
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
          right.padding = 0, key.right = 0, axis.key.padding = 0, 
          axis.right = 0, strip.left = 0, key.left = 0, axis.panel = 0
        ),
        layout.heights = list(
          axis.bottom = 0, bottom.padding = 2, axis.xlab.padding = 0, xlab = 0, 
          xlab.key.padding = xykey.cex, key.sub.padding = 0, 
          axis.top = 0, top.padding = 0, main = 0, main.key.padding = 0, 
          key.top = 0, key.axis.padding = 0, axis.panel = 0)
      ), 
      axis = function(side, ...) 
        plot.axis(
          side, ntick.x = ntick.x, ntick.y = ntick.y, xykey.cex = xykey.cex, 
          plot.type = "auxiliary", is.data = is.data, useLogs = useLogs, ...
        ),
      panel = function(x, y, subscripts) {
        type.id <- levels(type)[which.packet()]
        
        if (info$compartment == "sediment" && 
            (!info$normalise || is.na(series$normaliser)) && 
            type.id == "concentration"
        ) {
          grid.text("data not-normalised", 0.5, 0.5, gp = gpar(cex = xykey.cex))
        } else {
          if (any(duplicated(x))) {
            x <- jitter(x, amount = 0.1)
          }
          censoring <- tolower(as.character(censoring[subscripts]))
          censoring <- ifelse(censoring %in% "", "+", censoring)
          lpoints(x, y, pch = censoring, cex = 2.5, col = "black")
        }

        pushViewport(viewport(clip = "off"))
        ylabel <- switch(
          type.id, 
          concentration = switch(
            info$compartment, 
            biota = paste(
              series$determinand, 
              dplyr::case_when(
                series$group %in% "Effects"               ~ "",
                series$determinand %in% "VDS"             ~ "stage",
                series$determinand %in% c("IMPS", "INTS") ~ "",
                TRUE                                      ~ "concentration"
              )
            ),
            sediment = paste(series$determinand, "normalised")
          ),
          value = paste(series$determinand, "non-normalised"),
          LNMEA = {
            family <- ctsm_get_info(info$species, series$species, "species_group")
            unit <- ctsm_get_info(
              info$determinand, 
              type.id, "unit", 
              info$compartment, 
              sep = "_"
            )
            paste0(
              "Mean ", 
              switch(
                family, 
                Fish = "fish ", 
                Bivalve = "shell ", 
                Gastropod = "shell ", 
                ""
              ), 
              "length (", unit, ")", 
              sep = ""
            )
          },
          {
            unit <- ctsm_get_info(
              info$determinand, 
              type.id, 
              "unit", 
              info$compartment, 
              sep = "_"
            )
            paste0(
              ctsm_get_info(info$determinand, type.id, "common_name"), 
              " (", unit, ")"
            )
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

  plot_key(series, info, "auxiliary")
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



plot_multiassessment <- function(data, assessment, series, info, ...) {

  # silence non-standard evaluation warnings
  .data <- NULL

  is.data <- sapply(assessment, function(i) !is.null(i))
  
  is.pred <- sapply(assessment, function(i) !is.null(i) && !is.null(i$pred))
  is.AC <- sapply(assessment, function(i) !is.null(i) && !all(is.na(i$AC)))
  
  series_distribution <- ctsm_get_info(
    info$determinand, 
    series$determinand, 
    "distribution"
  )
  
  useLogs <- series_distribution %in% "lognormal"
  
  names(useLogs) <- series$seriesID
  

  # make data types compatible - i.e. raw data or assessment indices

  data <- sapply(series$seriesID, simplify = FALSE, function(i) {
    if (is.data[i]) {
      out <- assessment[[i]]$annualIndex
      names(out)[2] <- "concentration"
      out
    }
    else data.frame(
      year = info$max_year, 
      concentration = 0, 
      censoring = factor("", levels = c("", "<"))
    )
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

  ylim <- sapply(series$seriesID, simplify = FALSE, function(i) {
    # NB have taken logs above, so this is log concentration
    args.list <- list(data[[i]]$concentration)         
    if (is.pred[[i]]) {
      args.list <- c(
        args.list, 
        assessment[[i]]$pred$ci.lower, 
        assessment[[i]]$pred$ci.upper
      )
    }
    do.call("plot.data.ylim", args.list)
  })

  args.list <- sapply(
    series$seriesID[is.data], 
    simplify = FALSE, 
    function(i) data[[i]]$year
  )
  args.list <- c(args.list, info$recent_years)
  xlim <- do.call("plot.data.xlim", args.list)

  plot.formula <- data$concentration ~ data$year


  # ensures ylabels are formatted correctly and fit into viewport

  ntick.y <- ntick.x <- 3
  ykey <- sapply(series$seriesID, simplify = FALSE, function(i) 
    format(plot.scales(ylim[[i]], n = ntick.y, logData = useLogs[i])))
  key.ylab.padding <- max(sapply(ykey, function(i) max(nchar(i))))
  

  ndet <- length(series$seriesID)
  layout.row <- ceiling(sqrt(ndet))
  layout.col <- ceiling(ndet / layout.row)

  add.xlab = 1:ndet <= layout.col
  names(add.xlab) <- series$seriesID

  xykey.cex <- switch(layout.row, 1.4, 1.1, 0.9, 0.7, 0.6)


  # sets up viewport so that assessment concentrations and ylabel fit correctly
  
  xAC <- 0.5

  AC.width <- sapply(series$seriesID[is.data], simplify = FALSE, function(i) {
    
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


  data.plot <- sapply(series$seriesID, simplify = FALSE, function(i) {
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
            x, y, data[[i]]$censoring, 
            type = "index_mp",
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
        grid.text(series$plotNames$assessment[i], 0, unit(1, "npc") + unit(1, "char"), 
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
   
   plot_key(series, info, "index_mp")
}



plot_multidata <- function(data, series, info,  ...) {

  data <- subset(data, !is.na(concentration))

  series_distribution <- ctsm_get_info(info$determinand, data$determinand, "distribution")
  
  useLogs <- series_distribution %in% "lognormal"
  
  data <- within(data, concentration[useLogs] <- log(concentration[useLogs]))
  data <- data[c("year", "sample", "seriesID", "censoring", "concentration")]


  data <- reshape(data, direction = "wide", idvar = c("sample", "year"), timevar = "seriesID")

  # add in extra rows for recent years to ensure there is a sensible range of years and to cover 
  # situation where a determinand has only a single value and a range can't be calculated - could be done
  # much more elegantly

#  if (max(data$year) < max(info$recent_years) | min(data$year) > min(info$recent_years)) 
#  {
#    last.id <- nrow(data) + c(1, 2)
#    data[last.id,] <- NA
#    data[last.id, "year"] <- range(info$recent_years)
#  }

  # varnames doesn't appear to work in splom, so make sure we give the columns the names
  # we want printed
  
  plot.data <- data[c("year", paste("concentration", series$seriesID, sep = "."))]
  names(plot.data) <- varnames <- c("year", series$plotNames$data)
  
  colID <- c("year", series$seriesID)


  pscales <- lapply(names(plot.data), function(i) {
    if (i == "year") 
      limits <- plot.data.xlim(plot.data[i], info$recent_years) 
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
      layout.heights = list(
        bottom.padding = 0, axis.bottom = 0, axis.xlab.padding = 0, xlab = 0
      )
    ),
    panel = function(x, y, i, j) {
      censoring <- if (i > 1) {
        data[paste("censoring", colID[i], sep = ".")]
      } else {
        rep("", length(x))
      }
      lpoints(
        x, y, col = "black", pch = ifelse(censoring == "", "+", "<"), cex = 1.1
      )
    }  
  )

  plot.setup(newPage = TRUE)
  pushViewport(viewport(layout.pos.row = 1))
  print(data.plot, newpage = FALSE)
  upViewport()
  plot_key(series, info, "data_splom")
}


ctsm.panel.pairs <- function (z, panel = lattice::lattice.getOption("panel.splom"), lower.panel = panel, 
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



plot_key <- function(series, info, type) {

  type <- match.arg(
    type, 
    c("data", "index", "auxiliary", "data_splom", "index_mp", "ratio_mp")
  )
  
  key <- make_key(series, info, type)
  
  pushViewport(viewport(layout.pos.row = 2))
  
  grid.text(key$media, x = unit(1, "char"), y = unit(4, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  grid.text(key$station, x = unit(1, "char"), y = unit(3, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  grid.text(key$unit, x = unit(1, "char"), y = unit(2, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  grid.text(key$extraction, x = unit(1, "char"), y = unit(1, "lines"), gp = gpar(cex = 0.8), 
            just = c("left", "centre"))
  
  invisible()
}



plot_ratio_data <- function(data, numerator, denominator, type = c("logistic", "log")) {

  # silence non-standard evaluation warnings
  .data <- NULL

  type <- match.arg(type)
  
  id <- c(
    "year", 
    paste(c("concentration", "censoring"), numerator, sep = "."),
    paste(c("concentration", "censoring"), denominator, sep = ".")
  )

  
  # calculate ratios
    
  data <- data[id]
  names(data) <- c("year", "n_conc", "n_censoring", "d_conc", "d_censoring")
  
  data <- dplyr::mutate(
    data, 
    ratio = switch(
      type,
      logistic = .data$n_conc / (.data$n_conc + .data$d_conc),
      log = .data$n_conc / .data$d_conc
    ), 
    censoring = dplyr::case_when(
      .data$n_censoring %in% "" & .data$d_censoring %in% ""               ~ "+",
      .data$n_censoring %in% "" & .data$d_censoring %in% c("<", "D", "Q") ~ ">",
      .data$n_censoring %in% c("<", "D", "Q") & .data$d_censoring %in% "" ~ "<",
      !is.na(.data$n_censoring) & !is.na(.data$d_censoring)               ~ "?",
      TRUE                                                        ~ NA_character_
    )
  )
  
  data <- data[c("year", "ratio", "censoring")]
  
  data
}


plot_ratio_pred <- function(
  data, 
  type = c("logistic", "log"), 
  control = list(nyear = 5, prop_censoring = 0.1)
) {

  # silence non-standard evaluation warnings
  se <- NULL

  type <- match.arg(type)
  
  data <- na.omit(data)
  
  # only fit smoother if number of years >= 5 and number of 'less-thans' < 10%
  
  nyear <- length(unique(data$year))
  
  prop_censoring <- sum(!data$censoring %in% "+") / length(data$censoring) 
  
  if (nyear < control$nyear || prop_censoring >= control$prop_censoring ) {
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



plot_ratio <- function(data, series, info, ...) {
  
  # get working data 
  # sediment - use non-normalised concentrations
  # biota - could use values before conversion to target bases, but would need to 
  #   ensure comparability - maybe later
  # don't need to log transform, because looking at ratios - maybe in plots

  if (info$compartment == "sediment") {
    data$concentration <- data$concOriginal
    data$censoring <- data$censoringOriginal
  }
  
  # set up ratios
  
  det_id <- switch(
    series$group, 
    Metals = c("SE", "HG"),
    PAH_parent = c("ANT", "PA", "FLU", "PYR", "ICDP", "BGHIP", "BAA", "CHR"), 
    PBDEs = c("BDE47", "BD153"), 
    Organofluorines = c("PFNA", "PFOA", "PFUNDA", "PFDA", "PFTRDA", "PFDOA"),
    Organochlorines = c("DDEPP", "DDTPP", "DDTOP")
  )
  
  ratio_id <- switch(
    series$group, 
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
    series$group,
    Metals = "SE",
    PAH_parent = c("ANT", "FLU", "ICDP", "BAA"), 
    PBDEs = "BDE47",
    Organofluorines = c("PFNA", "PFUNDA", "PFTRDA"),
    Organochlorines = c("DDEPP", "DDTOP")
  )
  
  denominator_id <- switch(
    series$group,
    Metals = "HG",
    PAH_parent = c("PA", "PYR", "BGHIP", "CHR"),
    PBDEs = "BD153",
    Organofluorines = c("PFOA", "PFDA", "PFDOA"),
    Organochlorines = c("DDTPP", "DDTPP")
  )
  
  ratio_type <- switch(
    series$group,
    PAH_parent = "logistic",
    "log"
  )
  
  use_logs <- ratio_type %in% "log"
  
  ref_lines <- switch(
    series$group, 
    Metals = list(78.96 / 200.59),
    PAH_parent = list(0.1, c(0.4, 0.5), c(0.2, 0.5), c(0.2, 0.35)),
    PBDEs = list(NA), 
    Organofluorines = list(NA, NA, NA),
    Organochlorines = list(1, c(0.3, 1))
  )
  
  ref_txt <- switch(
    series$group, 
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
  
  data <- data[c("year", "sample", "determinand", "censoring", "concentration")]
  
  data <- reshape(
    data, 
    direction = "wide", 
    idvar = c("sample", "year"), 
    timevar = "determinand"
  )

  
  # add in dummy columns to deal with variables that are not reported
  
  if (length(missing_det) > 0) {
    new_id <- paste("censoring", missing_det, sep = ".")
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
    FUN = plot_ratio_data,
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
  
  
  pred <- lapply(data, plot_ratio_pred, type = ratio_type)
  
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
  args.list <- c(args.list, info$recent_years)
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
            x, y, data[[i]]$censoring,
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
  
  
  # hack of plot_key to allow for dimensionless units
  
  key <- make_key(series, info, type = "ratio_mp")
  
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
