#' @export
ctsm_uncrt_workup <- function(clean_data) {

  library(dplyr)
  
  # turn 'clean' data into uncertainty data
  
  # read in data
  
  data <- clean_data$data
  stations <- clean_data$stations
  compartment <- clean_data$info$compartment

  rm(clean_data)
  
  
  # link to country
  
  data$country <- stations[as.character(data$station), "country"]
  

  # get alabo and remove missing alabo
  
  data <- within(data, {
    alabo <- sapply(strsplit(as.character(qaID), "_"), "[", 3)
    alabo[alabo %in% "NA"] <- NA
    alabo <- factor(alabo)
  })
  
  data <- data[!is.na(data$alabo), ]
  

  # retain useful variables (with special case for sediment)
  # NB don't drop data with missing uncertainty here - might lose auxiliary 
  #   (CORG and AL) information

  id_aux <- c(
    "", ".uncertainty", ".censoring", ".limit_detection", ".limit_quantification"
  )
  
  id <- intersect(
    c("country", "alabo", "year", "sample", "group", "determinand", 
      "concentration", "uncertainty", 
      "censoring", "limit_detection", "limit_quantification", 
      paste0("AL", id_aux), 
      paste0("LI", id_aux),
      paste0("CORG", id_aux), 
      paste0("LOIGN", id_aux)),
    names(data)
  )
  data <- data[id]
  

  # sort out AL and CORG etc for sediment
  
  if (compartment == "sediment") {
    
    id <- c("country", "alabo", "year", "group", "sample", "determinand")

    out_names <- c(
      "concentration", "uncertainty", "censoring", "limit_detection", 
      "limit_quantification"
    )
    
    AL_names <- paste0("AL", id_aux)
    AL <- data[c(id, AL_names)]
    AL <- mutate(AL, determinand = "AL", group = "auxiliary")
    AL <- unique(AL)
    pos <- match(AL_names, names(AL))
    names(AL)[pos] <- out_names

    LI_names <- paste0("LI", id_aux)
    LI <- data[c(id, LI_names)]
    LI <- mutate(LI, determinand = "LI", group = "auxiliary")
    LI <- unique(LI)
    pos <- match(LI_names, names(LI))
    names(LI)[pos] <- out_names

    CORG_names <- paste0("CORG", id_aux)
    CORG <- data[c(id, CORG_names)]
    CORG <- mutate(CORG, determinand = "CORG", group = "auxiliary")
    CORG <- unique(CORG)
    pos <- match(CORG_names, names(CORG))
    names(CORG)[pos] <- out_names
    
    pos <- match(c(AL_names, LI_names, CORG_names), names(data))
    data <- data[-pos]

    is_LOIGN <- "LOIGN" %in% names(data)
    
    if (is_LOIGN) {
      LOIGN_names <- paste0("LOIGN", id_aux)
      LOIGN <- data[c(id, LOIGN_names)]
      LOIGN <- mutate(LOIGN, determinand = "LOIGN", group = "auxiliary")
      LOIGN <- unique(LOIGN)
      pos <- match(LOIGN_names, names(LOIGN))
      names(LOIGN)[pos] <- out_names
      
      pos <- match(LOIGN_names, names(data))
      data <- data[-pos]
    }    
      
    data <- rbind(data, AL, LI, CORG)
    
    if (is_LOIGN) {
      data <- rbind(data, LOIGN)
    }
  }

  
  # remove missing uncertainties
  
  data <- data[!is.na(data$uncertainty), ]


  # restrict to 'log-normally' distributed responses
  
  ok <- with(data, {
    dist <- ctsm_get_info(
      "determinand", determinand, "distribution", na_action = "output_ok"
    ) 
    dist %in% "lognormal" | determinand %in% c("CORG", "LOIGN")
  })
  
  data <- data[ok, ]


  # order groups and determinands within group
  
  det_list <- determinands[[stringr::str_to_title(compartment)]]

  data <- within(data, {
    group <- factor(as.character(group), levels = c(names(det_list), "auxiliary"))
    determinand <- factor(
      as.character(determinand), 
      levels = c(unlist(det_list), "AL", "LI", "CORG", "LOIGN"))
  })

  
  # calculate relative uncertainty
  
  data <- within(data, relative_u <- 100 * uncertainty / concentration)

  data <- droplevels(data)

  list(compartment = compartment, data = data)
}

#' @export
ctsm_uncrt_estimate <- function(data) {
  
  # initialise output with total number of values by determinand
  
  options(dplyr.summarise.inform = FALSE)
  on.exit(options(dplyr.summarise.inform = NULL))

  out <- data %>% 
    group_by(.data$determinand) %>% 
    summarise(n_values = n())

  
  # remove duplicate combinations of concentration and uncertainty (and associated censoring variables)
  # big problems when there are lots of less thans
  
  id <- c(
    "determinand", "concentration", "uncertainty", "censoring", 
    "limit_detection", "limit_quantification", "alabo"
  )
  
  dup <- duplicated(data[id])
  data <- data[!dup, ]  
  
  
  # get number of 'unique values
  
  out_unique <- data %>% 
    group_by(.data$determinand) %>% 
    summarise(n_unique = n(), n_alabo = n_distinct(.data$alabo))

  out <- left_join(out, out_unique, by = "determinand")
  

  # remove data that is way off
  
  data <- data %>% 
    filter(between(relative_u, 1, 100)) 
  
  
  # relative error
  # median relative_u for values above the detection level by alabo
  
  out_relative <- data %>% 
    filter(.data$censoring == "") %>% 
    group_by(.data$determinand, .data$alabo) %>% 
    summarise(sd_variable = median(.data$relative_u) / 100) 
  
  # now the median value across alabos
  
  out_relative <- out_relative %>% 
    group_by(.data$determinand) %>% 
    summarise(sd_variable = median(sd_variable))
  
  out <- left_join(out, out_relative, by = "determinand")
  
  
  # constant error
  # median limit_detection for values with censoring == D, Q or "" by alabo
  # don't use "<" because we can't trust any of the limit values
  
  out_constant <- data %>% 
    filter(.data$censoring %in% c("D", "Q", "")) %>% 
    drop_na(.data$limit_detection) %>% 
    group_by(.data$determinand, .data$alabo) %>% 
    summarise(sd_constant = median(.data$limit_detection) / 3) 
  
  # now the median value across alabos
  
  out_constant <- out_constant %>% 
    group_by(.data$determinand) %>% 
    summarise(sd_constant = median(sd_constant))
  
  out <- left_join(out, out_constant, by = "determinand")
  
  # tidy up
  
  out <- out %>% 
    as.data.frame() %>% 
    column_to_rownames("determinand") %>% 
    round(6)
    
  out
}

#' @export
ctsm_uncrt_plot_estimates <- function(uncrt_obj, old_estimates, group_id) {

  library(lattice)

  id <- with(uncrt_obj$data, group %in% group_id)  
  data <- uncrt_obj$data[id, ]
  
  data <- data[with(data, order(determinand, concentration)), ]
  
  ok <- with(data, relative_u >= 1 & relative_u <= 100)
  data <- data[ok, ]

  new <- uncrt_obj$estimates[c("sd_constant", "sd_variable")]
  names(new) <- c("sdC", "sdV")
  
  var_id <- paste(uncrt_obj$compartment, c("sd_constant", "sd_variable"), sep= ".")
  old <- old_estimates[var_id]
  names(old) <- c("sdC", "sdV")
  
  xyplot(
    relative_u ~ concentration | determinand, data = data, 
    aspect = 1,
    scales = list(alternating = FALSE, x = list(log = TRUE, relation = "free", equispaced = FALSE)), 
    as.table = TRUE,
    panel = function(x, y, subscripts) {
      data <- data[subscripts, ]
      lod <- with(data, limit_detection[censoring %in% c("D", "Q", "")])
      lod <- na.omit(lod)
      panel.rug(x = log10(lod), lwd = 2, col = "red")

      lpoints(x, y, col = grey(0.7))

      det_id <- unique(as.character(data$determinand))
      
      old <- old[det_id, ]
      new <- new[det_id, ]
      
      conc <- 10^x
      
      fit <- 100 * sqrt((old$sdC / conc) ^ 2 + (old$sdV ^ 2))
      llines(x, fit, lwd = 2, col = "blue")
      
      fit <- 100 * sqrt((new$sdC / conc) ^ 2 + (new$sdV ^ 2))
      llines(x, fit, lwd = 2, col = "red")
    }
  )
}


#' Generates a relative uncertainty plot
#' 
#' @param uncrt_obj an uncertainty object
#' @param det_id a list of determinand ids
#' @export
ctsm_uncrt_plot_data <- function(uncrt_obj, det_id) {
  
  library(lattice)

  id <- with(uncrt_obj$data, determinand %in% det_id)  
  data <- uncrt_obj$data[id, ]

  xyplot(
    relative_u ~ concentration | country, data = data, groups = alabo, 
    ylab = "relative uncertainty (%)",
    scales = list(alternating = FALSE, equispaced.log = FALSE, log = TRUE), 
    as.table = TRUE, aspect = 0.7,
    panel = function(x, y, subscripts, ...) {
      data <- data[subscripts, ]
      lod <- with(data, limit_detection[censoring %in% c("D", "Q", "")])
      lod <- na.omit(lod)
      panel.rug(x = log10(lod), lwd = 2)
      panel.superpose(x, y, subscripts, ..., panel.groups = "panel.xyplot")
    }
  )
}
