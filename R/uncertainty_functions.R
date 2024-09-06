#' @export
ctsm_uncrt_workup <- function(harsat_obj) {

  # silence non-standard evaluation warnings
  .data <- NULL
  
  # turn 'clean' data into uncertainty data
  
  # read in data
  
  data <- harsat_obj$data
  stations <- harsat_obj$stations
  info <- harsat_obj$info

  rm(harsat_obj)
  
  
  # link to country
  
  data <- dplyr::left_join(
    data, 
    stations[c("station_code", "country")], 
    by = "station_code"
  )
  

  # remove data with no analytical laboratory information
  
  data <- data[!is.na(data$alabo), ]
  

  # retain useful variables (with special case for sediment)
  # NB don't drop data with missing uncertainty here - might lose auxiliary 
  #   (CORG and AL) information

  id_aux <- c(
    "", ".uncertainty", ".censoring", ".limit_detection", ".limit_quantification"
  )

  
  id <- c(
    "country", "alabo", "year", "sample", "group", "determinand", 
    "concentration", "uncertainty", 
    "censoring", "limit_detection", "limit_quantification", 
    paste0("AL", id_aux), 
    paste0("LI", id_aux),
    paste0("CORG", id_aux), 
    paste0("LOIGN", id_aux)
  )
  
  data <- dplyr::select(data, dplyr::any_of(id))
  
  

  # sort out AL and CORG etc for sediment
  
  if (info$compartment == "sediment") {
    
    id <- c("country", "alabo", "year", "group", "sample", "determinand")

    out_names <- c(
      "concentration", "uncertainty", "censoring", "limit_detection", 
      "limit_quantification"
    )
    
    AL_names <- paste0("AL", id_aux)
    AL <- data[c(id, AL_names)]
    AL <- dplyr::mutate(AL, determinand = "AL", group = "auxiliary")
    AL <- unique(AL)
    pos <- match(AL_names, names(AL))
    names(AL)[pos] <- out_names

    LI_names <- paste0("LI", id_aux)
    LI <- data[c(id, LI_names)]
    LI <- dplyr::mutate(LI, determinand = "LI", group = "auxiliary")
    LI <- unique(LI)
    pos <- match(LI_names, names(LI))
    names(LI)[pos] <- out_names

    CORG_names <- paste0("CORG", id_aux)
    CORG <- data[c(id, CORG_names)]
    CORG <- dplyr::mutate(CORG, determinand = "CORG", group = "auxiliary")
    CORG <- unique(CORG)
    pos <- match(CORG_names, names(CORG))
    names(CORG)[pos] <- out_names
    
    pos <- match(c(AL_names, LI_names, CORG_names), names(data))
    data <- data[-pos]

    is_LOIGN <- "LOIGN" %in% names(data)
    
    if (is_LOIGN) {
      LOIGN_names <- paste0("LOIGN", id_aux)
      LOIGN <- data[c(id, LOIGN_names)]
      LOIGN <- dplyr::mutate(LOIGN, determinand = "LOIGN", group = "auxiliary")
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
  # keep explicit mention of CORG and LOIGN just in case
  
  data <- dplyr::mutate(
    data, 
    distribution = ctsm_get_info(
      info$determinand, 
      .data$determinand, 
      "distribution", 
      na_action = "output_ok"
    ) 
  )  

  data <- dplyr::filter(
    data, 
    .data$distribution %in% "lognormal" | .data$determinand %in% c("CORG", "LOIGN")
  )
  

  # order groups and determinands within group
  
  # det_list <- determinands[[stringr::str_to_title(compartment)]]
  # 
  # data <- within(data, {
  #   group <- factor(as.character(group), levels = c(names(det_list), "auxiliary"))
  #   determinand <- factor(
  #     as.character(determinand), 
  #     levels = c(unlist(det_list), "AL", "LI", "CORG", "LOIGN"))
  # })

  
  # calculate relative uncertainty
  
  data <- dplyr::mutate(
    data, 
    relative_u = 100 * .data$uncertainty / .data$concentration
  )

  list(compartment = info$compartment, data = data)
}


#' @export
ctsm_uncrt_estimate <- function(data) {
  
  # silence non-standard evaluation warnings
  .data <- NULL

  # initialise output with total number of values by determinand
  
  options(dplyr.summarise.inform = FALSE)
  on.exit(options(dplyr.summarise.inform = NULL))

  out <- data |> 
    dplyr::group_by(.data$determinand) |> 
    dplyr::summarise(n_values = dplyr::n())

  
  # remove duplicate combinations of concentration and uncertainty (and associated censoring variables)
  # big problems when there are lots of less thans
  
  id <- c(
    "determinand", "concentration", "uncertainty", "censoring", 
    "limit_detection", "limit_quantification", "alabo"
  )
  
  dup <- duplicated(data[id])
  data <- data[!dup, ]  
  
  
  # get number of 'unique values
  
  out_unique <- data |> 
    dplyr::group_by(.data$determinand) |> 
    dplyr::summarise(
      n_unique = dplyr::n(), 
      n_alabo = dplyr::n_distinct(.data$alabo)
    )

  out <- dplyr::left_join(out, out_unique, by = "determinand")
  

  # remove data that is way off
  
  data <- dplyr::filter(
    data, 
    dplyr::between(.data$relative_u, 1, 100)
  ) 
  
  
  # relative error
  # median relative_u for values above the detection level by alabo
  
  out_relative <- data |> 
    dplyr::filter(.data$censoring == "") |> 
    dplyr::group_by(.data$determinand, .data$alabo) |> 
    dplyr::summarise(sd_variable = median(.data$relative_u) / 100) 
  
  # now the median value across alabos
  
  out_relative <- out_relative |> 
    dplyr::group_by(.data$determinand) |> 
    dplyr::summarise(sd_variable = median(.data$sd_variable))
  
  out <- dplyr::left_join(out, out_relative, by = "determinand")
  
  
  # constant error
  # median limit_detection for values with censoring == D, Q or "" by alabo
  # don't use "<" because we can't trust any of the limit values
  
  out_constant <- data |> 
    dplyr::filter(.data$censoring %in% c("D", "Q", "")) |> 
    tidyr::drop_na(.data$limit_detection) |> 
    dplyr::group_by(.data$determinand, .data$alabo) |> 
    dplyr::summarise(sd_constant = median(.data$limit_detection) / 3) 
  
  # now the median value across alabos
  
  out_constant <- out_constant |> 
    dplyr::group_by(.data$determinand) |> 
    dplyr::summarise(sd_constant = median(.data$sd_constant))
  
  out <- dplyr::left_join(out, out_constant, by = "determinand")
  
  # tidy up
  
  out <- out |> 
    as.data.frame() |> 
    tibble::column_to_rownames("determinand") |> 
    round(6)
    
  out
}

#' @export
ctsm_uncrt_plot_estimates <- function(uncrt_obj, group_id) {

  # silence non-standard evaluation warnings
  .data <- NULL
  
  data <- dplyr::filter(uncrt_obj$data, .data$group %in% group_id)
  
  data <- dplyr::arrange(data, .data$determinand, .data$concentration)
  
  data <- dplyr::filter(data, .data$relative_u >= 1 & .data$relative_u <= 100)
  
  new <- uncrt_obj$estimates[c("sd_constant", "sd_variable")]
  names(new) <- c("sdC", "sdV")
  
  var_id <- paste0(uncrt_obj$compartment, c("_sd_constant", "_sd_variable"))
  old <- uncrt_obj$old_estimates[var_id]
  names(old) <- c("sdC", "sdV")
  
  xyplot(
    relative_u ~ concentration | determinand, data = data, 
    aspect = 1,
    scales = list(
      alternating = FALSE, 
      x = list(log = TRUE, relation = "free", equispaced = FALSE)
    ), 
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

  # silence non-standard evaluation warnings
  alabo <- NULL
  
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
