ctsm_check_variable <- function(data, var_id, info) {

  # location: import_check_functions.R
  # purpose: wrapper function for checking that a data variable (var_id) has 
  #   valid values, and for supplying / correcting values where possible 

  if (!(var_id %in% names(data))) {
    return(data)
  }
  
  
  # augment data with four variables: 
  # ok says whether original value is ok and should be retained
  # ok.delete says whether original value is valid but is not to be used in the  
  #   assessment
  # action says what we do - none, warning, error or delete: 
  #   initialised as NA to test whether all cases have been considered; 
  #   rows marked as error are deleted; 
  #   rows marked as delete are valid but are not to be used in the assessment 
  # new holds the revised version of var_id; this might happen even if action is
  #   none, so always check for this
  
  data_names <- names(data)
  new_names <- c("new", "ok", "ok.delete", "action")
  if (any(new_names %in% data_names)) { 
    stop("variable(s) 'new', 'ok', 'ok.delete' or 'action' already exist")
  }
  
  data$ok <- NA
  data$ok.delete <- NA
  data$action <- factor(NA, levels = c("none", "delete", "error", "warning"))
  data$new <- data[[var_id]]

  if (is.factor(data$new)) {
    data$new <- as.charcter(data$new)
  }
    

  # call checking function 
  
  data <- do.call(
    paste("ctsm.check", var_id, info$compartment, sep = "."), 
    list(data = data, info = info)
  )

  
  # output results 
  
  outfile_name <- paste0(var_id, "_queries.csv")
  outfile <- file.path(info$oddity_dir, info$compartment, outfile_name)
  

  # check all records have a valid action - if not, probably need to update
  # the check function for the variable in question
  
  if (any(is.na(data$action))) {
    oddities <- data[is.na(data$action), data_names]
    readr::write_excel_csv(oddities, outfile, na = "")
    stop(
      "Not all cases considered when checking '", var_id, "': see '", 
      out_file_name, "'\n", 
      "You might need to contact the HARSAT development team to fix this.", 
    )
  }


  # write out oddities
  
  id <- data$action %in% c("error", "warning")
  if (any(id)) {
    message(
      "  Unexpected or missing values for '", var_id, "': see ", outfile_name
    )
    oddities <- data[id, c("action", data_names)]
    readr::write_excel_csv(oddities, outfile, na = "")
  }
  
  
    
  # delete any errors or records that are not required

  data <- data[data$action %in% c("none", "warning"), ]
  
  
  # tidy up and return
    
  data[[var_id]] <- data$new
  data <- data[data_names]
    
  data
}



# NB I-RNC is for isotope ratios (used as an auxiliary) - check have correct pargroup

ctsm_is_contaminant <- function(pargroup, exclude = NULL) {
  ok <- c(
    "I-MET", "I-RNC", "O-MET", "O-BR", "O-FL", "O-PAH", "OC-CB", "OC-CL", "OC-CP", 
    "OC-DD", "OC-DN", "OC-DX", "OC-HC", "O-HER", "O-INS", "O-TRI"
  )
  ok <- setdiff(ok, exclude)
  pargroup %in% ok
}  


ctsm.check.basis.sediment <- function(data, info) {

  id <- ctsm_is_contaminant(data$pargroup) | 
    data$determinand %in% c("CORG", "LOIGN")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- basis %in% c("D", "W")
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "DRYWT%"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- basis %in% "W"
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "W"
    })
  
  data               
}

ctsm.check.basis.biota <- function(data, info) {

  id <- ctsm_is_contaminant(data$pargroup)
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- basis %in% c("D", "W", "L")
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% c("LNMEA", "AGMEA")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- is.na(basis)
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- NA
    })
  
  id <- data$determinand %in% "DRYWT%"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- basis %in% "W"
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "W"
    })

  id <- data$determinand %in% c("EXLIP%", "FATWT%", "LIPIDWT%")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- basis %in% c("W", "D")
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "W"
    })

  id <- data$determinand %in% c(
    "VDS", "IMPS", "INTS", "VDSI", "PCI", "INTSI", "%FEMALEPOP", "SURVT", "NRR", "LP", "%DNATAIL", 
    "MNC", "CMT-QC-NR", "MNC-QC-NR")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- is.na(basis)
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- NA
    })
  
  id <- data$group %in% "Metabolites" | 
    data$determinand %in% c("ALAD", "EROD", "SFG", "ACHE", "GST")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- basis %in% c("W", NA)
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- NA
    })
  
  data             
}             

ctsm.check.basis.water <- function(data, info) {
  
  data <- within(data, {
    ok <- basis %in% "W"
    action <- ifelse(ok, "none", ifelse(basis %in% NA, "warning", "error"))
    new[basis %in% NA] <- "W"
  })
  
  data               
}


ctsm.check.matrix.sediment <- function(data, info) {
  
  data <- within(data, {
    ok <- substr(matrix, 1, 3) %in% "SED"
    action <- ifelse(ok, "none", "error")
  })
  
  # rationalise sediment matrices: necessary since e.g. some measure organics in
  # SEDTOT and CORG in SED2000, which are effectively the same thing
  
  if (any(data$matrix %in% c("SED62", "SED500", "SED1000", "SED2000"))) {
    cat("   Relabelling matrix SED62 as SED63 and SED500, SED1000, SED2000 as SEDTOT\n")
    data <- within(data, {
      new[matrix %in% "SED62"] <- "SED63"
      new[matrix %in% c("SED500", "SED1000", "SED2000")] <- "SEDTOT"
    })
  }

  if (!all(data$new %in% c("SED20", "SED63", "SEDTOT")))
    stop("Unrecognised sediment matrix")
  
  data               
}

ctsm.check.matrix.biota <- function(data, info) {
  
  id <- ctsm_is_contaminant(data$pargroup) | 
    data$determinand %in% c("DRYWT%", "EXLIP%", "FATWT%", "LIPIDWT%")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- (species_group %in% "Fish" & matrix %in% c("MU", "LI", "MU&EP")) |
        (species_group %in% c("Bivalve", "Gastropod") & matrix %in% "SB") |
        (species_group %in% "Crustacean" & matrix %in% "TM") | 
        (species_group %in% "Bird" & matrix %in% c("EH", "FE", "LI", "MU", "BL", "ER")) | 
        (species_group %in% "Mammal" & matrix %in% c("BB", "HA", "KI", "LI", "MU", "EP"))
      change <- species_group %in% "Bird" & matrix %in% "EG"
      action <- ifelse(ok, "none", ifelse(change, "warning", "error"))
      new[change] <- "EH"
      rm(change)
    })
  
  # some ambiguity here about LNMEA for birds - could be WO (LI, MU, BL, ER, FE) or ES (EG, EH), so if LNMEA 
  # matrix is not one of these, throw an error
  # not clear if matrix should be ES, SH or EG for eggs, need to get clarification
  # actually LNMEA could be feather length (but have not allowed for this at present)
  # NB procedures for merging with LNMEA are similary complicated for birds

  id <- data$determinand %in% "LNMEA"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- (species_group %in% c("Fish", "Mammal") & matrix %in% "WO") |
        (species_group %in% "Bird" & matrix %in% c("WO", "ES")) |
        (species_group %in% c("Bivalve", "Gastropod", "Crustacean") & matrix %in% "SH")
      action <- ifelse(
        ok, "none", ifelse(
          species_group %in% "Bird" & ! (matrix %in% c("LI", "MU", "BL", "ER", "FE", "EH", "EG", "SH")), 
          "error", "warning"))
      new[!ok] <- ifelse(
        species_group[!ok] %in% c("Fish", "Mammal"), "WO", ifelse(
          species_group[!ok] %in% "Bird" & matrix[!ok] %in% c("LI", "MU", "BL", "ER", "FE"), "WO", ifelse(
            species_group[!ok] %in% "Bird" & matrix[!ok] %in% c("EG", "EH", "SH"), "ES", ifelse(
              species_group[!ok] %in% "Bird", NA, "SH"))))
    })

  # could maybe measure age of eggs as well, but have assumed always WO
  
  id <- data$determinand %in% "AGMEA"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- matrix %in% "WO"
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "WO"
    })

  id <- data$group %in% "Imposex"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% "SB"
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "SB"
    })
  
  id <- data$determinand %in% "%FEMALEPOP"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% "POP"
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "POP"
    })
  
  id <- data$determinand %in% "EROD"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("LIMIC", "LIS9")
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% "LP"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("LI")
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "LI"
    })
  
  id <- data$determinand %in% c("NAP2OH", "PYR1OH", "PYR1OHEQ", "PA1OH", "BAP3OH")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("BI")
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "ALAD"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("BL")
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "BL"
    })

  id <- data$determinand %in% c("MNC", "MNC-QC-NR")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("ER")
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "ER"
    })

  id <- data$determinand %in% c("SFG")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% "SB"
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "SB"
    })

  id <- data$determinand %in% "ACHE"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- (species_group %in% "Fish" & matrix %in% c("MU", "BR")) |
       (species_group %in% "Bivalve" & matrix %in% "GI")  
      action <- ifelse(
        ok, "none", 
        ifelse(species_group %in% "Fish", "error", "warning")
      )
      new[!ok & species_group %in% "Bivalve"] <- "GI" 
    })

  id <- data$determinand %in% "GST"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- (species_group %in% "Fish" & matrix %in% "LICYT") |
        (species_group %in% "Bivalve" & matrix %in% "SB")  
      action <- ifelse(
        ok, "none", 
        ifelse(species_group %in% "Fish", "error", "warning")
      )
      new[!ok & species_group %in% "Bivalve"] <- "SB" 
    })
  
  
  id <- data$determinand %in% c("%DNATAIL", "CMT-QC-NR", "NRR")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("HML")
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "HML"
    })

  id <- data$determinand %in% c("SURVT")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- matrix %in% c("WO")
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- "WO"
    })
  
  data             
}             

ctsm.check.matrix.water <- function(data, info) {
  
  data <- within(data, {
    ok <- matrix %in% c("WT", "AF", "BF")
    action <- ifelse(ok, "none", "error")
    new[matrix %in% c("AF", "BF")] <- "WT"
  })

  data               
}



ctsm.check.species_group.biota <- function(data, info) {  
  
  # check species_group appropriate for each determinand
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "O-PAH") | data$group %in% "Auxiliary"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- TRUE
      action <- "none"
    })
  
  id <- data$pargroup %in% "O-PAH"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !species_group %in% "Fish"
      action <- ifelse(ok, "none", "error")    
    })

  id <- data$group %in% "Imposex"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- species_group %in% "Gastropod"
      action <- ifelse(ok, "none", "error")
    })    
    
  id <- data$group %in% "Metabolites" | data$determinand %in% c("EROD", "ALAD", "LP", "MNC")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- species_group %in% "Fish"
      action <- ifelse(ok, "none", "error")
    })
       
  id <- data$determinand %in% c("ACHE", "GST")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- species_group %in% c("Bivalve", "Fish")
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("SFG", "%DNATAIL", "NRR", "SURVT")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- species_group %in% "Bivalve"
      action <- ifelse(ok, "none", "error")
    })
  
  data             
}             

ctsm.check.sex.biota <- function(data, info) {

  # NB any changes should really be made at the sub-sample level
  
  id <- ctsm_is_contaminant(data$pargroup) | 
    data$group %in% "Metabolites" | 
    data$determinand %in% c("AGMEA", "LNMEA", "DRYWT%", "EXLIP%", "FATWT%", "LIPIDWT%") | 
    data$determinand %in% c("ALAD", "SFG", "ACHE", "GST", "SURVT", "NRR", "LP", "%DNATAIL", 
                            "MNC", "CMT-QC-NR", "MNC-QC-NR")

  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- sex %in% c("F", "I", "M", "U", "X", NA)
      action <- ifelse(ok, "none", "warning")
      new[!ok] <- NA
    })


  id <- data$determinand %in% c("VDS", "IMPS", "INTS")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- sex %in% "F"
      action <- ifelse(ok, "none", ifelse(sex %in% NA, "warning", "error"))
      new[sex %in% NA] <- "F"
    })
  
  id <- data$determinand %in% c("VDSI", "PCI", "INTSI", "%FEMALEPOP")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- sex %in% c("X", "F", NA)
      action <- ifelse(ok, "none", "error")
      new[sex %in% NA] <- "X"
    })

  id <- data$determinand %in% "EROD"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- sex %in% c("F", "M")
      ok.delete <- sex %in% c("U", "I", "X")
      action <- ifelse(ok, "none", ifelse(ok.delete, "delete", "error"))
      if (any(ok.delete))
        cat("   Dropping EROD data with immature or unidentifiable sex\n")
    })

  data             
}             


ctsm.check.unit.biota <- function(data, info) {

  standard_unit <- c(
    "g/g", "mg/mg", "ug/ug", "ng/ng", "pg/pg", "mg/g", "ug/g", "ng/g", "pg/g", "g/kg", "mg/kg", 
    "ug/kg", "ng/kg", "pg/kg")
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "I-RNC") & data$determinand != "TEQDFP"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% standard_unit
      action <- ifelse(ok, "none", "error")
    })
    
  id <- data$determinand %in% "TEQDFP"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% paste("TEQ", standard_unit)
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("LNMEA")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% c("km", "m", "cm", "mm")
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("AGMEA")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% "y"
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% c("DRYWT%", "EXLIP%", "FATWT%", "LIPIDWT%")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "%"
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("C13D", "N15D")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "ppt"
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("VDS", "IMPS", "INTS")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "st"
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% c("VDSI", "PCI", "INTSI")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "idx"
      action <- ifelse(ok, "none", "error")
      message("   imposex index units changed from 'idx' to 'st' before merging with individual data")
      new[ok] <- "st"
    })
  
  id <- data$determinand %in% c("%FEMALEPOP", "%DNATAIL")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "%"
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "EROD"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "pmol/min/mg protein"
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$group %in% "Metabolites"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% c("ng/g", "ng/ml", "ug/kg", "ug/l", "ug/ml")
      #    ok.delete <- unit %in% "ng g-1 a660-1"
      ok.delete <- FALSE
      action <- ifelse(ok, "none", ifelse(ok.delete, "delete", "error"))
      
      if (any(ok.delete))
        cat("   Dropping bile metabolites with unit 'ng g-1 a660-1'\n")
      
      if (any(unit %in% c("ng/g", "ug/kg"))) {
        cat ("   Bile metabolite units changed from 'ng/g' to 'ng/ml' and from",
             "'ug/kg' to 'ug/l'\n")
        new[unit %in% "ng/g"] <- "ng/ml"
        new[unit %in% "ug/kg"] <- "ug/l"
      }
    })

  id <- data$determinand %in% "ALAD"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "ng/min/mg protein"
      action <- ifelse(ok, "none", "error")
    })

  id <- with(data, determinand %in% c("ACHE", "GST"))
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% c("umol/min/mg protein", "nmol/min/mg protein", "pmol/min/mg protein")
      action <- ifelse(ok, "none", "error")
    })
                
  id <- with(data, determinand %in% "SFG")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "j/h/g"
      action <- ifelse(ok, "none", "error")
    })

  id <- with(data, determinand %in% "SURVT")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "d"
      action <- ifelse(ok, "none", "error")
    })
  
  id <- with(data, determinand %in% c("LP", "NRR"))
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "mins"
      action <- ifelse(ok, "none", "error")
    })

  id <- with(data, determinand %in% c("CMT-QC-NR", "MNC-QC-NR"))
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "nr"
      action <- ifelse(ok, "none", "error")
    })
  
  id <- with(data, determinand %in% "MNC")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- unit %in% "nr/1000 cells"
      action <- ifelse(ok, "none", "error")
    })

  data             
}             

ctsm.check.unit.sediment <- function(data, info) {

  standard_unit <- c(
    "g/g", "mg/mg", "ug/ug", "ng/ng", "pg/pg", "mg/g", "ug/g", "ng/g", "pg/g", "g/kg", "mg/kg", 
    "ug/kg", "ng/kg", "pg/kg")
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "I-RNC") & !data$determinand %in% c("AL", "LI")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% standard_unit
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% c("AL", "LI", "CORG", "LOIGN")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% c(standard_unit, "%")
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "DRYWT%"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% "%"
      action <- ifelse(ok, "none", "error")
    })
  
  data             
}             

ctsm.check.unit.water <- function(data, info) {
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "I-RNC")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% c("mg/l", "ug/l", "ng/l", "pg/l")
      action <- ifelse(ok, "none", "error")
  })

  data             
}             


ctsm.check.method_analysis.sediment <- function(data, info) {

  data <- within(data, {
    ok <- TRUE
    action <- "none"
  })
  
  data
}  

ctsm.check.method_analysis.water <- function(data, info) {
  
  data <- within(data, {
    ok <- TRUE
    action <- "none"
  })
  
  data
}  

ctsm.check.method_analysis.biota <- function(data, info) {

  id <- data$group != "Metabolites"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- TRUE
      action <- "none"
    })
  
  id <- data$group %in% "Metabolites"
  if (any(id))
    data[id,] <- within(data[id,], {    
      # ok <- method_analysis %in% c("FLM-SS", "HPLC-FD", "GC-MS", "GC-MS-MS", "GC-MS-SIM")
      ok <- !is.na(method_analysis)
      action <- ifelse(ok, "none", "error")
      
      if (any(method_analysis %in% c("GC-MS-MS", "GC-MS-SIM"))) {
        cat ("   Bile metabolite method_analysis GC-MS-MS and GC-MS-SIM changed to GC-MS\n")
        new[method_analysis %in% c("GC-MS-MS", "GC-MS-SIM")] <- "GC-MS"
      }
    })
  
  data             
}             


ctsm.check.value.biota <- function(data, info) {

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  
    abs(x - round(x)) < tol
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "I-MTC") | 
    data$group %in% "Metabolites" | 
    data$determinand %in% c("EROD", "ALAD", "ACHE", "GST", "AGMEA", "LNMEA")

  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value > 0
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "SFG"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value)
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("DRYWT%", "EXLIP%", "FATWT%", "LIPIDWT%")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value > 0 & value < 100
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("C13D", "N15D")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) 
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("%FEMALEPOP", "%DNATAIL")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value >= 0 & value <= 100
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% c("VDSI", "PCI", "INTSI")
  if (any(id)) {
    min.imposex <- with(
      data[id,], 
      get.info.imposex(species, determinand, info$imposex, "min_value", na.action = "ok")
    )
    max.imposex <- with(
      data[id,], 
      get.info.imposex(species, determinand, info$imposex, "max_value", na.action = "ok")
    )
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & !is.na(min.imposex) & !is.na(max.imposex) &
        value >= min.imposex & value <= max.imposex
      ok.delete <- !is.na(value) & (is.na(min.imposex) | is.na(max.imposex))
      action <- ifelse(ok, "none", ifelse(ok.delete, "delete", "error"))
    })
  }
  
  id <- data$determinand %in% c("VDS", "IMPS", "INTS")
  if (any(id)) {
    min.imposex <- with(
      data[id,], 
      get.info.imposex(species, determinand, info$imposex, "min_value", na.action = "ok")
    )
    max.imposex <- with(
      data[id,], 
      get.info.imposex(species, determinand, info$imposex, "max_value", na.action = "ok")
    )
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & !is.na(min.imposex) & !is.na(max.imposex) &
        value >= min.imposex & value <= max.imposex & is.wholenumber(value)
      ok.delete <- !is.na(value) & (is.na(min.imposex) | is.na(max.imposex))
      action <- ifelse(ok, "none", ifelse(ok.delete, "delete", "error"))
    })
  }
  
  id <- data$determinand %in% c("SURVT", "CMT-QC-NR", "MNC-QC-NR")
  if (any(id)) {
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & is.wholenumber(value) & value >= 1 
      action <- ifelse(ok, "none", "error")
    })
  }

  id <- data$determinand %in% "MNC"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value >= 0 & value <= 1000
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "NRR"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- value %in% c(0, 15, 30, 60, 90, 120, 150, 180)
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% "LP"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- value %in% c(0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50)
      action <- ifelse(ok, "none", "error")
    })

  data             
}             

ctsm.check.value.sediment <- function(data, info) {

  id <- ctsm_is_contaminant(data$pargroup) & !data$determinand %in% c("AL", "LI")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value > 0
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% c("AL", "LI", "CORG", "LOIGN")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value > 0 
      ok <- ifelse(unit %in% "%", ok & value < 100, ok)
      action <- ifelse(ok, "none", "error")
    })
  
  id <- data$determinand %in% "DRYWT%"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & value > 0 & value < 100
      action <- ifelse(ok, "none", "error")
    })
  
  data             
}             

ctsm.check.value.water <- function(data, info) {
  
  data <- within(data, {
    ok <- !is.na(value) & value > 0
    action <- ifelse(ok, "none", "error")
  })

  data             
}             


ctsm.check.n_individual.biota <- function(data, info) {

  id <- data$group != "Imposex"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- is.na(n_individual) | n_individual > 0L
      action <- ifelse(ok, "none", "warning")
    })
  
  id <- data$determinand %in% c("VDS", "IMPS", "INTS")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- !is.na(n_individual) & n_individual == 1L
      action <- ifelse(ok, "none", ifelse(is.na(n_individual), "warning", "error"))
      new[is.na(n_individual)] <- 1L
    })
    
  id <- data$determinand %in% c("VDSI", "PCI", "INTSI")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- !is.na(n_individual) & n_individual >= 1L
      action <- ifelse(ok, "none", "error")
    })
  
  data             
}             
