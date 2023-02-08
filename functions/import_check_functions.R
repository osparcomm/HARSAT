

ctsm.check0 <- function(data, type, compartment, message, fileName) {

  # wrapper function for checking routines
  # type is the variable within data that is going to be checked
  # compartment is one of biota, sediment
 

  if (!(type %in% names(data))) return(data)
  
  
  # augment data with three variables: 
  # ok says whether original value is ok and should be retained
  # ok.delete says whether original value is correct but is not suitable for assessment
  # action says what we do - none, warning, error or delete: initialised as NA to 
  #   allow test of whether all cases have been considered (see ctsm.check2)
  # delete means data are correct, but not suitable for inclusion in assessment
  # new holds the revised version of type
  
  newNames <- c("new", "ok", "ok.delete", "action")
  if (any(newNames %in% names(data))) 
    stop("variable(s) 'new', 'ok', 'ok.delete' or 'action' already exist")
  
  data <- within(data, {
    ok <- logical()
    ok.delete <- logical()
    action <- factor(NA, levels = c("none", "delete", "error", "warning"))
  })
  
  data[["new"]] <- data[[type]]
  typeIsFactor <- is.factor(data[[type]])
  if (typeIsFactor) data <- within(data, new <- as.character(new))
  
  
  # call checking function 
  
  data <- do.call(paste("ctsm.check", type, compartment, sep = "."), list(data = data))
  

  # check all cases considered, print out data that have a warning or error and
  # delete data if required
  
  summary <- data[newNames]
  data <- data[setdiff(names(data), newNames)]
  
  if (missing(fileName)) fileName <- paste(type, "queries")
  if (missing(message)) message <- paste("Unexpected or missing values for", type)
  
  data <- ctsm.check2(data, summary$action, message, fileName)

  
  # now replace type with new and return

  data[[type]] <- with(summary, {
    out <- new[action %in% c("none", "warning")]
    if (typeIsFactor) out <- factor(out)
    out
  })
  
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


ctsm.check.basis.sediment <- function(data) {

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

ctsm.check.basis.biota <- function(data) {

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

ctsm.check.basis.water <- function(data) {
  
  data <- within(data, {
    ok <- basis %in% "W"
    action <- ifelse(ok, "none", ifelse(basis %in% NA, "warning", "error"))
    new[basis %in% NA] <- "W"
  })
  
  data               
}


ctsm.check.matrix.sediment <- function(data) {
  
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

ctsm.check.matrix.biota <- function(data) {
  
  id <- ctsm_is_contaminant(data$pargroup) | 
    data$determinand %in% c("DRYWT%", "EXLIP%", "FATWT%", "LIPIDWT%")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- (family %in% "Fish" & matrix %in% c("MU", "LI", "MU&EP")) |
        (family %in% c("Bivalvia", "Gastropoda") & matrix %in% "SB") |
        (family %in% "Crustacea" & matrix %in% "TM") | 
        (family %in% "Bird" & matrix %in% c("EH", "FE", "LI", "MU", "BL", "ER")) | 
        (family %in% "Mammal" & matrix %in% c("BB", "HA", "KI", "LI", "MU", "EP"))
      change <- family %in% "Bird" & matrix %in% "EG"
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
      ok <- (family %in% c("Fish", "Mammal") & matrix %in% "WO") |
        (family %in% "Bird" & matrix %in% c("WO", "ES")) |
        (family %in% c("Bivalvia", "Gastropoda", "Crustacea") & matrix %in% "SH")
      action <- ifelse(
        ok, "none", ifelse(
          family %in% "Bird" & ! (matrix %in% c("LI", "MU", "BL", "ER", "FE", "EH", "EG", "SH")), 
          "error", "warning"))
      new[!ok] <- ifelse(
        family[!ok] %in% c("Fish", "Mammal"), "WO", ifelse(
          family[!ok] %in% "Bird" & matrix[!ok] %in% c("LI", "MU", "BL", "ER", "FE"), "WO", ifelse(
            family[!ok] %in% "Bird" & matrix[!ok] %in% c("EG", "EH", "SH"), "ES", ifelse(
              family[!ok] %in% "Bird", NA, "SH"))))
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
      ok <- (family %in% "Fish" & matrix %in% c("MU", "BR")) |
       (family %in% "Bivalvia" & matrix %in% "GI")  
      action <- ifelse(
        ok, "none", 
        ifelse(family %in% "Fish", "error", "warning")
      )
      new[!ok & family %in% "Bivalvia"] <- "GI" 
    })

  id <- data$determinand %in% "GST"
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- (family %in% "Fish" & matrix %in% "LICYT") |
        (family %in% "Bivalvia" & matrix %in% "SB")  
      action <- ifelse(
        ok, "none", 
        ifelse(family %in% "Fish", "error", "warning")
      )
      new[!ok & family %in% "Bivalvia"] <- "SB" 
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

ctsm.check.matrix.water <- function(data) {
  
  data <- within(data, {
    ok <- matrix %in% c("WT", "AF", "BF")
    action <- ifelse(ok, "none", "error")
    new[matrix %in% c("AF", "BF")] <- "WT"
  })

  data               
}


ctsm.check.species.biota <- function(data) {

  # only assess species for which info.species$assess is TRUE
  
  data <- within(data, {
    ok <- ctsm_get_info("species", species, "assess")
    action <- ifelse(ok, "none", "delete")
  })

  with(data, if (!all(ok)) 
    cat("   Dropping following species:", 
        paste(unique(species[!ok]), collapse = ", "), "\n"))

  data
}
  
ctsm.check.family.biota <- function(data) {  
  
  # check family appropriate for each determinand
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "O-PAH") | data$group %in% "Auxiliary"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- TRUE
      action <- "none"
    })
  
  id <- data$pargroup %in% "O-PAH"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- !family %in% "Fish"
      action <- ifelse(ok, "none", "error")    
    })

  id <- data$group %in% "Imposex"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- family %in% "Gastropoda"
      action <- ifelse(ok, "none", "error")
    })    
    
  id <- data$group %in% "Metabolites" | data$determinand %in% c("EROD", "ALAD", "LP", "MNC")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- family %in% "Fish"
      action <- ifelse(ok, "none", "error")
    })
       
  id <- data$determinand %in% c("ACHE", "GST")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- family %in% c("Bivalvia", "Fish")
      action <- ifelse(ok, "none", "error")
    })

  id <- data$determinand %in% c("SFG", "%DNATAIL", "NRR", "SURVT")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- family %in% "Bivalvia"
      action <- ifelse(ok, "none", "error")
    })
  
  data             
}             

ctsm.check.sex.biota <- function(data) {

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


ctsm.check.unit.biota <- function(data) {

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

ctsm.check.unit.sediment <- function(data) {

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

ctsm.check.unit.water <- function(data) {
  
  id <- ctsm_is_contaminant(data$pargroup, exclude = "I-RNC")
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- unit %in% c("mg/l", "ug/l", "ng/l", "pg/l")
      action <- ifelse(ok, "none", "error")
  })

  data             
}             


ctsm.check.method_analysis.sediment <- function(data) {

  data <- within(data, {
    ok <- TRUE
    action <- "none"
  })
  
  data
}  

ctsm.check.method_analysis.water <- function(data) {
  
  data <- within(data, {
    ok <- TRUE
    action <- "none"
  })
  
  data
}  

ctsm.check.method_analysis.biota <- function(data) {

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


ctsm.check.value.biota <- function(data) {

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
    min.imposex <- with(data[id,], get.info.imposex(species, determinand, "min_value", na.action = "ok"))
    max.imposex <- with(data[id,], get.info.imposex(species, determinand, "max_value", na.action = "ok"))
    data[id,] <- within(data[id,], {
      ok <- !is.na(value) & !is.na(min.imposex) & !is.na(max.imposex) &
        value >= min.imposex & value <= max.imposex
      ok.delete <- !is.na(value) & (is.na(min.imposex) | is.na(max.imposex))
      action <- ifelse(ok, "none", ifelse(ok.delete, "delete", "error"))
    })
  }
  
  id <- data$determinand %in% c("VDS", "IMPS", "INTS")
  if (any(id)) {
    min.imposex <- with(data[id,], get.info.imposex(species, determinand, "min_value", na.action = "ok"))
    max.imposex <- with(data[id,], get.info.imposex(species, determinand, "max_value", na.action = "ok"))
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

ctsm.check.value.sediment <- function(data) {

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

ctsm.check.value.water <- function(data) {
  
  data <- within(data, {
    ok <- !is.na(value) & value > 0
    action <- ifelse(ok, "none", "error")
  })

  data             
}             


ctsm.check.noinp.biota <- function(data) {

  id <- data$group != "Imposex"
  if (any(id))
    data[id,] <- within(data[id,], {
      ok <- is.na(noinp) | noinp > 0
      action <- ifelse(ok, "none", "warning")
    })
  
  id <- data$determinand %in% c("VDS", "IMPS", "INTS")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- !is.na(noinp) & noinp == 1
      action <- ifelse(ok, "none", ifelse(is.na(noinp), "warning", "error"))
      new[is.na(noinp)] <- 1
    })
    
  id <- data$determinand %in% c("VDSI", "PCI", "INTSI")
  if (any(id))
    data[id,] <- within(data[id,], {    
      ok <- !is.na(noinp) & noinp >= 1
      action <- ifelse(ok, "none", "error")
    })
  
  data             
}             
