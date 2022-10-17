# 12/07/16 pastOmitNA extend to trim spaces of more than two blanks
# 08/12/17 pastOmitNA copes with strings with two spaces


swap.names <- function(data, old.names, new.names) {
  stopifnot(length(old.names) == length(new.names))
  pos <- match(old.names, names(data))
  if (any(is.na(pos))) stop(paste(old.names[is.na(pos)], collapse = ", "), " are not valid names")
  names(data)[pos] <- new.names
  if (any(duplicated(names(data)))) stop("renaming variables has created duplicate column names")
  data
}



convert.html.characters <- function(x) {
  for (i in 1:nrow(info.html)) x <- gsub(info.html$value[i], info.html$html1[i], x, fixed = TRUE)
  x
}


pasteOmitNA <- function(..., sep = " ") {
  
  # paste but omit missing characters (rather than having "NA")
  # if everything is missing, will end up with "", rather than NA
  
  # convert to character
  # replace missing values by null character string, 
  # check separator doesn't begin, end or have double values in middle of any of the arguments 
  #   (treat "" as special case) - needed to ensure that the gsub commands below don't change 
  #   original arguments
  # paste

  # use a dummy separator because some stations (Sweden!!!!) have double spaces
  
  true_sep <- sep
  
  sep <- "____"
    
  sepDouble <- paste0(sep, sep)
  sepEnd <- paste0("(^", sep, "|", sep, "$)")
  
  L <- list(...)
  L <- lapply(L,function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    if (sep != "" && any(grepl(sepDouble, x) | grepl(sepEnd, x)))
      stop("choice of separator means that pasting will change original strings")
    x
  })
  out <- do.call(paste, c(L, list(sep=sep)))
  
  if (sep == "") return(out)
  
  # replace double occurences of sep (due to missing elements) with single occurrence
  
  while (any(grepl(sepDouble, out))) out <- gsub(sepDouble, sep, out)
  
  # replace separator with " "
  
  out <- gsub(sep, true_sep, out)
  
  # trim leading and trailing blanks
  
  out <- trimws(out)
  
  out
}


stripBlankRows <- function(x, message = TRUE) {
  y <- is.na(x)
  ok <- ! apply(y, 1, all)
  if (sum(ok) == nrow(y)) 
    return(x)
  
  message("dropping ", nrow(y) - sum(ok), " blank rows")
  x[ok, ]
}

