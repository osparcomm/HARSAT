file = 'information/assessment criteria biota.csv'

#Libraries
library(magrittr)
library(purrr)

#Make a vector of all the encodings supported by R
encodings <- set_names(iconvlist(), iconvlist())
#Make a simple reader function
reader <- function(encoding, file) {
  read.csv(file, fileEncoding = encoding, nrows = 3, header = TRUE)
}
#Create a "safe" version so we only get warnings, but errors don't stop it
# (May not always be necessary)
safe_reader <- safely(reader)

#Use the safe function with the encodings and the file being interrogated
map(encodings, safe_reader, file) %>%
  #Return just the results
  map("result") %>%
  #Keep only results that are dataframes
  keep(is.data.frame) %>%
  #Keep only results with more than one column
  #This predicate will need to change with the data
  keep(~ ncol(.x) > 1) %>%
  #Return the names of the encodings
  names()



