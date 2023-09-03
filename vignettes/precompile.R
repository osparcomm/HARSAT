# Vignettes that are compute-intensive have been precompiled:
library(devtools)
load_all('.')

library(harsat)
library(knitr)
setwd('vignettes')
# knit("example_OSPAR.Rmd.orig", "example_OSPAR.Rmd")
# knit("example_HELCOM.Rmd.orig", "example_HELCOM.Rmd")
# knit("example_external_data.Rmd.orig", "example_external_data.Rmd")
knit("harsat.Rmd.orig", "harsat.Rmd")

setwd('..')
devtools::build_vignettes()
