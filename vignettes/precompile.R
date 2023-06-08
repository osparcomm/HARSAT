# Vignettes that are compute-intensive have been precompiled:
library(devtools)
load_all('.')

library(harsat)
library(knitr)
setwd('vignettes')
knit("example_HELCOM.Rmd.orig", "example_HELCOM.Rmd")
knit("example_external_data.Rmd.orig", "example_external_data.Rmd")
knit("example_simple_OSPAR.Rmd.orig", "example_simple_OSPAR.Rmd")

setwd('..')
devtools::build_vignettes()
