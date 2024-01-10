
# harsat

<!-- badges: start -->
[![lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) ![development build](https://github.com/osparcomm/harsat/actions/workflows/build.yml/badge.svg?branch=develop) ![stable build](https://github.com/osparcomm/harsat/actions/workflows/build.yml/badge.svg?branch=main)
<!-- badges: end -->

## What is HARSAT?

HARSAT (Harmonised Regional Seas Assessment Tool), is a tool that is applied 
by The Arctic Monitoring and Assessment Programme (AMAP), the Helsinki 
Commission (HELCOM) and the OSPAR Commission (OSPAR) to support their 
assessments of data concerning contaminants (hazardous substances) and 
their effects in the marine environment. 

HARSAT code includes tools for pre-processing data, statistical trend 
analysis and comparison with threshold values, and post-processing for 
archiving and reporting.

HARSAT is developed in statistical computing language R. R is available 
for most operating systems and can be downloaded from the R-project website. 

Disclaimer: The HARSAT tool is made available under an Open Source licence. 
While every attempt has been made to ensure that the HARSAT version(s) 
developed and supported by AMAP/HELCOM/OSPAR are free of errors, any use 
of the tool by third parties, and the quality of products of third party 
use is the responsibility of the third party concerned.

## System requirements

-	R programming language (version 4.2.1 or later). Additional R packages 
to those which come with the standard R installation may be 
needed -- they will normally be installed automatically, but you may need 
permissions or tools to do that

- RStudio (version 2023.03.1 or later, recommended). The HARSAT developers recommend 
running the HARSAT code using the RStudio integrated development 
environment. Although R can be run independently of RStudio, some of the 
examples presented are easier within RStudio. 

## Installation

For full installation details, see the [Getting Started](https://osparcomm.github.io/HARSAT/articles/harsat.html) guide.
We recommend installing `harsat` from a packaged bundle, either using RStudio or the R command line, 
as this ensures that all dependencies are up-to-date, properly downloaded and installed.

### In RStudio

To do this, you will need a downloaded package bundle, typically a file 
called something like `harsat_0.1.2.tar.gz`.

From the `Packages` tab, choose `Install`, make sure you have selected
to install from a Package Archive File, then press the `Browse...` button
and locate your bundled package file. Then finally press the `Install` 
button.

### From the R command line

To do this, you will need a downloaded package bundle, typically a file 
called something like `harsat_0.1.2.tar.gz`.

In the R command line, use a command (passing the filename of wherever 
you have downloaded the file to):

```r
install.packages(remotes) -- if needed
library(remotes)
remotes::install_local("~/Downloads/harsat_0.1.2.tar")
```
### Directly from Github

Alternatively, you can also install `harsat` directly from GitHub. 
``` r
library(remotes)
remotes::install_github("osparcomm/HARSAT@main")
```

The development version is similar, but changes more often, so we only recommend this if you
enjoy a more exciting time for your analysis.

``` r
library(remotes)
remotes::install_github("osparcomm/HARSAT", auth_token = 'XXXX')
```

## Additional funding

Additional resource that contributed to the development of HARSAT from: 
[Baltic Data Flows project](https://balticdataflows.helcom.fi/),
the [PreEMPT project](https://helcom.fi/helcom-at-work/projects/pre-empting-pollution-by-screening-for-possible-risks-preempt/),
Nordic Council of Ministers (Harmonized Regional Seas Assessment Tool: NKE 2023-011) 
and via the standing working process OSPAR has with ICES.

## Example usage

We have prepared zip files containing all the other files you need, for both
OSPAR (a subset of OSPAR 2022), and HELCOM (based on HELCOM 2023). These zip files
contain two directories: a data directory and an information directory. You can
unzip these anywhere you like on your system.

* <a href="https://osparcomm.github.io/HARSAT/ospar.zip" download>Download the OSPAR archive</a>
* <a href="https://osparcomm.github.io/HARSAT/helcom.zip" download>Download the HELCOM archive</a>
* <a href="https://osparcomm.github.io/HARSAT/external.zip" download>Download the external data (AMAP) archive</a>

Let's say you have unzipped these to a file on your system, such as `C:\Users\test\ospar`. So, 
the data is now in `C:\Users\test\ospar\data`, and the reference files are in `C:\Users\test\ospar\information`,
although can rename and move these anywhere you like. 

To read the data, you will then do something like this (we're using water in OSPAR as an example here -- 
for the complete example, have a look at [the full OSPAR example](https://osparcomm.github.io/HARSAT/articles/example_OSPAR.html)):

```r
water_data <- read_data(
  compartment = "water", 
  purpose = "OSPAR",                               
  contaminants = "water.txt", 
  stations = "stations.txt", 
  data_dir = file.path("C:", "users", "test", "ospar", "data"),         ## i.e., C:\Users\test\ospar\data
  info_dir = file.path("C:", "users", "test", "ospar", "information"),  ## i.e., C:\Users\test\ospar\information
  extraction = "2023/08/23"
)
```

And follow the rest of the process as shown in the [Getting Started guide](https://osparcomm.github.io/HARSAT/articles/harsat.html) or the
[OSPAR example](https://osparcomm.github.io/HARSAT/articles/example_OSPAR.html).

You can, of course, start to edit the files in these directories as you choose. Note that there are
some specific naming conventions for the files in the reference file directory especially. 
Find out more in [the documentation page for file management](https://osparcomm.github.io/HARSAT/articles/file-management.html).

## More information

For more information, take a look at the [Getting Started guide](https://osparcomm.github.io/HARSAT/articles/harsat.html).

We welcome any other contributions you can make. Check out the [Contributor's guide](https://osparcomm.github.io/HARSAT/CONTRIBUTING.html) for more.
