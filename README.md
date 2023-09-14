
# harsat

<!-- badges: start -->
[![lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) ![development build](https://github.com/osparcomm/harsat/actions/workflows/build.yml/badge.svg?branch=develop) ![stable build](https://github.com/osparcomm/harsat/actions/workflows/build.yml/badge.svg?branch=main)
<!-- badges: end -->

## Requirements

-	R programming language (version 4.2.1)
- Additional R packages to those which come with the standard R installation -- they will be installed automatically, but you may need permissions or tools to do that

## Installation

For full installation details, see the [Getting Started](https://osparcomm.github.io/HARSAT/articles/harsat.html) guide.
We recommend installing `harsat` from a packaged bundle, either using RStudio or the R command line, 
as this ensures that all dependencies are up-to-date, properly downloaded and installed.

### In RStudio

To do this, you will need a downloaded package bundle, typically a file 
called something like `harsat_0.1.1.tar.gz`.

From the `Packages` tab, choose `Install`, make sure you have selected
to install from a Package Archive File, then press the `Browse...` button
and locate your bundled package file. Then finally press the `Install` 
button.

### From the R command line

To do this, you will need a downloaded package bundle, typically a file 
called something like `harsat_0.1.1.tar.gz`.

In the R command line, use a command (passing the filename of wherever 
you have downloaded the file to):

```r
install.packages("~/Downloads/harsat_0.1.1.tar", repos = NULL)
```
### Directly from Github

Alternatively, you can also install `harsat` directly from GitHub. For this to work, at the moment, you will
need a Github personal access token (the "classic" kind is best), because the repository will be private until
`harsat` is finally released. As wth bundles, this ensures that all dependencies are up-to-date, properly 
downloaded and installed.

The short version is as follows:

> The `XXXX` is a Github personal access token. You only need this optional parameter while
> the `harsat` package is private because it is still under development. 
> The [Getting Started](https://osparcomm.github.io/HARSAT/articles/harsat.html) guide
> has more information on how to create a personal access token. 

``` r
library(remotes)
remotes::install_github("osparcomm/HARSAT@main", auth_token = 'XXXX')
```

The development version is similar, but changes more often, so we only recommend this if you
enjoy a more exciting time for your analysis.

``` r
library(remotes)
remotes::install_github("osparcomm/HARSAT", auth_token = 'XXXX')
```

## Example usage

We have prepared zip files containing all the other files you need, for both
OSPAR (a subset of OSPAR 2022), and HELCOM (based on HELCOM 2023). These zip files
contain two directories: a data directory and an information directory. You can
unzip these anywhere you like on your system.

* <a href="https://osparcomm.github.io/HARSAT/ospar.zip" download>Download the OSPAR archive</a>
* <a href="https://osparcomm.github.io/HARSAT/helcom.zip" download>Download the HELCOM archive</a>

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
