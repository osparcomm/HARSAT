
# harsat

## Requirements

- R programming language (version 4.2.1)
- Additional R packages to those which come with the standard R installation -- they will be installed automatically, but you may need permissions or tools to do that

## Installation

For full installation details, see the [Getting Started](./articles/harsat.html) guide.
We strongly recommend installing `harsat` from GitHub, as this ensures that all dependencies are up-to-date, properly 
downloaded and installed.

The short version is as follows:

> The `XXXX` is a Github personal access token. You only need this optional parameter while
> the `harsat` package is private because it is still under development. 
> The [Getting Started](./articles/harsat.html) guide
> has more information on how to create a personal access token. 

``` r
library(remotes)
remotes::install_github("osparcomm/HARSAT", auth_token = 'XXXX')
```

The stable version is similar:

``` r
library(remotes)
remotes::install_github("osparcomm/HARSAT@main", auth_token = 'XXXX')
```

## Example usage

We have prepared zip files containing all the other files you need, for both
OSPAR (a subset of OSPAR 2022), and HELCOM (based on HELCOM 2023). These zip files
contain two directories: a data directory and an information directory. You can
unzip these anywhere you like on your system.

* <a href="./ospar.zip" download>Download the OSPAR archive</a>
* <a href="./helcom.zip" download>Download the HELCOM archive</a>

Let's say you have unzipped these to a file on your system, such as `C:\Users\test\ospar`. So, 
the data is now in `C:\Users\test\ospar\data`, and the reference files are in `C:\Users\test\ospar\information`,
although can rename and move these anywhere you like. 

To read the data, you will then do something like this (we're using water in OSPAR as an example here -- 
for the complete example, have a look at [the full OSPAR example](./articles/example_OSPAR.html)):

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

And follow the rest of the process as shown in the [Getting Started guide](./articles/harsat.html) or the
[OSPAR example](./articles/example_OSPAR.html).

You can, of course, start to edit the files in these directories as you choose. Note that there are
some specific naming conventions for the files in the reference file directory especially. 
Find out more in [the documentation page for file management](./articles/file-management.html).

## More information

For more information, take a look at the [Getting Started guide](./articles/harsat.html).

We welcome any other contributions you can make. Check out the [Contributor's guide](./CONTRIBUTING.html) for more.
