# Reads station dictionary

Reads station dictionary

## Usage

``` r
read_stations(file, data_dir = ".", info)
```

## Arguments

- file:

  A file reference for the station dictionary.

- data_dir:

  A path to the directory holding the station dictionary. Defaults to
  the working directory.

- info:

  A list containing at least the following two elements:

  - compartment: `"biota"`, `"sediment"` or `"water"`

  - data_format: `"ICES"` or `"external"`

## Value

A data frame containing the station dictionary.
