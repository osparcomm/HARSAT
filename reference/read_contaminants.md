# Reads contaminant data

A quick way of reading the contaminant data (which also avoids having to
do any station matching if `info$data_format == "ICES`).

## Usage

``` r
read_contaminants(file, data_dir = ".", info)
```

## Arguments

- file:

  A file reference for the contaminant data.

- data_dir:

  A path to the directory holding the contaminant data. Defaults to the
  working directory.

- info:

  A list containing at least the following two elements:

  - compartment: `"biota"`, `"sediment"` or `"water"`

  - data_format: `"ICES"` or `"external"`

  An ICES biota extraction also requires:

  - use_stage: `TRUE` or `FALSE`

## Value

A data frame containing the contaminant data.
