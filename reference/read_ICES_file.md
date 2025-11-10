# Reads an ICES extraction

Utility function to read an ICES extraction (either data or stations)
having specified the required column names and types. Prints a warning
if required column names are missing.

## Usage

``` r
read_ICES_file(infile, var_id)
```

## Arguments

- infile:

  The input file name

- var_id:

  A named character vector giving the name and type (character, integer,
  logical, numeric, etc.) of each column in the file; these do not need
  to be in the same order as in the file.

## Value

A data frame containing the input file
