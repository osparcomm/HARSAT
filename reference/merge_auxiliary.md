# Merge auxiliary variables with data

Merge auxiliary variables with data

## Usage

``` r
merge_auxiliary(data, info, determinands)
```

## Arguments

- data:

  a data frame containing the contaminant data in long format, both the
  contaminants to be assessed and their auxiliary variables

- info:

  a harsat info object

- determinands:

  a character string given the identifiers of the determinands that are
  to be assessed

## Value

a data frame containing the contaminant data in wide format, with the
auxiliary variables pivoted to match the determinands they are linked to

## Details

`info$determinand` identifies which auxiliary variables are linked to
each determinand; `info$auxiliary` also allows the user (currently
limited) some control over how the auxiliary variables are linked to the
determinand data

Some variables can both be determinands to be assessed and auxiliary
variables (for example, CORG or DRYWT%)
