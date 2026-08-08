# Create a time series

Cleans the data and turns it into time series structures ready for
assessment

## Usage

``` r
create_timeseries(
  ctsm.obj,
  determinands = ctsm_get_determinands(ctsm.obj$info),
  determinands.control = NULL,
  oddity_path = "oddities",
  return_early = FALSE,
  print_code_warnings = FALSE,
  get_basis = get_basis_default,
  normalise = FALSE,
  normalise.control = list()
)
```

## Arguments

- ctsm.obj:

  the harsat object, as returned from `tidy_data`

- determinands:

  the determinands to use, by default derived by calling
  `ctsm_get_determinands`, which takes values from the determinand
  reference table

- determinands.control:

  determinands control values, when needed

- oddity_path:

  a directory to write the oddities to

- return_early:

  a boolean (default `FALSE`), if `TRUE`, returns early

- print_code_warnings:

  a boolean (default `FALSE`)

- get_basis:

  a basis function, by default
  [`get_basis_default`](http://osparcomm.github.io/HARSAT/reference/get_basis_default.md)

- normalise:

  boolean or function, if TRUE, uses a default normalization

- normalise.control:

  a list of control data for the normalization function

## Value

a completed timeseries object, which can be used for assessments
