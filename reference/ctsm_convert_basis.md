# Basis and matrix information and basis conversion

Converts concentrations between wet, dry and lipid bases

## Usage

``` r
ctsm_convert_basis(
  conc,
  from,
  to,
  drywt = NA_real_,
  lipidwt = NA_real_,
  drywt_censoring = ifelse(is.na(drywt), NA_character_, ""),
  lipidwt_censoring = ifelse(is.na(lipidwt), NA_character_, ""),
  exclude = FALSE,
  print_warning = TRUE
)
```

## Arguments

- conc:

  the concentration

- from:

  the current basis

- to:

  the target basis

- drywt:

  assumed to be a percentage taking values in \[0, 100\]

- lipidwt:

  assumed to be a percentage taking values in \[0, 100\]; assumed to be
  on a wet weight basis

- drywt_censoring:

  required when `drywt` is present, character: must be ”, '\<', 'D', 'Q'
  or NA

- lipidwt_censoring:

  required when `lipidwt` is present, character: must be ”, '\<', 'D',
  'Q' or NA

- exclude:

  a logical identifying records that do not need to be converted, useful
  when the data contain biological effects measurements

- print_warning:

  a boolean, gives the number of records lost during conversion
