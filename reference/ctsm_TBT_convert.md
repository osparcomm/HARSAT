# Convert tin concentrations

Convert tin concentrations to cation concentrations.

## Usage

``` r
ctsm_TBT_convert(
  data,
  subset,
  action,
  from = c("tin", "cation"),
  convert_var = c("value", "limit_detection", "limit_quantification", "uncertainty")
)
```

## Arguments

- data:

  the data

- subset:

  optional, a subset expression

- action:

  one of `"relabel"`, `"convert"`, `"change_unit"` – note that
  `change_unit` moves units from tin units to conventional units

- from:

  either `"tin"` or `"cation"`

- convert_var:

  one of `"value"`, `"limit_detection"`, `"limit_quantification"`,
  `"uncertainty"`
