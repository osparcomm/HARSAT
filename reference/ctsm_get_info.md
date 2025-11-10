# Gets data from the standard reference tables

Gets data from the standard reference tables

## Usage

``` r
ctsm_get_info(
  ref_table,
  input,
  output,
  compartment = NULL,
  na_action = c("fail", "input_ok", "output_ok", "ok"),
  info_type = NULL,
  sep = "."
)
```

## Arguments

- ref_table:

  the reference table

- input:

  the input

- output:

  the output

- compartment:

  the compartment

- na_action:

  character, how to handle missing values, one of: `fail` (no missing
  values allowed), `input_ok` (allows missing values in input, but all
  non-missing values must be recognised, and must have output),
  ``` output_ok`` requires all input values to be present, but allows missing values in output (for e.g. dryweight by species), or  ```ok\`
  (allows missing values everywhere)

- info_type:

  the information type, by default picked up from the function call

- sep:

  character, the separator, defaulting to '.'
