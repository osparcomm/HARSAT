# Read a reference table

Utility function for reading reference tables:

## Usage

``` r
read_info_engine(id, files, compartment, file_diagnostics = FALSE)
```

## Arguments

- id:

  Character string identifying the reference table; currently must be
  one of `"determinand"`, `"species"`, `"thresholds"`, `"matrix"`,
  `"imposex"`,`"pivot_values", `"method_extraction"\`.

- files:

  Named list of file paths, with one of the names matching `id`

- compartment:

  Character string: must be one of `"biota"`, `"sediment"` or `"water"`

- file_diagnostics:

  Logical. The default (`FALSE`) is to print out a simple message
  confirming which file has been read. When `TRUE`, the MD5 digest and
  full path is also printed.

## Value

A reference table, or `NULL` if it is not required or cannot be found.

## Details

- checks whether the table is required

- checks whether the specified file exists

- calls the appropriate function to read the file

- generates messages, warnings and errors as appropriate
