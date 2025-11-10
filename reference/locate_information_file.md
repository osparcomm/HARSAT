# Find the path for an information file

Locates a requested information file, searching the information file
path. If the requested file cannot be found and `required` is not false,
stops with an error.

## Usage

``` r
locate_information_file(name, path)
```

## Arguments

- name:

  A string: the name of the file, e.g., `thresholds_biota.csv`

- path:

  A vector of strings, directories to search. The information directory
  for the package is automatically searched if we haven't found the file
  anywhere else

## Value

A string, the absolute path for the file, or `NULL` if the file cannot
be found anywhere.
