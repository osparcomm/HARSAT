# Find the path for an information file

Locates a requested information file, searching the information file
path, and then looking in the package itself.

## Usage

``` r
locate_information_file(path, name)
```

## Arguments

- path:

  A string giving the directory to search. The information directory for
  the package is automatically searched if we haven't found the file in
  this directory

- name:

  A string: the name of the file, e.g., `thresholds_biota.csv`

## Value

A string, the path to the file, or `""` if the file cannot be found.
