# Reads an input file, given a particular encoding, and possibly additional hints

Reads an input file, given a particular encoding, and possibly
additional hints

## Usage

``` r
safe_read_file(
  file,
  header = TRUE,
  sep = ",",
  quote = "\"",
  dec = ".",
  fill = TRUE,
  comment.char = "",
  strip.white = TRUE,
  ...
)
```

## Arguments

- file:

  the file name

- header:

  logical (default `TRUE`), whether or not a header line is expected

- sep:

  character (defailt `,`), field separator character

- quote:

  character (default `"`), field quote character

- dec:

  character (default `.`), decimal character

- fill:

  logical (default `TRUE`), if rows have unequal length, additional
  blank fields are added

- comment.char:

  character (default is empty), if not empty, a comment character

- strip.white:

  local (default `TRUE`), strips leading and trailing white s pace from
  fields

- ...:

  any additional parameters to
  [`utils::read.table`](https://rdrr.io/r/utils/read.table.html)

## Value

a data frame
