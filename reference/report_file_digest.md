# Generates and logs a file digest

Writes a line on the console containing the file name and its MD5
digest. This can be used to track changes in input files, for
reproducibility reasons. MD5 is good enough for this as we aren't after
cryptograohic levels of security, and its cheaper to compute.

## Usage

``` r
report_file_digest(infile)
```

## Arguments

- infile:

  the input file name
