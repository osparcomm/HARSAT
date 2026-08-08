# extract determinands from control structure

gets the names of all determinands involved in argument
`determinands.control`

## Usage

``` r
get_control_dets(control, .names = TRUE)
```

## Arguments

- control:

  list of control structures passed into `create_timeseries` by argument
  `determinand.control`

- .names:

  logical determining whether the names of the list are included (`TRUE`
  default) or not (`FALSE`)

## Value

charcater string of determinands
