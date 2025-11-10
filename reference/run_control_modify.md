# Modifies control parameters for `run_assessment`

Undates default control parameters with user specification and does
basis error checking.

## Usage

``` r
run_control_modify(control_default, control = list())
```

## Arguments

- control:

  List of replacement control parameters; defaults to an empty list.

- run_control_default:

  List of default control parameters produced by a call to
  `run_control_default`

## Value

List of updated control parameters
