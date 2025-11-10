# Detects the environment from the call

This is a utility function that detects the package environment. If you
have imported `harsat` as a package, it returns the package environment.
Otherwise, it returns the global environment. You can safely export
functions from the result of this, for example, when you are setting up
a cluster of child processes for parallel computation.

## Usage

``` r
cstm.VDS.environment()
```
