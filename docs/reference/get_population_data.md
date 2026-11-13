# get_population_data

Extracts population level data from a squid object

## Usage

``` r
get_population_data(x, list = FALSE, ...)
```

## Arguments

- x:

  an R object of class 'squid'

- list:

  Logical - whether to return data as a list or data.frame (FALSE;
  default).

- ...:

  further arguments passed to or from other methods.

## Value

A data.frame with the response(s), predictors and data structure for all
populations, or a list of data.frames, one for each population, if
list=TRUE.

## Examples

``` r
squid_data <- simulate_population(n=10, parameters=list(residual=list(vcov=1)))
get_population_data(squid_data)
#> Error in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 10, 0, 1
```
