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
#>              y    residual squid_pop
#> 1   1.58404814  1.58404814         1
#> 2   0.11770964  0.11770964         1
#> 3  -0.21549876 -0.21549876         1
#> 4   0.11059731  0.11059731         1
#> 5   0.01965524  0.01965524         1
#> 6  -1.49344407 -1.49344407         1
#> 7   0.72746304  0.72746304         1
#> 8   1.14817156  1.14817156         1
#> 9   0.30232378  0.30232378         1
#> 10  0.69576900  0.69576900         1
```
