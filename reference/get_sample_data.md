# get_sample_data

Extracts sampled data from a squid object

## Usage

``` r
get_sample_data(x, sample_set = 1, list = FALSE, ...)
```

## Arguments

- x:

  an R object of class 'squid'

- sample_set:

  Integer - which sample set to return. Defaults to 1

- list:

  Logical - whether to return data as a list or data.frame (FALSE;
  default).

- ...:

  further arguments passed to or from other methods.

## Value

A data.frame with the sampled data for all populations, or a list of
data.frames, one for each population, if list=TRUE.

## Examples

``` r
squid_data <- simulate_population(
  data_structure = make_structure("individual(10)", repeat_obs=2),
  parameters = list(individual=list(vcov=1), residual=list(vcov=1)),
  sample_type = "nested",
  sample_param = matrix(5, dimnames=list(NULL,"individual"))
)
get_sample_data(squid_data)
#>             y individual_effect    residual individual squid_pop
#> 5  -2.7736273        -2.2147291 -0.55889829          3         1
#> 6  -0.3217541        -2.2147291  1.89297496          3         1
#> 9  -0.2039293         0.6986377 -0.90256692          5         1
#> 10  0.3165975         0.6986377 -0.38204014          5         1
#> 15  0.3074246        -0.3911199  0.69854448          8         1
#> 16 -2.4287862        -0.3911199 -2.03766634          8         1
#> 17 -0.2629960        -0.4343846  0.17138861          9         1
#> 18  0.4423013        -0.4343846  0.87668584          9         1
#> 19 -0.6111450        -0.6334419  0.02229697         10         1
#> 20  0.6658022        -0.6334419  1.29924414         10         1
```
