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
#>             y individual_effect   residual individual squid_pop
#> 1   0.5438577         0.4860238  0.0578339          1         1
#> 2  -0.5256584         0.4860238 -1.0116822          1         1
#> 3  -1.2242032        -0.9495541 -0.2746491          2         1
#> 4  -0.3345312        -0.9495541  0.6150229          2         1
#> 7   1.1681690         0.3960297  0.7721393          4         1
#> 8   0.9287238         0.3960297  0.5326940          4         1
#> 13 -0.8620939        -0.5659450 -0.2961489          7         1
#> 14 -0.4155502        -0.5659450  0.1503948          7         1
#> 19 -2.3305708        -1.3454021 -0.9851688         10         1
#> 20 -1.7387584        -1.3454021 -0.3933563         10         1
```
