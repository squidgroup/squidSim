# get_parameters

Extracts the simulation parameters from a squid object

## Usage

``` r
get_parameters(x)
```

## Arguments

- x:

  an R object of class 'squid'

## Value

A named vector of the means, variances, covariances and betas used in
the simulation.

## Examples

``` r
squid_data <- simulate_population(n=10, parameters=list(residual=list(vcov=1)))
get_parameters(squid_data)
#>     intercept residual_mean  residual_var residual_beta 
#>             0             0             1             1 
```
