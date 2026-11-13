# lat2exp

Transforms means and covariance matrix from latent (normal) to expected
(log-normal) scale

## Usage

``` r
lat2exp(mean, cov)
```

## Arguments

- mean:

  a vector of means

- cov:

  a single variance or a covariance matrix

## Value

a list, which an element 'mean' for the transformed mean(s) and another
element 'cov' for the transformed covariance matrix

## Details

Transforms expected mean(s) and (co)variance(s) of a multivariate normal
distribution to expected mean(s) and (co)variance(s) of a multivariate
log-normal distribution. This is useful when trying to go from latent
(normal) to expected (log-normal) scales in a Poisson GLM with a log
link function.

## Author

Joel Pick - joel.l.pick@gmail.com

## Examples

``` r
lat2exp(1,1)
#> $mean
#> [1] 4.481689
#> 
#> $cov
#>          [,1]
#> [1,] 34.51261
#> 
```
