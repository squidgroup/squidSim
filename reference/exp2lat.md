# exp2lat

Transforms means and covariance matrix from expected (log-normal) to
latent (normal) scale

## Usage

``` r
exp2lat(mean, cov)
```

## Arguments

- mean:

  a vector of means. These need to all be \>0.

- cov:

  a single variance or a covariance matrix

## Value

a list, which an element 'mean' for the transformed mean(s) and another
element 'cov' for the transformed covariance matrix

## Details

Transforms expected mean(s) and (co)variance(s) of a multivariate
log-normal distribution to expected mean(s) and (co)variance(s) of a
multivariate normal distribution. This is useful when trying to go from
expected (log-normal) to latent (normal) scales in a Poisson GLM with a
log link function.

## Author

Joel Pick - joel.l.pick@gmail.com

## Examples

``` r
exp2lat(1,1)
#> $mean
#> [1] -0.3465736
#> 
#> $cov
#>           [,1]
#> [1,] 0.6931472
#> 
```
