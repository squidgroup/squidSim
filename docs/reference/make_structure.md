# make_structure

Make nested and crossed balanced hierarchical structures

## Usage

``` r
make_structure(structure, repeat_obs = 1, level_names, int_names, ...)
```

## Arguments

- structure:

  A formula specifying the structure and sample sizes at each level. See
  details.

- repeat_obs:

  Number of repeated observations at the lowest level

- level_names:

  An optional list, containing names for the levels of different
  grouping factors

- int_names:

  An optional vector, containing names for the interactions (see
  details). These will be applied to the interaction terms in the same
  order they are specified in the structure argument

- ...:

  Further arguments passed to or from other methods.

## Value

A data.frame with the data structure

## Details

Factors are input as a text string. The name of each factor is followed
by the number of levels in that factor in brackets e.g.
"individual(100)". Nested factors can be specified using "/", e.g.
"population(2)/individual(2)", the lower levels being specified after
the higher levels, and the sample sizes of the lower levels. Crossed
factors are indicated using "+". Factor levels can be combined to make
new combinations, with unique levels for each combination, as an
'interaction' using ':', with factors already in the character existing
(e.g. "population(2)/individual(2) + population:individual").

## Author

Joel Pick - joel.l.pick@gmail.com

## Examples

``` r
# simple data structure with 5 'individuals' and 2 observations per individual
make_structure(structure="individual(5)", repeat_obs=2)
#>    individual
#> 1           1
#> 2           1
#> 3           2
#> 4           2
#> 5           3
#> 6           3
#> 7           4
#> 8           4
#> 9           5
#> 10          5

# nested data structure with 2 sexes, 5 individuals per sex
# and 2 observations per individual
make_structure(structure="sex(2)/individual(5)", repeat_obs=2)
#>    sex individual
#> 1    1          1
#> 2    1          1
#> 3    1          2
#> 4    1          2
#> 5    1          3
#> 6    1          3
#> 7    1          4
#> 8    1          4
#> 9    1          5
#> 10   1          5
#> 11   2          6
#> 12   2          6
#> 13   2          7
#> 14   2          7
#> 15   2          8
#> 16   2          8
#> 17   2          9
#> 18   2          9
#> 19   2         10
#> 20   2         10

# Naming the sexes in the data_structure
make_structure(structure="sex(2)/individual(5)", repeat_obs=2,
 level_names=list(sex=c("female","male")))
#>       sex individual
#> 1  female          1
#> 2  female          1
#> 3  female          2
#> 4  female          2
#> 5  female          3
#> 6  female          3
#> 7  female          4
#> 8  female          4
#> 9  female          5
#> 10 female          5
#> 11   male          6
#> 12   male          6
#> 13   male          7
#> 14   male          7
#> 15   male          8
#> 16   male          8
#> 17   male          9
#> 18   male          9
#> 19   male         10
#> 20   male         10

# crossed data structure with 5 individuals in 2 treatments 
# and 2 observations per individual and treatment combination
make_structure(structure="treatment(2) + individual(5)", repeat_obs=1)
#>    treatment individual
#> 1          1          1
#> 2          1          2
#> 3          1          3
#> 4          1          4
#> 5          1          5
#> 6          2          1
#> 7          2          2
#> 8          2          3
#> 9          2          4
#> 10         2          5
```
