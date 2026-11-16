# simulate_population

Simulate population level data

## Usage

``` r
simulate_population(
  data_structure,
  n,
  parameters,
  n_response = 1,
  response_names,
  known_predictors,
  model,
  index_link,
  family = "gaussian",
  link = "identity",
  pedigree,
  phylogeny,
  cov_str,
  sample_type,
  sample_param,
  n_pop = 1,
  seed,
  verbose = FALSE,
  suppress_index_warning = FALSE
)
```

## Arguments

- data_structure:

  A matrix or data.frame with a named column for each grouping factor,
  including the levels

- n:

  Sample size when data_structure is not specified

- parameters:

  A list of parameters for each hierarchical level. See details.

- n_response:

  The number of response variables, defaults to 1.

- response_names:

  Names given to response variables. Defaults to 'y', or
  c('y1','y2',...) if n_response\>1. Not used if model is specified.

- known_predictors:

  This argument provides a way of inputting existing predictor
  variables. This argument takes a list, with items 'predictors' and
  'beta', where 'predictors' is a matrix or data.frame, the number of
  rows of which MUST equal 'n' or the number of rows in
  'data_structure'. The column names of 'predictors' are used as the
  variable names, and the 'beta's are assumed to be in the same order as
  the predictors.

- model:

  Optional. A equation description of the simulation model given as a
  character string. See details.

- index_link:

  Optional. Make new factors in the data structure indexed by other
  factors. This is only really used in the context of the model
  argument. Takes a list, the names of which indicate the new grouping
  factor name and the elements of the list represent which grouping
  factors should be linked, for example "mother-ID" would index mother
  with the same indexes as ID.

- family:

  A description of the error distribution. Options are 'gaussian'
  (default), 'poisson' and 'binomial'. 'binomial' generates a binary
  response variable.

- link:

  A description of the link function distribution. Options are
  'identity' (default),'log', 'inverse', 'sqrt', 'logit', 'probit' and
  'cloglog'.

- pedigree:

  A list of pedigrees for each hierarchical level. Each pedigree must be
  matrix or data.frame, that is at least 3 columns, which correspond to
  ID, dam and sire. The name in the pedigree list must match a name in
  the parameter list.

- phylogeny:

  A list of phylogenies for each hierarchical level. Each phylogeny
  should be phylo class. The name in the phylogeny list must match a
  name in the parameter list.

- cov_str:

  A list of covariance structures for each hierarchical level. The name
  in the cov_str list must match a name in the parameter list.

- sample_type:

  Type of sampling, must be one of 'nested', 'missing', 'survival' or
  'temporal.' If not specified, then no sampling is done. See details

- sample_param:

  A set of parameters, specific to the sampling type. See details.

- n_pop:

  Number of populations. Default = 1

- seed:

  Seed for random number generation. Setting seed allows for simulations
  to be easily reproduced. If seed is not set, one is set internally,
  which is then saved and returned in the squidSim object

- verbose:

  Logical. Whether to print diagnostics. Useful for debugging. Defaults
  to FALSE

- suppress_index_warning:

  Logical. Whether to print warnings relating to the index-link
  argument. Useful to switch off if using in large number of
  simulations. Defaults to FALSE.

## Value

a squid object, which is a list including all inputs and simulated data.

## Details

A detailed vignette can be found at
https://squidgroup.org/squidSim_vignette/

The parameters list contains one (or more) list for each hierarchical
level that you want to simulate at. A residual list is always needed,
specifying variances/covariances for the residual. Additionally, the
parameter list can also be provided with an intercept vector and
interactions list. For each item in the parameter list (excluding
intercept, interactions, and residual), the following can be specified
(but all have default values): names - vector containing the names of
predictors simulated at this level group - character string relating to
the data_structure mean - vector of means for the predictor variables
vcov or vcorr - Either a vector of variances, or a
variance-covariance/correlation matrix, for the predictor variables
beta - vector of effect sizes (or matrix with n_response columns when
n_response\>1) fixed - Logical, indicating whether the effects for the
levels are fixed or to be simulated covariate - Logical, indicating
whether the indexes in the data structure are to be used as a continuous
variable rather than simulating one functions - vector - transformation
to be applied to the response variable. Defaults to ‘identity’. A more
detailed explanation can be found at
https://squidgroup.org/squidSim_vignette/1.9-parameter-list-summary.html

The model argument is character string which explicitly tells the
simulate_population function how to put the simulated predictors
together to form the response variable. For example, if the predictors
temperature and rainfall had been specified in the parameter list,
providing the model argument with "y = temperature + rainfall +
residual" would result in generating a response variable 'y' in the same
way the simulate_population function does by default. For more detailed
information see http://squidgroup.org/squidSim_vignette/1.7-modeleq.html

Different sampling schemes can be implemented, (sample_type can be
'nested', 'missing', 'survival' or 'temporal'). The sample_param takes a
different form depending on the sample_type. See
https://squidgroup.org/squidSim_vignette/7-sampling.html for full
details.

## Author

Joel Pick - joel.l.pick@gmail.com

## Examples

``` r
# simple linear model with three predictors variables
squid_data <- simulate_population(
  n=50,
  parameters = list(
    observation = list(
      names = c("temperature","rainfall", "wind"),
      beta = c(0.5,-0.3, 0.4)    
  ),
    residual = list(
      vcov = 1
    )
  )
)
```
