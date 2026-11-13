
squidSim 0.1.0
===========

-   First release.


squidSim 0.2.0
===========

-   Simulation of additive genetic effects now uses `MCMCglmm::rbv` internally. It is substantially faster!
-   Removed nadiv dependency until it back on CRAN, so temporarily doesn't simulate dominance or epigenetic effects internally. This can still be done by including matrices generated from nadiv as covariance structures


squidSim 0.2.1
===========

-   Bug fixes with simulating additive genetic effects, and added index_link argument to help with indexing in the model argument. This allows new factors in the data structure to be made which are indexed by other factors. Previously simulating maternal genetic effects didn't work properly, because of the indexing, and this solves that problem.


squidSim 0.2.2
===========

-   New survival sampling functionality. If you simulate binomial data with a age structure, then you can sample so that you only have the ages up to (and including) the first 0 or 1 (depending on whether you are simulating survival or mortality)


squidSim 0.2.3
===========

- Bug with 'fixed' factor is fixed. Now if the levels in the data structure have names, the same names have to be specified in the parameter list. The order doesn't matter (this will be matched internally), but the names have to match exactly. If the factor levels in the data structure are numeric then these do not need to match; the order of the effects will match the order of the numeric factor levels.


squidSim 0.2.4
===========
- Added 'int_names' argument to `make_structure` function, that allows combination/interaction terms to be named. 


squidSim 0.2.5
===========
- Added cloglog link function to `simulate_population`
- Added functions `exp2lat` and `lat2exp`, which convert means and (co)variances between normal and log normal distributions, and so can be used to transform between latent and expected scales with GLMs
- Fixed phylogenetic effect simulation


squidSim 0.2.6
===========
- Optimising sparse matrices to speed up simulating with `cov_str`
- Updating internal error checking


squidSim 0.3.0
===========
- This version follows a comprehensive bug and documentation check using Claude Opus 5.5 medium. All bugs found and changes suggested by Claude were checked thoroughly by maintainers, and no automated bug correction was used.
- Without a `seed`, repeated calls to `simulate_population` no longer cycle through a small set of datasets. `simulate_population` now leaves the user's random number stream as it found it
- `link="sqrt"` now applies the inverse link (squares the linear predictor)
- Predictors with 0 variance now equal their mean, rather than 0
- `mean` is now used with `pedigree` and `phylogeny`, and `mean` is no longer transformed by `cov_str`
- `covariate=TRUE` now uses the values in the data structure, rather than their ranks
- Fixed factors can have NAs in the data structure
- `cov_str` can have more levels than the data structure, and is checked for being symmetric and positive definite
- Singular `vcov` matrices can now be simulated
- Singular `vcov` matrices and pedigrees with more than 3 columns now work with `pedigree`
- `get_parameters` works with interactions and includes known predictors
- Names on the intercept and column names on betas are now ignored
- tibbles and data.tables can be used as the data structure
- Interactions now work with predictors with numeric names
- Sampling fixes: Nested sampling works with one grouping factor, or a group with a single level and with factors, and without a data structure. Now gives an error if more levels are requested than any group has, and a warning if only some groups are short. Temporal sampling works with non-consecutive and named time points, when variance is 0, and no longer excludes first time point. Missing sampling work when groups end up with equal sized samples.
- `simulated_variance` now includes covariances involving two way interactions in the total variance
- More informative errors in `make_structure` and `get_sample_data`
- Added tests. Now requires R (>= 4.1.0)
- `simulated_variance` removed from this version. Is now in the development branch, and will be reinstated in a later version
