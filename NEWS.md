
squidSim 0.1.0
===========

-   First release.



squidSim 0.2.0
===========

-   Simulation of additive genetic effects now uses MCMCglmm::rbv internally. It is substantially faster!
-   Removed nadiv dependency until it back on CRAN, so temporarily doesn't simulate dominance or epigenetic effects internally. This can still be done by including matrices generated from nadiv as covariance structures



squidSim 0.2.1
===========

-   Bug fixes with simulating additive genetic effects, and added index_link argument to help with indexing in the model argument. This allows new factors in the data structure to be made which are indexed by other factors. Previously simulating maternal genetic effects didn't work properly, because of the indexing, and this solves that problem.



squidSim 0.2.2
===========

-   New survival sampling functionality. If you simulate binomial data with a age structure, then you can sample so that you only have the ages up to (and including) the first 0 or 1 (depending on whether you are simulating survival or mortality)