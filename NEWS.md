
squidSim 0.1.0
===========

-   First release.



squidSim 0.2.0
===========

-   Simulation of additive genetic effects now uses MCMCglmm::rbv internally. It is substantially faster!
-   Removed nadiv dependency until it back on CRAN, so temporarily doesn't simulate dominance or epigenetic effects internally. This can still be done by including matrices generated from nadiv as covariance structures
