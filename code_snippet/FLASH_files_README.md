Contained in this directory are files that use the ADM quantities to evolve the hydrodynamics.

In `hy_rk_updateSoln.F90` look for !!~~ symbols in the code that will highlight areas of interest for you.  Specifically,
* `genrelSources()` subroutine calculates the source terms according to Eqn. (17) of Mosta+ 2013 (https://arxiv.org/pdf/1304.5544.pdf)
* `calc_metric_derivatives()` subroutine calculates the (relevant) Christoffel symbols and metric derivatives according to Table A.1 of O'Connor & Ott 2011 (https://arxiv.org/pdf/0912.2393.pdf) 
  * This is where we will need spatial and time derivatives of metric quantities
* The bottom of the `prim2conGR.F90` file calculates the stress energy tensor
* `evolve_spacetime.F90` evolves the metric quantities (radial gauge polar slicing coordinates only so far)
  * This will eventually be moved externally and replaced by the spacetime solver
* At the bottom of `hy_rk_updateSoln()`, the metric potential is porperly shifted (post integration), finally calculating the lapse and g00 covariant metric quantities
