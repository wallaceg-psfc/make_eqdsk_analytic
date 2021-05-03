# make_eqdsk_analytic
Analytic EQDSK generator

MAKE_ANALYTIC_EQDSK(eq_option_in,R0_in,eps_in,delta_X_in,kappa_X_in,delta_in,kappa_in,nu_in,B0)
generates an analytic GEQDSK equilibrium file based on the method of Guazzotto and Freidberg (JPP 2021).  

Input arguments:
eq_option_in = 1 for elongated limiter equilibrium, 2 for double null, 3 for single null (lower)
R0_in = major radius of plasma center
eps_in = inverse aspect ratio
delta_X_in = triangularity of plasma near x-point (eq_options 2 and 3)
kappa_X_in = elongation of plasma near x-point (eq_options 2 and 3)
delta_in = triangularity of plasma away from x-point (eq_options 1, 2, and 3)
kappa_in = elongation of plasma away from x-point (eq_options 1, 2, and 3)
nu_in = beta_p like parameter (see JOPP 2021 paper for details)
B0 = magnetic field on axis in Tesla

The script generates a geqdsk output file with the name 'geqdsk_analytic' upon completion.

The script can be buggy if input parameters go outside typical bounds for tokamak equilibria.  Users should experiment with input parameters to get a feel for the code before using in an automated manner.

The heart of the code was originally written by Luca Guazzotto (Auburn U.), and adapted for the GEQDSK format by Greg Wallace (MIT).

Requires a recent release of MATLAB.  Tested with R2020b by the author.  Earlier versions (e.g. R2015b) may not work.
