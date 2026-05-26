README

MATLAB Codes for Peridynamic Dispersion Analyses

This repository contains MATLAB codes for peridynamic dispersion analyses related to the paper:

"Spatio-temporal dispersion analysis of dissipative time integration schemes for coupled peridynamic wave propagation"
by Minchul Yu and Gunwoo Noh

This repository will be made public after acceptance of the paper.

-------------------------------------------------------------------------------
1. Overview
-------------------------------------------------------------------------------

The codes in this repository are used to perform one- and two-dimensional
peridynamic dispersion analyses. The repository includes:

- one-dimensional peridynamic spatial and spatio-temporal dispersion analysis,
- two-dimensional peridynamic spatial and spatio-temporal dispersion analysis,
- verification code for the one-dimensional formulas presented in Appendix A,
- local stiffness matrix construction functions for one- and two-dimensional
  bond-based peridynamics,
- a utility function for identifying the discarded-mode boundary of the
  Noh-Bathe method.

The main purpose of these codes is to evaluate how peridynamic spatial
discretization and time integration schemes affect numerical wave propagation
behavior. In particular, the codes compute spatial dispersion relationships,
relative wave-speed errors, and numerical damping ratios for the central
difference method (CDM) and the Noh-Bathe (NB) method.

-------------------------------------------------------------------------------
2. File Description
-------------------------------------------------------------------------------

AppendixA_Formula_Verification.m

This script verifies the one-dimensional spatial and spatio-temporal dispersion
relationships presented in Appendix A of the paper. The user can select or
substitute explicit spatial dispersion formulas, such as FE or PD with different
horizon ratios, and combine them with CDM or NB spatio-temporal dispersion
formulas.

PD_1d_SpatioTemporal_dispersion.m

This script performs one-dimensional peridynamic spatial and spatio-temporal
dispersion analyses. It constructs the one-dimensional PD stiffness and mass
matrices, solves the spatial dispersion relationship between the numerical
wavenumber k and the exact wavenumber k0, and then computes the relative
wave-speed error and numerical damping ratio for CDM or the NB method.

PD_2d_SpatioTemporal_dispersion.m

This script performs two-dimensional peridynamic spatial and spatio-temporal
dispersion analyses. It constructs the two-dimensional PD stiffness and mass
matrices, evaluates P-wave and S-wave spatial dispersion curves for selected
propagation angles, and computes the corresponding relative wave-speed errors
and numerical damping ratios for CDM or the NB method.

Peri_const_mk_1d.m

This function returns the local stiffness contribution for one-dimensional
bond-based peridynamics.

Peri_const_mk_2d.m

This function returns the local stiffness contribution for two-dimensional
bond-based peridynamics. The output corresponds to the local contribution
associated with [u_i, u_j, v_i, v_j]. The function also includes the
distance-dependent influence function and partial-volume correction near the
horizon boundary.

NB_find_discard_kappa.m

This function estimates the normalized wavenumber at which the discarded-mode
criterion of the NB method is reached. The criterion is based on Delta t/T > 0.3
and is used as a practical guideline for identifying high-frequency modes that
are strongly attenuated by numerical dissipation.

-------------------------------------------------------------------------------
3. Notation
-------------------------------------------------------------------------------

The following notation is used consistently in the codes:

k       : numerical wavenumber
k0      : exact wavenumber
omega0  : exact angular frequency obtained from the discrete system
c       : numerical wave speed
c0      : exact wave speed
e       : relative wave-speed error, e = (c - c0)/c0
kappa   : normalized numerical wavenumber, kappa = k*Delta x
kappa0  : normalized exact wavenumber, kappa0 = k0*Delta x
CFL     : Courant-Friedrichs-Lewy number, CFL = c0*Delta t/Delta x
p       : splitting ratio of the NB method
m       : horizon ratio, m = delta/Delta x

In the generated figures, the horizontal axis is generally written as

    k*Delta x/pi

and the spatial dispersion result is written as

    k0*Delta x/pi.

-------------------------------------------------------------------------------
4. Main Parameters
-------------------------------------------------------------------------------

The main user-defined parameters are:

h              : nodal spacing, Delta x
m              : horizon ratio, delta/Delta x
kappa_range    : maximum value of k*Delta x
kappa_num      : number of wavenumber sampling points
CFL_set        : set of CFL numbers
p              : splitting ratio of the NB method
theta_set      : propagation angles for the two-dimensional analysis

Representative values used in the paper include:

CDM:
    CFL = 0.4, 0.7, 1.0

NB method:
    p = 0.54
    p = 2 - sqrt(2)

-------------------------------------------------------------------------------
5. How to Run
-------------------------------------------------------------------------------

1. Open MATLAB.
2. Add the repository folder to the MATLAB path.
3. Run the desired main script:
   - AppendixA_Formula_Verification.m
   - PD_1d_SpatioTemporal_dispersion.m
   - PD_2d_SpatioTemporal_dispersion.m
4. Modify parameters such as h, m, CFL_set, p, and theta_set if needed.
5. Check the generated figures for spatial dispersion, relative wave-speed
   error, and numerical damping ratio.

The helper functions Peri_const_mk_1d.m, Peri_const_mk_2d.m, and
NB_find_discard_kappa.m should be located in the same folder or included in
the MATLAB path.

-------------------------------------------------------------------------------
6. Important Notes
-------------------------------------------------------------------------------

- The codes are written for dispersion analysis and are not intended to be
  general-purpose peridynamic solvers.

- The variable k denotes the numerical wavenumber, while k0 denotes the exact
  wavenumber. This convention is used throughout the codes.

- The normalized wavenumbers are defined as:

      kappa  = k*Delta x
      kappa0 = k0*Delta x

- In the two-dimensional analysis, P-wave and S-wave branches are obtained
  from a reduced two-by-two eigenvalue problem based on Rayleigh quotient
  terms.

- The distance-dependent influence function used in the two-dimensional
  peridynamic stiffness calculation is w = 1 - |xi|/delta. Therefore, the
  m = 1 case is generally not recommended because nearest-neighbor bonds can
  lie on the horizon boundary, where the influence function becomes zero.
  For the dispersion analyses reported in the paper, m >= 2 is recommended.

- The discarded-mode boundary of the NB method is used as a practical
  diagnostic criterion based on Delta t/T > 0.3. It is not intended to define
  an exact analytical cut-off boundary.

-------------------------------------------------------------------------------
7. Requirements
-------------------------------------------------------------------------------

The codes were written for MATLAB R2023b. The Symbolic Math Toolbox may be required
for scripts that use symbolic variables or formula verification.

No external third-party MATLAB package is required.

-------------------------------------------------------------------------------
8. Citation
-------------------------------------------------------------------------------

If you use these codes, please cite the following paper after publication:

Minchul Yu and Gunwoo Noh,
"Spatio-temporal dispersion analysis of dissipative time integration schemes
for coupled peridynamic wave propagation."

-------------------------------------------------------------------------------
9. Contact
-------------------------------------------------------------------------------

For questions regarding the codes, please contact the corresponding author of
the paper.
