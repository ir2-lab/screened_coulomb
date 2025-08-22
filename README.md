# screened_coulomb

Header-only C++ library for Screened Coulomb scattering calculations.

### Features

- All code in a single header file.
- Different types of screened potentials implemented
  - Ziegler-Biersack-Littmark (ZBL) Universal potential
  - Bohr
  - Kr-C
  - Moliere
- Classical calculation of scattering angle, scattering cross-section, stopping cross-section in both center-of-mass and lab systems
- Different types of numerical evaluation of scattering integrals
  - Gauss-Chebyshev quadrature
  - 4-th order Lobatto quadrature
  - Approximate analytic MAGIC formula
  - Exact results for the un-screened potential

### Installation

Just copy and include `include/screened_coulomb.h` in your project.

### Usage

`test/test.cpp` contains some tests and example code.



