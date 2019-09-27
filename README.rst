partdens
========
Partitioning schemes for molecular electron densities.

Created 2018 -- Ryan J. MacDonell and Michael S. Schuurman (uOttawa, NRC)

Scripts
-------
bader.f90
  Partial atomic charges calculated using a grid-based implementation of
  Bader's Quantum Theory of Atoms in Molecules method.

gridchg.f90
  Partial atomic charges calculated on a Becke grid with options for
  Becke, Voronoi and Hirshfeld weighting function, with or without
  iterative deformation densities.

nci.f90
  Non-covalent interaction analysis on a cartesian grid (work in progress).

nat_ibo.f90
  Natural orbital extension of Knizia's Intrinsic Atomic Orbitals and
  Intrinsic Bond Orbitals (work in progress).

newton.f90
  Partial atomic charges calculated by fitting atomic densities to the
  molecular density using the Newton-Raphson method (work in progress).

Libraries
---------
Fortran libraries by Serguei Patchkovskii (MBI), Michael Spanner (NRC)
and Sergei N. Yurchenko (UCL):

- accuracy.f90
- atoms.f90
- dgedi.f
- dgefa.f
- gamess_internal.f90
- import_gamess*.f90
- lapack.f90
- lebedev.f90
- math.f90
- molecular_grid.f90
- os_integral_operators.f90
- timer.f90

Fortran library by Simon Neville (NRC):

- molden.f90
