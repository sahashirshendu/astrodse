This folder contains the modules and main routine for StatStar, written in Fortran 95.  The modules must be compiled first, before they can be linked to the main routine.  You will also need Constants.F90 which is available in the Appendix I folder.

The list of required files to compile and link StatStar are:

StatStar.f90  (Main Driver routine)
Boundary Conditions.f90
Composition.f90
ODE Integrator.f90
Physics.f90
Stellar Structure Equations.f90
User IO.f90
Zone Quantities.f90

and from Appendix I:

Constants.f90




StatStar is described in Chapter 10 and Appendix L of:

An Introduction to Modern Astrophysics
Bradley W. Carroll and Dale A. Ostlie
Addison Wesley
Copyright 2007.

