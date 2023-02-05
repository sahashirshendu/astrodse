gfortran -c Constants.f90 Composition.f90 Zone_Quantities.f90 Physics.f90 User_IO.f90 ODE_Integrator.f90 Stellar_Structure_Equations.f90 Boundary_Conditions.f90 StatStar.f90
gfortran Constants.o Composition.o Zone_Quantities.o Physics.o User_IO.o ODE_Integrator.o Stellar_Structure_Equations.o Boundary_Conditions.o StatStar.o
./a.out
