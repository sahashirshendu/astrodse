g++ -c Physics.cpp User_IO.cpp ODE_Integrator.cpp Stellar_Structure_Equations.cpp Boundary_Conditions.cpp StatStar.cpp
g++ -o starmodl Physics.o User_IO.o ODE_Integrator.o Stellar_Structure_Equations.o Boundary_Conditions.o StatStar.o
./starmodl
