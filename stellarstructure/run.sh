mkdir -p data
c++ -c SSFunctions.cpp
c++ -c Initialize.cpp
c++ -c main.cpp
c++ main.o SSFunctions.o Initialize.o
