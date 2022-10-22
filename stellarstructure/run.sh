mkdir -p data
g++ -c SSFunctions.cpp
g++ -c main.cpp
g++ -c Initialize.cpp
g++ main.o SSFunctions.o Initialize.o
