mkdir -p data
c++ -c -g main.cpp
c++ -c -o SSFunctions.o SSFunctions.cpp
c++ -c -g Initialize.cpp
c++ -g main.o SSFunctions.o Initialize.o -o stellarStructure
./stellarStructure
