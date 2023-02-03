CC = g++
DEBUG = -g
CFLAGS = -Wall -c ${DEBUG}
LFLAGS = -Wall ${DEBUG}
SRC = Physics.cpp User_IO.cpp ODE_Integrator.cpp Stellar_Structure_Equations.cpp Boundary_Conditions.cpp StatStar.cpp
OBJ = ${SRC:.cpp=.o}

%.o : %.cpp %.h
	${CC} -o $@ -c $<

statStar : ${OBJ}
	${CC} ${LFLAGS} ${OBJ} -o $@
	./$@

clean:
	@rm -rf *.o
