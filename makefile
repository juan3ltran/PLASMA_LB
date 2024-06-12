# Definición de variables para simplificar cambios futuros
CXX=g++
CXXFLAGS=-O3
dbFLAGS= -g -Wall -fsanitize=address -fsanitize=undefined -fsanitize=leak
simulation_SOURCE=MHDM.cpp
simulation_EXECUTABLE=MHDM.x

.PHONY: all simulation  clean

# Objetivo por defecto que construye la simulacion
all: simulation 

# Objetivo para construir y ejecutar la parte de simulation
simulation: $(simulation_SOURCE)	vector.o
	$(CXX) $(CXXFLAGS) $(simulation_SOURCE) vector.o -o $(simulation_EXECUTABLE)

debug: 
	$(CXX) $(dbFLAGS) $(simulation_SOURCE) vector.o -o $(simulation_EXECUTABLE)

vector.o: vector.cpp
	g++ -c vector.cpp -o vector.o

# Objetivo para limpiar los archivos generados
clean:
	rm -f *.o *.x *.out