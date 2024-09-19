# Definición de variables para simplificar cambios futuros
CXX=g++
CXXFLAGS=-O3
SOURCES = main.cpp Lattice_dynamics.cpp Lattice_Init.cpp vector.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = Hartmann.x

all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ 


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(OBJECTS) $(EXECUTABLE) *.dat 