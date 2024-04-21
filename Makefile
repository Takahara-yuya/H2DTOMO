CXX = mpicxx
CXXFLAGS = -O3
LDFLAGS = -L/opt/mkl/mkl/lib/intel64
LDLIBS = -lmkl_rt
INCLUDES = -I/opt/mkl/mkl/include

SRC_DIR = src
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = H2DTOMO

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
