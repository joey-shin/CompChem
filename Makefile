# Define the C++ compiler
CXX = g++

# Compilation flags: -Wall for compiler warnings, -g for debugging information
CXXFLAGS = -Wall -std=c++14

# Include directories for header files
INCLUDES = -I/opt/homebrew/include -I/opt/homebrew/opt/openblas/include

# Library paths, if Armadillo is not in a standard location
LFLAGS = -L/opt/homebrew/lib -L/opt/homebrew/opt/openblas/lib

# Libraries to link against, including Armadillo
LIBS = -larmadillo -lopenblas

# Source files
SRCS = main.cpp input/input.cpp basis/basis.cpp

# Object files, derived from SRCS
OBJS = $(SRCS:.cpp=.o)

# Executable name
MAIN = chem

# Targets
.PHONY: depend clean

all: $(MAIN)
	@echo EXECUTABLE COMPILED

$(MAIN): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# Rule to convert .cpp to .o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^


