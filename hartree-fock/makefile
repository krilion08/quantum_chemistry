# Compiler
CXX = g++
# Compiler flags
CXXFLAGS = -I. -I./Eigen -std=c++11 -w
# Executable name
EXE = hf.x
# Source files
SRCS = main.cc molecule.cc masses.cc
# Object files
OBJS = $(SRCS:.cc=.o)

# Default target
all: $(EXE)

# Rule to link the executable
$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Rule to compile source files into object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(EXE) $(OBJS)

# Phony targets
.PHONY: all clean test

# Rule to run the test
test: $(EXE)
	./$(EXE) ./examples/sto-3g_water/
	rm -f $(OBJS)