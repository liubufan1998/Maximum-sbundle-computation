# ********   Author: Yang Liu    ******
# ******** Email:  ******
#
CXX = g++
CXXFLAGS = -std=c++11 -flto -O3
INCLUDES = -I.
SOURCES = main.cpp Graph.cpp sbundle_tool.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = MSBP_SymBD_H

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clean-all:
	rm -rf CMakeFiles cmake_install.cmake CMakeCache.txt Makefile *.cmake
	find . -type f -name '*.o' -delete
	find . -type f -name '$(EXECUTABLE)' -delete

.PHONY: all clean clean-all
