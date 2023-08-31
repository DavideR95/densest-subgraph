.DEFAULT_GOAL := main

CXX = icpx
CXXFLAGS = -Wall -pedantic -std=c++17 -g -O3
LDFLAGS =

all: clean main 

main: main.o graph.hpp Makefile main.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

main.o: main.cpp graph.hpp Makefile
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf *.o main
