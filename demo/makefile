CXXFLAGS = -Wall -g -std=c++17
CXX = g++
IFLAGS = -I../include

default: simplicial_complex.out field.out vector.out filtration.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS)

clean:
	rm *.out
