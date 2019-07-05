CXXFLAGS = -Wall -g
CXX = g++
IFLAGS = -I../include

default: simplicial_complex.out field.out vector.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS)

clean:
	rm *.out
