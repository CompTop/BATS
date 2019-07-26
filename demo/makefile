CXXFLAGS = -Wall -g -std=c++17
CXX = g++
IFLAGS = -I../include

default: flag.out rips.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS)

clean:
	rm *.out
