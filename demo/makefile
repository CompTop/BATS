CXXFLAGS = -Wall -g -std=c++17 -O3
TESTFLAGS = -pg
CXX = g++
IFLAGS = -I../include

default: rips.test rips.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS)

%.test: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(TESTFLAGS)

clean:
	rm *.out *.test
