CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3 -fopenmp
PROFFLAGS = -pg -O3
DBFLAGS = -g -Og
OPFLAGS = -O3 -DBATS_OPCOUNT
CXX ?= g++
IFLAGS = -I../include

default: divide_conquer.out zigzag_circle.out zigzag_cube.out scaling.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.op: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OPFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

.PHONY: clean
clean:
	rm -f *.out *.prof *.db
	rm -rf *.dgm
