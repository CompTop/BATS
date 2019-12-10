CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3 -fopenmp
PROFFLAGS = -pg -O3
DBFLAGS = -g -O0 
CXX = g++
IFLAGS = -I../include

default: rips.out rips.prof rips.db

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

.PHONY: clean
clean:
	rm *.out *.prof *.db
