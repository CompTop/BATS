CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3
PROFFLAGS = -pg -O3
DBFLAGS = -g -Og
CXX = g++
IFLAGS = -I../include

default: field.out dense_tests.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

.PHONY: clean
clean:
	rm *.out
