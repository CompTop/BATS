CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3 -fopenmp
PROFFLAGS = -pg -O3
DBFLAGS = -g -Og
CXX = g++
IFLAGS = -I../include

default: rips.out serpinski.out zigzag_codensity.out zigzag_cover.out zigzag_rips_samples.out zigzag_witness.out rips_filtration.out witness_filtration.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

.PHONY: clean
clean:
	rm -f *.out *.prof *.db
	rm -rf *.dgm
