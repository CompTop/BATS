CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3 -fopenmp -lomp
PROFFLAGS = -pg -O3
DBFLAGS = -g -O0
OPFLAGS = -O3 -DBATS_OPCOUNT
CXX ?= g++
IFLAGS = -I../include

default: rips.out serpinski.out zigzag_codensity.out zigzag_cover.out zigzag_rips_samples.out zigzag_witness.out rips_filtration.out witness_filtration.out

zigzag: zigzag_rips_samples.out

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
