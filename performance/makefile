CPPFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3 -fopenmp
PROFFLAGS = -pg -O3
DBFLAGS = -g -Og
OPFLAGS = -O3 -DBATS_OPCOUNT
CPP ?= g++
IFLAGS = -I../include

default: divide_conquer.out

%.out: %.cpp
	$(CPP) -o $@ $< $(IFLAGS) $(CPPFLAGS) $(OUTFLAGS)

%.op: %.cpp
	$(CPP) -o $@ $< $(IFLAGS) $(CPPFLAGS) $(OPFLAGS)

%.prof: %.cpp
	$(CPP) -o $@ $< $(IFLAGS) $(CPPFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CPP) -o $@ $< $(IFLAGS) $(CPPFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

.PHONY: clean
clean:
	rm -f *.out *.prof *.db
	rm -rf *.dgm
