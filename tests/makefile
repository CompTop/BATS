CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3
PROFFLAGS = -pg -O3
DBFLAGS = -g -Og
CXX ?= g++
IFLAGS = -I../include

ifeq ($(CXX),g++)
  OUTFLAGS += -fopenmp
else ifeq ($(CXX),clang++)
  OUTFLAGS += -lomp
else
	OUTFLAGS += -lomp -DMESSAGE='"Compiled with an unknown compiler"'
endif


default: field.out col_matrix.out sparse_vector.out util.out reduction.out quiver.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

test: field.out col_matrix.out sparse_vector.out util.out reduction.out quiver.out filtration.out
	./field.out
	./col_matrix.out
	./sparse_vector.out
	./util.out
	./reduction.out
	./quiver.out
	./filtration.out

.PHONY: clean
clean:
	rm *.out
