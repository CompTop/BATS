CXXFLAGS = -Wall -std=c++17 -Wextra
OUTFLAGS = -O3
PROFFLAGS = -pg -O3
DBFLAGS = -g -Og
CXX ?= g++
IFLAGS = -I../include

default: field.out dense_tests.out col_matrix.out sparse_vector.out util.out reduction.out quiver.out

%.out: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(OUTFLAGS)

%.prof: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(PROFFLAGS)

%.db: %.cpp
	$(CXX) -o $@ $< $(IFLAGS) $(CXXFLAGS) $(DBFLAGS)

%: %.out %.prof %.db

test: field.out dense_tests.out col_matrix.out sparse_vector.out util.out reduction.out quiver.out
	./field.out
	./dense_tests.out
	./col_matrix.out
	./sparse_vector.out
	./util.out
	./reduction.out
	./quiver.out

.PHONY: clean
clean:
	rm *.out
