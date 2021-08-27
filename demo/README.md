# Demo Files

This folder contains a variety of files which demonstrate various capabilities of BATS.

**Warning** Some may be out of data and may not compile.

## Basic compilation
You can make a standard executable compiled with OpenMP and GCC optimizations by specifying the `*.out` extension.  For example
```
make  simplicial_complex.out
```

## Debugging

You can compile an executable for debugging with `gdb` by specifying a `*.db` extension
```
make simplicial_complex.db
gdb ./simplicial_complex.db # run gdb
```

## Profiling

You can compile an executable for profiling with `gprof` by specifying a `*.prof` extension
```
make simplicial_complex.prof
./simplicial_complex.prof # generates gmon.out
gprof simplicial_complex.prof gmon.out > prof.txt # writes to prof.txt
gprof -q simplicial_complex.prof gmon.out > prof.txt # gall graph analysis
```

## Static analysis

An ongoing effort is to do static code analysis using `clang-tidy`.

The following will perform the static code analysis on `quickstart.cpp`.
```
clang-tidy quickstart.cpp -checks=clang-analyzer-*,openmp-*,bugprone-*,performance-*,cppcoreguidelines- -- -I../include
```
