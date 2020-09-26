# Performance

The code in this folder is used to compute performance information.

## Homology Computation

The computation of homology is embarassingly parallel.

## Divide and Conquer Parallelization

In many situations, the performance bottleneck in computing persistent or zigzag homology is the performance bottleneck.

To test whether we get a parallel speedup from the divide and conquer implementation (implemented by the funciton `barcode_sparse_divide_conquer`), we perform two tests:
1. A randomly generated quiver representation (`barcode_sparse_divide_conquer.cpp`)
2. We generate a barcode with fairly high dimensional homology in each vector space.
