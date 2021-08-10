# Performance

The code in this folder is used to compute performance information.

## SLURM Jobs

You can submit a job on a slurm cluster using `bats_perf.sbatch`:
```
sbatch bats_perf.sbatch
```

The included script is setup to run on the UChicago midway3 cluster.  [User Guide](https://mdw3.rcc.uchicago.edu/user-guide/).
You will likely need to modify to run on another cluster.

## Run locally

You can `make -j` and then run any of the executables.

## Homology Computation

The computation of homology is embarassingly parallel.

## Divide and Conquer Parallelization

In many situations, the performance bottleneck in computing persistent or zigzag homology is the performance bottleneck.

To test whether we get a parallel speedup from the divide and conquer implementation (implemented by the funciton `barcode_sparse_divide_conquer`), we perform two tests:
1. A randomly generated quiver representation (`barcode_sparse_divide_conquer.cpp`)
2. We generate a barcode with fairly high dimensional homology in each vector space.

