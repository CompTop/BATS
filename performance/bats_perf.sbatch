#!/bin/sh

#SBATCH --job-name=bats-perf
#SBATCH --account=pi-bradnelson
#SBATCH --output=bats_perf.out
#SBATCH --error=bats_perf.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

module load gcc
make -j

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./scaling.out
./divide_conquer.out
./zigzag_circle.out
./zigzag_cube.out