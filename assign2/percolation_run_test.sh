#!/bin/bash
#SBATCH --job-name=randpi      # Job name
#SBATCH --partition=fsr8602       # Partition (queue) name
#SBATCH --ntasks=1              # Single task job
#SBATCH --cpus-per-task=1       # Number of cores per task
#SBATCH --mem=2gb               # Total memory for job
#SBATCH --time=00:10:00         # Time limit hrs:min:sec
#SBATCH --output=log.%j         # Standard output and error log
#SBATCH --mail-user=Yue.Wu@uga.edu # mailing address
#SBATCH --mail-type=ALL         # mailing information
cd $SLURM_SUBMIT_DIR
module load icc/2018.1.163-GCC-6.4.0-2.28
# icc -O2 percolation.cpp -o percolation

# tests
# ./percolation 4 0 1 1 0
# ./percolation 4 1 1 1 0
# ./percolation 4 0.5 1 1 0
# ./percolation 4 0.5 2 1 0
# ./percolation 4 0.75 1 1 0
# ./percolation 4 0.75 2 1 0
# ./percolation 10 0 1 1 0
# ./percolation 10 1 1 1 0
# ./percolation 10 0.1 1 1 0
# ./percolation 10 0.9 1 1 0
# ./percolation 10 0.5 1 1 0
# ./percolation 10 0.7 1 1 0
# ./percolation 10 0.7 2 1 0

# ./percolation 10 0.7 2 1 1

# test time and L
# time ./percolation 1000 0.5 0 0 0
# time ./percolation 10000 0.5 0 0 0

# test p range
# randi=(1 2 3 4 5)
# for rand in "${randi[@]}"
# do
#   echo $rand
#   ./percolation 1000 0.6 ${rand} 0 0
# done

# randi=(1 2 3 4 5)
# for rand in "${randi[@]}"
# do
#   echo $rand
#   ./percolation 1000 0.59 ${rand} 0 0
# done

produce result needed for 1
produce result needed for 2
produce result needed for 3
