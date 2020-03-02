#!/bin/bash
#SBATCH --job-name=percolation      # Job name
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
rm -r *.tab

SECONDS=0
#range plot also to calculate M
Lrange=(100 200 500)
prange=(0.0 0.4 0.5 0.55 0.58 0.585 0.59 0.595 0.6 0.7 0.8 1.0)
randseq=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
for L in "${Lrange[@]}"
do
  for p in "${prange[@]}"
  do
    for randi in "${randseq[@]}"
    do
      # echo $randi
      ./percolation ${L} ${p} ${randi} 0 0 >> rangeplot.tab
    done
  done
done
