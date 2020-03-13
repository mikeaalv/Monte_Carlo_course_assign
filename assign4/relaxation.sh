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
Lrange=(4 32)
Trange=(2.3 5.0)
for L in "${Lrange[@]}"
do
  for T in "${Trange[@]}"
  do
    for randi in "${nreplic[@]}"
    do
      # echo $randi
      ./metropois ${L} ${T} 50099 ${randi} 0.5 1 > output.L${L}.T${T}.rand${randi}.gener.tab
    done
  done
done


##tested on approximation of tau
