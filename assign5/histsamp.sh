#!/bin/bash
#SBATCH --job-name=metropolis      # Job name
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
# icc -O2 metropois.cpp -o metropois
# icc -O2 randomsampling.cpp -o randomsampling
# test L T N randseeds inip
# ./metropois 60 10 30000 1 1 10 > output

SECONDS=0
L=8
## simple random sampling at infinite temperature
nreplic=(1 2 3)
for randi in "${nreplic[@]}"
do
  # echo $randi
  ./randomsampling ${L} 5000099 ${randi} 0.5 1 > output.simplesampling.hist.rand${randi}.gener.tab
done
# echo $SECONDS
## metropolis at selected temperature
Trange=(1 2 2.2 2.5 2.7 3 4)
nreplic=(1 2 3 4 5 6 7 8 9)
for T in "${Trange[@]}"
do
  for randi in "${nreplic[@]}"
  do
    # echo $randi
    ./metropois ${L} ${T} 50099 ${randi} 0.5 1 > output.metropois.T${T}.rand${randi}.gener.tab
  done
done

## Metropolis at Tc
Tc=3.0
nreplic=(1 2 3 4 5 6 7 8 9)
for randi in "${nreplic[@]}"
do
  # echo $randi
  ./metropois ${L} ${Tc} 300099 ${randi} 0.5 1 > output.metropois.Tc.hist.rand${randi}.gener.tab
done
##

echo $SECONDS
