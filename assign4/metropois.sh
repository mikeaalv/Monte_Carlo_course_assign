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
# test L T N randseeds inip
# ./metropois 60 10 30000 1 1 10 > output

SECONDS=0
#different lattice size(L) and temperature(T)
Lrange=(5 20 60)
Trange=(1 2 2.5 3 4 10 19 19.5 20)
nreplic=(1 2 3)
for L in "${Lrange[@]}"
do
  for T in "${Trange[@]}"
  do
    for randi in "${nreplic[@]}"
    do
      # echo $randi
      ./metropois ${L} ${T} 5099 ${randi} 0.5 1 > output.L${L}.T${T}.rand${randi}.gener.tab
    done
  done
done

#different length run
L=20
Trange=(1 2 2.5 3 4 10 20)
nreplic=(1 2 3)
lenmc=(599 50099)
for len in "${lenmc[@]}"
do
  for T in "${Trange[@]}"
  do
    for randi in "${nreplic[@]}"
    do
      # echo $randi
      ./metropois ${L} ${T} ${len} ${randi} 0.5 1 > output.L${L}.T${T}.rand${randi}.mc${len}.mccomp.tab
    done
  done
done

# critical temperature estimation
Lrange=(20 30 40 50 60)
Trange=(2.4 2.45 2.5 2.55 2.6)
nreplic=(1 2 3 4 5 6 7 8 9 10)
for L in "${Lrange[@]}"
do
  for T in "${Trange[@]}"
  do
    # echo $T
    for randi in "${nreplic[@]}"
    do
      # echo $randi
      ./metropois ${L} ${T} 5099 ${randi} 0.5 1 > output.L${L}.T${T}.rand${randi}.Tc.tab
    done
  done
done
echo $SECONDS
