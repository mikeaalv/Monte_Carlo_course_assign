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

# uncertainty estimation
Lrange_sub=(100 500)
prange_sub=(0.59 0.6)
randseq=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
for L in "${Lrange_sub[@]}"
do
  for p in "${prange_sub[@]}"
  do
    for repi in {1..4}
    do
      for randi in "${randseq[@]}"
      do
        let randreal=$randi+$repi*20
        # echo $randreal
        echo -n $repi" " >> rangeplot_uncernty.tab
        ./percolation ${L} ${p} ${randreal} 0 0 >> rangeplot_uncernty.tab
      done
    done
  done
done

#example cluster
Lrange_sub=(100 500)
prange_sub=(0.5 0.8)
randi=1
for L in "${Lrange_sub[@]}"
do
  for p in "${prange_sub[@]}"
  do
    # echo $randi
    ./percolation ${L} ${p} ${randi} 1 0 >> examp_plot${L}${p}.tab
  done
done

#range plot also to calculate M
Lrange=(100 300 500 700)
randseq=(1 2 3 4 5)
for L in "${Lrange[@]}"
do
  for p in `seq 0.58 0.0005 0.595`
  do
    for randi in "${randseq[@]}"
    do
      # echo $randreal
      echo -n $randi" " >> save_pc_calc.tab
      ./percolation ${L} ${p} ${randi} 0 0 >> save_pc_calc.tab
    done
  done
done

#histogram of cluster size
p=0.59
Lrange_sub=(100 500)
randseq=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
for L in "${Lrange_sub[@]}"
do
  for randi in "${randseq[@]}"
  do
    # echo $randi
    ./percolation ${L} ${p} ${randi} 0 1 >> tao_esti${L}.tab
  done
done

echo $SECONDS
