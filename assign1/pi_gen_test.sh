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
# icc -O2 pi_gen.cpp -o pi_gen
SECONDS=0
randseedstart=0
# different sample size
N_sample_vec=(100 1000 10000 100000 1000000 10000000 100000000 1000000000)
for N_sample_ele in "${N_sample_vec[@]}"
do
  echo $N_sample_ele
  ./pi_gen ${N_sample_ele} ${randseedstart} > output_${N_sample_ele}.tab
done
echo $SECONDS
