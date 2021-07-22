#!/bin/sh
#SBATCH --job-name=BiSSECetacean  #Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, FAIL)
#SBATCH --mail-user=nhans@ufl.edu # Where to send mail
#SBATCH --nodes=1  # Run all processes on a single node	
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=8 # Number of cores: Can also use -c=4
#SBATCH --mem-per-cpu=1gb # Per processor memory
#SBATCH -t 5:00:00     # Walltime
#SBATCH -o BiSSECetacean.%j.out # Name output file
#SBATCH
#SBATCH --account=burleigh
#SBATCH --qos=burleigh-b

pwd; hostname; date

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# change the location here
cd /blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Run2/

module load R

# change the script name and path
Rscript /blue/burleigh/nhans/GitHubNH/DataCloning/BiSSE/Scripts/Cetaceanfit_forcluster_CORE.R 

date 
