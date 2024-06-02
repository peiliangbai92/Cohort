#!/bin/bash
#SBATCH --account=statistics
#SBATCH --qos=statistics
#SBATCH --job-name=S1                 # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=baipl92@ufl.edu   # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=144:00:00              # Time limit hrs:min:sec
#SBATCH --output=S1.log               # Standard output and error log
pwd; hostname; date

module load R

echo "Running R script on a single CPU core"

R CMD BATCH S1.R

date
