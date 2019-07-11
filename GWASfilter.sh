#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --output=GWASfilter.out
#SBATCH --job-name="GWASfilter"
#SBATCH --mail-user=danev1@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -p batch

cd $SLURM_SUBMIT_DIR

hostname

echo "rmarkdown::render('GWASfilter.Rmd')" | R --slave

