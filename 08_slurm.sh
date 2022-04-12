#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=08_slurm
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --time=96:00:00
#SBATCH --mem=64gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user="bjl34716@uga.edu"
#SBATCH --mail-type=BEGIN,END,FAIL

singularity run docker://lorentzb/r_08 bash 08_make_report.sh
