#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=Generate_sklearn_classifiers
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem=64gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user="bjl34716@uga.edu"
#SBATCH --mail-type=BEGIN,END,FAIL

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

echo "this script will take approx 11:16:54 to generate 2 classifiers"
singularity run docker://lorentzb/automate_16_nf bash generate_sklearn_classifiers.sh
