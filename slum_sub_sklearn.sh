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

bash generate_sklearn_classifiers.sh
