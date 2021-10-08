#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=yadav_analysis
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --time=96:00:00
#SBATCH --mem=64gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user="bjl34716@uga.edu"
#SBATCH --mail-type=BEGIN,END,FAIL

module load Nextflow/20.04.1
nextflow run main.nf --input seqs --metadata metadata_1.tsv --manifest mapping.tsv --itemOfInterest treatment --outdir yadav