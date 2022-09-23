#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=pull_singularity_images
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --mem=8gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user="${USER}@uga.edu"
#SBATCH --mail-type=BEGIN,END,FAIL

SINGULARITY_CACHEDIR="/scratch/${USER}/singularity"

singularity pull docker://lorentzb/automate_16_nf:2.0
singularity pull docker://lorentzb/py2_env:2.0
singularity pull docker://lorentzb/py2_test:2.0
singularity pull docker://lorentzb/r_01:2.0
singularity pull docker://lorentzb/r_02:2.0
singularity pull docker://lorentzb/r_03:2.0
singularity pull docker://lorentzb/r_04:2.0
singularity pull docker://lorentzb/r_05:2.0
singularity pull docker://lorentzb/r_06:2.0
singularity pull docker://lorentzb/r_07:2.0
singularity pull docker://lorentzb/r_08:2.0
singularity pull docker://lorentzb/r_09:2.0
singularity pull docker://lorentzb/r_10:2.0
singularity pull docker://lorentzb/r_11:2.0
singularity pull docker://lorentzb/r_12:2.0
singularity pull docker://lorentzb/r_13:2.0
singularity pull docker://lorentzb/r_14:2.0
singularity pull docker://lorentzb/r_latest_2:2.0
singularity pull docker://lorentzb/r_latest:2.0
singularity pull docker://lorentzb/qiime2lefse:1.0
singularity pull docker://lorentzb/py2_env:1.0
