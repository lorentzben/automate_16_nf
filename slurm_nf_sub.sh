#!/bin/bash

#SBATCH --partition=batch
#SBATCH --job-name=automate_16_nf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem=8gb

#Replace this with your UGA email to get notified on completion
#SBATCH --mail-user="${USER}@uga.edu"
#SBATCH --mail-type=BEGIN,END,FAIL

#metadata/manifest name and outdir names
CURR_META="all_days_sbm_cec_metadata"
CURR_MANI="all_days_sbm_cec_raw_manifest"
ANALYSIS="clean_dir"
ME=$(whoami)

#create workdir and copy files over if required 
if [[ ! -e /scratch/${USER}/nf_dev/"$ANALYSIS" ]]; then
	mkdir /scratch/${USER}/nf_dev/"$ANALYSIS"
	cp /home/${USER}/applegate/villegas/compare_methods/"$CURR_MANI".txt /scratch/${USER}/nf_dev/"$ANALYSIS"/"$CURR_MANI".tsv
	cp /home/${USER}/applegate/villegas/compare_methods/"$CURR_META".txt /scratch/${USER}/nf_dev/"$ANALYSIS"/"$CURR_META".tsv
fi

#move into work dir and start analysis
cd /scratch/${USER}/nf_dev/"$ANALYSIS"

#Requires nextflow to be installed or loaded
module load Nextflow/22.04.5

#nextflow command, revision, forward, reverse and classifiers should be updated.
nextflow run lorentzben/automate_16_nf \
	-r 2.1.0 \
	--input /work/sealab/$(echo $USER)/applegate/villegas_data/RAW/raw_files \
	--metadata $CURR_META.tsv \
	--manifest $CURR_MANI.tsv \
	--itemOfInterest day \
	--forward 276 \
	--rev 217 \
	--classify_515 /home/$(echo $USER)/auto_16_nf/classifiers/515-806-classifier.qza \
  	--classify_full /home/$(echo $USER)/auto_16_nf/classifiers/16s-whole-seq-classifier.qza \
	--outdir /work/sealab/${USER}/nf_dev/"$ANALYSIS" \
	-c /home/$(echo $USER)/applegate/villegas/compare_methods/.nextflow/config/gacrc.config \
	-with-tower \
	-profile slurm,singularity \
	-resume


