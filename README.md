Automate 16s NF
-------------------------------------------------
This project's aim is focused on automating processes of 16s rRNA analysis and learn the pipeline development tool Nextflow. The ultimate aim is to make the process more accessible. 

## Prerequisities
* Linux based system (can be HPC)
* Nextflow
* docker or singularity installed

## Install

```shell
$ git clone https://github.com/lorentzben/automate_16s_nf.git
$ python3 -m pip install PySimpleGUI
```
## Slurm Scripts

I have included two slurm scripts: 

[One pulls all singularity images since nextflow can time out](slurm_pull_singularity.sh)

[One sets up work directory and then starts the analysis](slurm_nf_sub.sh)

The nextflow run command can look like this:

```bash
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
```

### Note:
TODO check to see if these links are correct for the version of Scikit that Qiime is using.
On setup, you will need to download one of the classifiers I generated below, or genertate your own. 

[515-806 Classifier](https://outlookuga-my.sharepoint.com/:u:/g/personal/bjl34716_uga_edu/ETXcJb8cC1VNnZUn2HGxAEcBCQZSet4635ZJENEd0TrDXA?e=YzPQcf)

[515-805 ALT](https://lorentzvault.quickconnect.to/d/s/lzaQnpI1s6n10gDy21VEJhTiloqNXuM1/X12oymNQrJQPRo3hmnf5qr1ydrFwk3Uj-Tr7A0QUqBAk)

[Whole 16s Classifier](https://outlookuga-my.sharepoint.com/:u:/g/personal/bjl34716_uga_edu/EfjwPSMRsuNJuJc4vsFsq48BmRg4Y7el79hdQpPB4KEGyQ?e=Nu7rZE)

[Whole 16s ALT](https://lorentzvault.quickconnect.to/d/s/lzaNTHGucSeauP4cgjJRBQoUlNg1coH7/jJMrGBAJd2cx2Ma841fxOhUIZxwx00aB-X77AxCEqBAk)

[Qiime Guide to Train Custom Classifer](https://docs.qiime2.org/2021.2/tutorials/feature-classifier/)

Also see generate_sklearn_classifers.sh for a slurm script.

I have a script in this repo that will also generate the files above, however it takes about 12 hours to run. 

For the item of interest and related column in the metadata file, QIIME will not let filenames use the '_' character so they should be replaced by the '-' or alternate character. 


### If you want to check the files for Depth Parameters

The file for cutoff is in:
$RESULT_DIRECTORY/qiime/demux_summary.qzv

The file for Sampling Depth and Rarefaction Depth is:
$RESULT_DIRECTORY/qiime/table.qzv

Both of these files can be unzipped or opened on the internet at view.qiime2.org

