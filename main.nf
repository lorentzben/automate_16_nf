#!/usr/bin/env nextflow

def helpMessage(){
	log.info"""
	Usage:

    This script should be called as follows:
    nextflow run automate_16_nf/analysis_16 --input seqs --metadata metadata.tsv --manifest manifest.tsv --item_of_interest "Day"

    Main arguments:
        --input [path/to/folder]      Folder containing demultiplexed fastq.gz files
        --metadata [path/to/file]     Path to metadata sheet in tsv format see EXAMPLE_METADATA.tsv
        --manifest [path/to/file]     Path to mapping file in tsv format see EXAMPLE_MAPPING.tsv 
        --itemOfInterest [str]        Item of interest, group defining treatment vs control or longitudinal variable
        --name [str]                   Name for the analysis run, if not provided nextflow will generate one
        --sampDepth [str]             Value of sampling depth derived from demux_summary.qzv
        --rareDepth [str]             Value of rarefaction depth derived from table.qzv
        --outdir [file]               The output directory where the results will be saved 


    """.stripIndent()
}

//show help message 
if (params.help){
    helpMessage()
    exit 0 
}


if(params.manifest) {
    tsvFile = file(params.manifest).getName()
    Channel
        .fromPath(params.manifest)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .into{ ch_single_pair ; ch_make_qiime ; ch_mani_veri }
}

if(params.input){
    Channel
        .fromPath(params.input)
        .ifEmpty {exit 1, log.info "Cannot find path file ${input}"}
        .into{ ch_make_qiime_seq ; ch_seqs_veri }
}

if(params.metadata) {
    tsvFile = file(params.metadata).getName()
    Channel
        .fromPath(params.metadata)
        .ifEmpty { exit 1, log.info "Cannot find path file ${tsvFile}"}
        .into{ ch_meta_veri ; ch_meta_feature_viz; ch_alpha_metadata ; ch_metadata_rare_curve ; ch_metadata_alpha_sig ; ch_metadata_beta_sig ; ch_metadata_phylo_tree ; ch_metadata_lefse ; ch_metadata_finalize}
}

if(params.sampDepth){
    Channel
        .from(params.sampDepth)
        .set{ ch_user_sample_depth }
}

if(params.rareDepth){
    Channel
        .from(params.rareDepth)
        .set{ ch_user_rarefaction_depth }
}

if(!params.sampDepth){
    Channel
        .from(0)
        .set{ ch_user_sample_depth }
}

if(!params.rareDepth){
    Channel
        .from(0)
        .set{ ch_user_rarefaction_depth }
}

Channel
    .fromPath("${baseDir}/plot_cladogram.py")
    .ifEmpty {exit 1, log.info "Cannot find file plot_cladogram.py!"}
    .set{ ch_clado_file }

Channel
    .fromPath("${baseDir}/plot_res.py")
    .ifEmpty {exit 1, log.info "Cannot find file plot_res.py!"}
    .set{ ch_plot_res }

Channel
    .fromPath("${baseDir}/16s-whole-seq-classifier.qza")
    .ifEmpty {exit 1, log.info "Cannot find the classifier!"}
    .set{ ch_whole_classifier}
Channel
    .fromPath("${baseDir}/515-806-classifier.qza")
    .ifEmpty {exit 1, log.info "Cannot find the classifier!"}
    .set{ ch_515_classifier }

Channel
    .from(params.itemOfInterest)
    .ifEmpty {exit 1, log.info "Cannot find Item of interest"}
    .into{ ch_ioi_veri ; ch_ioi_beta_sig ; ch_ioi_phylo_tree ; ch_ioi_lefse ; ch_ioi_denoise_to_file }

Channel
    .fromPath("${baseDir}/graph.sh")
    .set{ ch_graph_script } 

Channel
    .fromPath("${baseDir}/qiime_to_lefse.R")
    .set { ch_lefse_format_script }

Channel
    .fromPath("${baseDir}/filter_samples.py")
    .set{ ch_filter_script }

Channel 
    .fromPath("${baseDir}/lefse_analysis.sh")
    .set{ ch_lefse_analysis_script }
Channel
    .fromPath("${baseDir}/report.Rmd")
    .set{ ch_report_outline }

Channel
    .fromPath("${baseDir}/make_report.sh")
    .set{ ch_report_bash_script }

Channel
    .fromPath("${baseDir}/init_and_refresh.R")
    .set{ ch_r_init }

Channel
    .fromPath("${baseDir}/renv.lock")
    .set{ ch_r_lock }


/*
process SetupPy2CondaEnv{
    //conda "${projectDir}/python2_env.yml"
    //conda "python2_env.yml"
    container "lorentzb/py2_env"

    input:
    file plot_clado from ch_clado_file
    file plot_res from ch_plot_res
    

    shell:
    '''
    #!/usr/bin/env bash
    LEFSE=$(which lefse.py)
    LEFSE_DIR=${LEFSE::-8}
    cp plot_res.py $LEFSE_DIR
    cp plot_cladogram.py $LEFSE_DIR
    '''

}
*/

//I removed Biocmanager and Microbiome, so if a function breaks, thats why. 
process SetupRPackages{
    //conda "${projectDir}/r_env.yml"
    //conda "r_env.yml"
    //label 'r'
    container "docker://lorentzb/r_latest"

    output:
    file "set.txt" into ch_r_wait

    script:
    """
    #!/usr/bin/env Rscript --vanilla
    #Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
    if(!require(rmarkdown)) {install.packages("rmarkdown", repos="http://cran.us.r-project.org")}
    if(!require(renv)) {install.packages("renv",repos="http://cran.us.r-project.org")}
    renv::init()
    if(!require(remotes)){install.packages("remotes",repos="http://cran.us.r-project.org")}
    if(!require(devtools)){install.packages("devtools",repos="http://cran.us.r-project.org")}
    if(!require(jamba)){remotes::install_github("jmw86069/jamba@0.0.6.900")}
    remotes::install_github("tidyverse/ggplot2@v3.3.2")
    remotes::install_github("vegandevs/vegan@v2.5-7")
    if(!require(ampvis2)){remotes::install_github("MadsAlbertsen/ampvis2@2.6.8")}
    if(!require(ggvegan)){remotes::install_github("gavinsimpson/ggvegan@4bc6ee9945dd9229ed486409c0acab9413b8c9af")}
    if(!require(ggConvexHull)){remotes::install_github("cmartin/ggConvexHull@660f4094da44dd500c3c0684b9c5c20c21ee823a")}
    
    cat('done',file='set.txt', sep='\n')
    """

}

//TODO write item of interest into csv/txt for r script
process VerifyManifest{
    publishDir "${params.outdir}", mode: 'copy'
    input:
    file manifest from ch_mani_veri
    path seqs_dir from ch_seqs_veri
    file metadata from ch_meta_veri
    val ioi from ch_ioi_veri

    output:

    file "order_item_of_interest.csv" into ch_format_ioi_order

    /*this is in place for local deployment, but the server does not give access to the dir for some reason
    The change is nessecary to do nextflow run -r main lorentzben/automate_16_nf
    */
    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    script:
    """
    #!/usr/bin/env python3
    import os 
    from pathlib import Path
    from pathlib import PurePath
    import pandas as pd 
    import csv 

    try:
        read_metadata = pd.read_table('${metadata}', index_col=0, sep='\t')
    except FileNotFoundError:
        exit(1)

    try:
        read_order = pd.read_table('${baseDir}/order_item_of_interest.csv', index_col=0, sep=',')
    except FileNotFoundError:
        iois = list(pd.Series.unique(read_metadata['${ioi}']))
        ioisdf = pd.DataFrame(iois[1:])
        ioisdf.columns = ['${ioi}']
        pd.DataFrame.to_csv(ioisdf, 'order_item_of_interest.csv', index=False)

    seq_dir = '${seqs_dir}'
    try:
        read_manifest = pd.read_table('${manifest}', index_col=0, sep='\t')
    except FileNotFoundError:
        exit(1)

    # sets current dir and finds the fastq and fastq.gz files in the current directory
    p = Path.cwd()
    list_of_fastq = list(p.glob(seq_dir + '/*.fastq'))
    list_of_gz = list(p.glob(seq_dir+'/*.fastq.gz'))

    fastq_files = []
    gz_files = []
    found = []
    missing = []

    # pulls only the filename and saves to a list
    for item in list_of_fastq:
        filename = os.path.split(item)[1]
        fastq_files.append(filename)

    for item in list_of_gz:
        filename = os.path.split(item)[1]
        gz_files.append(filename)
    if read_manifest.columns[0] == 'forward-absolute-filepath':
        # iterates over the forward reads and then the reverse reads to check to make sure they are all accounted for
        try:
            for item in read_manifest['forward-absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    
                    found.append(filename)
                else:
                    if filename in gz_files:
                        
                        found.append(filename)
                    else:
                        
                        missing.append(filename)
        except KeyError:
            print('single read project')

        # try except in the case that the user only has single end reads.
        try:
            for item in read_manifest['reverse-absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    
                    found.append(filename)
                else:
                    if filename in gz_files:
                        
                        found.append(filename)
                    else:
                        
                        missing.append(filename)

        except KeyError:
            print("looking for forward only reads")
    else:
        # this case is if there are only single reads and after which we can figure that the manifest file is wrong
        try:
            for item in read_manifest['absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    
                    found.append(filename)
                else:
                    if filename in gz_files:
                        
                        found.append(filename)
                    else:
                        
                        missing.append(filename)
        except KeyError:
            print("headings in the manifest appear to be incorrect")
            exit(1)

    # if missing is not an empty list, i.e. a file listed in the manifest is not detected, it raises an error and
    # creates a list for the user
    if missing != []:

        print("files are missing, please see missing.csv to correct")

        with open('missing.csv', 'w', newline='') as csvfile:
            fieldnames = ['filename']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for filename in missing:
                writer.writerow({'filename': filename})

        exit(0)


    print("the manifest called: " + '${manifest}' +
                 " is valid and ready to go")
    """

}

process CheckSinglePaired { 

    publishDir "${params.outdir}", mode: 'copy'
    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"
    
    input: 
    file manifest from ch_single_pair
    val ioi from ch_ioi_denoise_to_file

    output: 
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType
    file "item_of_interest.csv" into ch_ioi_file_out


    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import os 
    import csv

    
    with open('item_of_interest.csv', 'w', newline='') as csvfile:
        fieldnames = ['item name']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'item name': '${ioi}'})

    read_manifest = pd.read_table('${manifest}', index_col=0, sep='\t+', engine='python')

    if read_manifest.columns[0] == 'absolute-filepath':
        print("single end analysis")
        format = "SingleEndFastqManifestPhred33V2"
        data = "SampleData[SequencesWithQuality]"

        with open("manifest_format.txt", "w") as file:
            file.write(format)

        with open("data_type.txt", "w") as d_file:
            d_file.write(data)
        print(format + " " + data)

        
    elif read_manifest.columns[0] == 'forward-absolute-filepath':
        print("paired end analysis")
        format = "PairedEndFastqManifestPhred33V2"
        data = "SampleData[PairedEndSequencesWithQuality]"
        with open("manifest_format.txt", "w") as file:
            file.write(format)

        with open("data_type.txt", "w") as d_file:
            d_file.write(data)
        print(format + " " + data)
    else:
        print(
            "cannot determine if paired or single end, check manifest file")
        exit(1)

    
    """
    
}

process GenerateSeqObject{

    publishDir "${params.outdir}/qiime", mode: 'copy'
    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input: 
    file manifest from ch_make_qiime
    file manifest_format from manifest_type
    file data_type from dataType
    path seqs from ch_make_qiime_seq

    output: 
    file 'demux.qza' into ch_qiime_obj
    file manifest_format into ch_manifest_type
    

    shell:
    '''
    DAT=$(head !{data_type})
    MANI=$(head !{manifest_format})
    #module load  QIIME2/2020.11
    qiime tools import \
    --type $DAT\
    --input-path !{manifest} \
    --output-path demux.qza \
    --input-format $MANI
    '''
}

process QualControl{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    file seq_obj from ch_qiime_obj
    
    output: 
    file('demux_summary/*') into ch_qiime_qual
    file seq_obj into ch_qiime_denoise

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    script:
    """
    #!/usr/bin/env bash

    qiime demux summarize \
    --i-data ${seq_obj} \
    --o-visualization demux_summary.qzv

    qiime tools export \
    --input-path demux_summary.qzv \
    --output-path demux_summary/
    """

}

process FindCutoffs{
    //TODO Remove this call after debugging
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    file 'manifest_format.txt' from ch_manifest_type
    file('demux_summary/*') from ch_qiime_qual
    
    output: 
    file("cutoffs.csv") into ch_cutoff_vals
    file("manifest_format.txt") into ch_manifest_type_denoise

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    from pathlib import Path
    import numpy as np 
    import csv 

    wd = Path.cwd()

    seq_file = pd.read_table("manifest_format.txt")
    if seq_file.columns[0] == "SingleEndFastqManifestPhred33V2":
        seq_format = "single"
    else:
        seq_format = "paired"

    def find_cutoffs(dataframe):
        mean_qual = dataframe[4:5]
        mean_count = dataframe[0:1]

        average_qual = np.round(mean_qual.mean(axis=1), 0)
        average_count = np.round(mean_count.mean(axis=1), 0)

        mean_qual_vals = np.array(mean_qual)[0]
        mean_count_vals = np.array(mean_count)[0]

        if int(average_qual) < 30:
            print(
                "The Average Quality of these sequences may be a concern would you like to continue?")
            exit(0)

        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[i] >= int(average_qual):
                if mean_count_vals[i] >= int(average_count):
                    left_cutoff = i+1
                    break
        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
                if mean_count_vals[len(mean_count_vals)-1-i] >= int(average_count):
                    right_cutoff = len(mean_qual_vals)-i
                    break
        return(left_cutoff, right_cutoff)

    def find_rev_cutoffs(dataframe):
        mean_qual = dataframe[4:5]
        mean_count = dataframe[0:1]

        average_qual = np.round(mean_qual.mean(axis=1), 0)+2
        average_count = np.round(mean_count.mean(axis=1), 0)

        mean_qual_vals = np.array(mean_qual)[0]
        mean_count_vals = np.array(mean_count)[0]

        if int(average_qual) < 30:
            print(
                "The Average Quality of these sequences may be a concern would you like to continue?")
            exit(0)

        left_cutoff = 0
        right_cutoff = len(mean_qual_vals)-1

        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[i] >= int(average_qual):
                if mean_count_vals[i] >= int(average_count):
                    left_cutoff = i+1
                    break

        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
                if mean_count_vals[len(mean_count_vals)-1-i] >= int(average_count):
                    right_cutoff = len(mean_qual_vals)-i
                    break

        return(left_cutoff, right_cutoff)

    if seq_format == "single":
        print("determining left and right cutoffs based on qual score")

        input_file = str(wd)+"/demux_summary/forward-seven-number-summaries.tsv"

        summary = pd.read_table(input_file, index_col=0, sep='\t')
        left_cutoff, right_cutoff = find_cutoffs(summary)

        print("right cutoff: "+str(right_cutoff))
        print("left cutoff: " + str(left_cutoff))

        with open('cutoffs.csv', 'w', newline='') as csvfile:
            fieldnames = ['cutoff', 'value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow({'cutoff': 'left', 'value': left_cutoff})
            writer.writerow({'cutoff': 'right', 'value': right_cutoff})
            writer.writerow({'cutoff': 'filename', 'value': input_file})

        print(left_cutoff, right_cutoff)

    elif seq_format == "paired":
        print("determining forward and revese, left and right cutoffs based on qual score")
        forward_file = str(wd)+"/demux_summary/forward-seven-number-summaries.tsv"
        fr_summary = pd.read_table(forward_file, index_col=0, sep='\t')

        forward = find_cutoffs(fr_summary)

        reverse_file = str(wd)+"/demux_summary/reverse-seven-number-summaries.tsv"
        rev_summary = pd.read_table(reverse_file, index_col=0, sep='\t')

        reverse = find_rev_cutoffs(rev_summary)

        print("forward cutoffs: "+str(forward))
        print("reverse cutoffs: " + str(reverse))

        with open('cutoffs.csv', 'w', newline='') as csvfile:
            fieldnames = ['cutoff', 'value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow({'cutoff': 'forward left', 'value': forward[0]})
            writer.writerow({'cutoff': 'forward right', 'value': forward[1]})
            writer.writerow({'cutoff': 'reverse left', 'value': reverse[0]})
            writer.writerow({'cutoff': 'reverse right', 'value': reverse[1]})
            writer.writerow({'cutoff': 'filename', 'value': forward_file})
            writer.writerow({'cutoff': 'filename', 'value': reverse_file})

    """
}


process Denoise {
    publishDir "${params.outdir}/qiime", mode: 'copy'
    input:
    file seq_object from ch_qiime_denoise
    file("cutoffs.csv") from ch_cutoff_vals
    file("manifest_format.txt") from ch_manifest_type_denoise
    
    
    output:
    file "rep-seqs-dada2.qza" into ch_rep_seqs
    file "table-dada2.qza" into ch_table
    file "stats-dada2.qza" into ch_dada2_stats
    
    

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd 
    from pathlib import Path
    import numpy as np 
    import csv 
    import subprocess

    wd = Path.cwd()


    seq_file = pd.read_table("manifest_format.txt")
    if seq_file.columns[0] == "SingleEndFastqManifestPhred33V2":
        seq_format = "single"
    else:
        seq_format = "paired"
    cutoff = pd.read_table("cutoffs.csv", sep=",")    
    if seq_format == 'single':
        left = cutoff['value'][0]
        right = cutoff['value'][1]
        command = "qiime dada2 denoise-single \
            --i-demultiplexed-seqs demux.qza \
            --p-trim-left " + str(left)+" \
            --p-trunc-len " + str(right) + " \
            --o-representative-sequences rep-seqs-dada2.qza \
            --o-table table-dada2.qza \
            --o-denoising-stats stats-dada2.qza"
    elif seq_format == 'paired':
        forward_left = cutoff['value'][0]
        forward_right = cutoff['value'][1]
        rev_left = cutoff['value'][2]
        rev_right = cutoff['value'][3]
        command = "qiime dada2 denoise-paired \
            --i-demultiplexed-seqs demux.qza \
            --p-trunc-len-f " + str(forward_right)+" \
            --p-trunc-len-r " + str(rev_right) + " \
            --p-trim-left-f " + str(forward_left)+" \
            --p-trim-left-r " + str(rev_left) + " \
            --o-representative-sequences rep-seqs-dada2.qza \
            --o-table table-dada2.qza \
            --o-denoising-stats stats-dada2.qza"

    subprocess.run([command], shell=True)
    
    """
}

process FeatureVisualization{
    publishDir "${params.outdir}/qiime", mode: 'copy'
    input:
    file "stats-dada2.qza" from ch_dada2_stats
    file "table-dada2.qza" from ch_table
    file metadata_file from ch_meta_feature_viz
    file "rep-seqs-dada2.qza" from ch_rep_seqs

    output:
    file "stats-dada2.qzv" into ch_dada_stats_export
    file "table.qzv" into ch_table_viz_obj
    file "rep-seqs.qzv" into ch_req_seq_vis_obj
    file "rep-seqs-dada2.qza" into ch_rep_seq_tree_gen
    file "table-dada2.qza" into ch_alpha_div_table

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    script:
    """
    #!/usr/bin/env bash

    qiime metadata tabulate \
    --m-input-file stats-dada2.qza \
    --o-visualization stats-dada2.qzv

    qiime feature-table summarize \
    --i-table table-dada2.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file ${metadata_file}

    qiime feature-table tabulate-seqs \
    --i-data rep-seqs-dada2.qza \
    --o-visualization rep-seqs.qzv
    """
}

process TreeConstruction{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"
    
    input:
    file "rep-seqs-dada2.qza" from ch_rep_seq_tree_gen

    output:
    file "aligned-rep-seqs.qza" into ch_aligned_rep_seqs
    file "masked-aligned-rep-seqs.qza" into ch_mask_align_rep_seq
    file "unrooted-tree.qza" into ch_unrooted_tree
    file "rooted-tree.qza" into ch_rooted_tree
    file "rep-seqs-dada2.qza" into ch_rep_seq_classify

    script:
    """
    #!/usr/bin/env bash 

    qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs-dada2.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza 
    """
}

process ExportTable{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"
    
    input:
    file "table.qzv" from ch_table_viz_obj

    output:
    path "table_viz/*" into ch_table_viz_dir
    path "table_viz/*" into ch_table_viz_dir_rare

    script:
    """
    #!/usr/bin/env bash

    qiime tools export \
    --input-path table.qzv \
    --output-path table_viz
    """

}


process DetermineDepth{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    path "table_viz/*" from ch_table_viz_dir

    output:
    file "sampling_depth.csv" into ch_sampling_depth_csv
    file "samp_depth_simple.txt" into ch_depth


    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np 
    import csv

    input_file = "table_viz/sample-frequency-detail.csv"

    features = pd.read_csv(input_file, index_col=0, header=None)

    total_count = sum(features[1])
    sampling_depth = 0
    perc_features_retain = 0.0

    print("total count: " + str(total_count))

    feature_array = np.array(features)

    for i in range(len(feature_array)-1, -1, -1):
        sampling_depth = feature_array[i][0]
        perc_features_retain = ((sampling_depth * (i+1))/total_count)
        if perc_features_retain > .22:
            sampling_depth = feature_array[i][0]
            print("sampling depth: " + str(sampling_depth) + " % features retained: " +
                  str(round(perc_features_retain, 3)) + " samples retained: " + str(i))
            break
    print("writing dept out to file")
    with open('sampling_depth.csv', 'w', newline='') as csvfile:
        fieldnames = ['stat', 'value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow({'stat': 'sampling_depth', 'value': sampling_depth})
        writer.writerow({'stat': '%_features_retained',
                         'value': perc_features_retain})
        writer.writerow({'stat': 'filename', 'value': input_file})

    print("sampling_depth: "+str(sampling_depth))
    print("%_features_retained: " + str(perc_features_retain))
    
    with open('samp_depth_simple.txt', 'w') as file:
        file.write(str(int(sampling_depth)))
    """

}

process AlphaDiversityMeasure{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    file metadata from ch_alpha_metadata
    file "table-dada2.qza" from ch_alpha_div_table
    file "rooted-tree.qza" from ch_rooted_tree
    file "samp_depth_simple.txt" from ch_depth
    val user_depth from ch_user_sample_depth

    output:
    path "core-metric-results/*" into ch_core_beta_significance 
    path "core-metric-results/*" into ch_core_report
    file "shannon.qza" into ch_shannon_qza
    file "simpson.qza" into ch_simpson_qza 
    file "chao1.qza" into ch_chao_qza
    file "ace.qza" into ch_ace_qza
    file "obs.qza" into ch_obs_qza
    file "faith_pd.qza" into ch_faith_qza
    file "table-dada2.qza" into ch_table_rare_curve
    file "rooted-tree.qza" into ch_tree_rare_curve
    

    

    shell:
    '''
    #!/usr/bin/env bash

    if [ !{user_depth} == 0 ];then
        SAMP_DEPTH=$(head samp_depth_simple.txt)
    fi
    
    if [ !{user_depth} != 0 ];then
        SAMP_DEPTH=!{user_depth}
    fi
  
    qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table-dada2.qza \
    --p-sampling-depth $SAMP_DEPTH \
    --m-metadata-file !{metadata} \
    --output-dir core-metric-results 

    qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric shannon \
    --o-alpha-diversity shannon.qza

    qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric simpson \
    --o-alpha-diversity simpson.qza 

    qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric chao1 \
    --o-alpha-diversity chao1.qza

    qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric ace \
    --o-alpha-diversity ace.qza

    qiime diversity alpha \
    --i-table table-dada2.qza \
    --p-metric observed_features \
    --o-alpha-diversity obs.qza 

    qiime diversity alpha-phylogenetic \
    --i-table table-dada2.qza \
    --i-phylogeny rooted-tree.qza \
    --p-metric faith_pd \
    --o-alpha-diversity faith_pd.qza 
    '''
}

process AssignTaxonomy{
    //TODO change out the classifier for the 515 only one 
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    file "rep-seqs-dada2.qza" from ch_rep_seq_classify
    file "16s-whole-seq-classifier.qza" from ch_whole_classifier
    file "515-806-classifier.qza" from ch_515_classifier

    output:
    file "taxonomy.qza" into ch_taxonomy_phylo_tree
    file "taxonomy.qzv" into ch_classified_qzv
    

    script:
    """
    #!/usr/bin/env bash

    if [ ! -f "16s-whole-seq-classifier.qza" ]; then 
    echo "Error, download the classifier from readme"
    exit 1
    fi
    if [ ! -f "515-806-classifier.qza" ]; then 
    echo "Error, download the classifier from readme"
    exit 1
    fi

    qiime feature-classifier classify-sklearn \
    --i-classifier 16s-whole-seq-classifier.qza \
    --i-reads rep-seqs-dada2.qza \
    --p-confidence 0.6 \
    --o-classification taxonomy.qza

    qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv
    """
}

process CalcRareDepth{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    path "table_viz/*" from ch_table_viz_dir_rare

    output:
    file "rare_depth.txt" into ch_rare_curve_depth

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np 

    sample_freq = pd.read_csv("table_viz/sample-frequency-detail.csv")
    depth = sample_freq.median()[0]

    with open("rare_depth.txt",'w') as file:
        file.write(str(int(depth)))
    
    """

}

process RareCurveCalc{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    file "rare_depth.txt" from ch_rare_curve_depth
    file metadata from ch_metadata_rare_curve
    file "table-dada2.qza" from ch_table_rare_curve
    file "rooted-tree.qza" from ch_tree_rare_curve
    val user_rare_depth from ch_user_rarefaction_depth


    output:
    file "alpha-rarefaction.qzv" into ch_alpha_rare_obj
    path "alpha-rareplot/*" into ch_alpha_rare_viz
    file "table-dada2.qza" into ch_table_phylo_tree
    file "rooted-tree.qza" into ch_tree_lefse
    

    shell:
    '''
    #!/usr/bin/env bash

    if [ !{user_rare_depth} == 0 ];then
        DEPTH=$(head rare_depth.txt)
    fi
    
    if [ !{user_rare_depth} != 0 ];then
        DEPTH=!{user_rare_depth}
    fi

    
    
    qiime diversity alpha-rarefaction \
    --i-table table-dada2.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth $DEPTH \
    --m-metadata-file !{metadata} \
    --o-visualization alpha-rarefaction.qzv 

    qiime tools export \
    --input-path alpha-rarefaction.qzv \
    --output-path alpha-rareplot

    '''
}

process AlphaDiversitySignificance{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    file metadata from ch_metadata_alpha_sig
    file "shannon.qza" from ch_shannon_qza
    file "simpson.qza" from ch_simpson_qza
    file "chao1.qza" from ch_chao_qza
    file "ace.qza" from ch_ace_qza
    file "obs.qza" from ch_obs_qza
    file "faith_pd.qza" from ch_faith_qza

    output:
    path "shannon/*" into ch_shannon_path
    path "simpson/*" into ch_simpson_path
    path "chao1/*" into ch_chao_path
    path "ace/*" into ch_ace_path
    path "obs/*" into ch_obs_path
    path "faith_pd/*" into ch_faith_path

    script:
    """
    #!/usr/bin/env bash

    qiime diversity alpha-group-significance \
    --i-alpha-diversity shannon.qza \
    --m-metadata-file ${metadata} \
    --o-visualization shannon.qzv 

    qiime diversity alpha-group-significance \
    --i-alpha-diversity simpson.qza \
    --m-metadata-file ${metadata} \
    --o-visualization simpson.qzv

    qiime diversity alpha-group-significance \
    --i-alpha-diversity chao1.qza \
    --m-metadata-file ${metadata} \
    --o-visualization chao1.qzv 

    qiime diversity alpha-group-significance \
    --i-alpha-diversity ace.qza \
    --m-metadata-file ${metadata} \
    --o-visualization ace.qzv 

    qiime diversity alpha-group-significance \
    --i-alpha-diversity obs.qza \
    --m-metadata-file ${metadata} \
    --o-visualization obs.qzv 

    qiime diversity alpha-group-significance \
    --i-alpha-diversity faith_pd.qza \
    --m-metadata-file ${metadata} \
    --o-visualization faith_pd.qzv 

    qiime tools export \
    --input-path shannon.qzv \
    --output-path shannon

    qiime tools export \
    --input-path simpson.qzv \
    --output-path simpson

    qiime tools export \
    --input-path chao1.qzv \
    --output-path chao1

    qiime tools export \
    --input-path ace.qzv \
    --output-path ace 

    qiime tools export \
    --input-path obs.qzv \
    --output-path obs

    qiime tools export \
    --input-path faith_pd.qzv \
    --output-path faith_pd
    """
}

process BetaDiversitySignificance{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    val ioi from ch_ioi_beta_sig
    file metadata from ch_metadata_beta_sig
    path "core-metric-results/*" from ch_core_beta_significance 

    output:
    path "unweighted-sig/*" into ch_u_unifrac_beta_path
    path "weighted-sig/*" into ch_w_unifrac_beta_path

    script:
    """
    #!/usr/bin/env bash

    qiime diversity beta-group-significance \
    --i-distance-matrix core-metric-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file ${metadata} \
    --m-metadata-column ${ioi} \
    --o-visualization unweighted-unifrac-${ioi}-significance.qzv \
    --p-pairwise

    qiime tools export \
    --input-path unweighted-unifrac-${ioi}-significance.qzv \
    --output-path unweighted-sig/

    qiime diversity beta-group-significance \
    --i-distance-matrix core-metric-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file ${metadata} \
    --m-metadata-column ${ioi} \
    --o-visualization  weighted-unifrac-${ioi}-significance.qzv \
    --p-pairwise

    qiime tools export \
    --input-path weighted-unifrac-${ioi}-significance.qzv \
    --output-path weighted-sig/
    """
}

process GeneratePhylogeneticTrees{
    publishDir "${params.outdir}/graphlan", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    file metadata from ch_metadata_phylo_tree
    val ioi from ch_ioi_phylo_tree
    file "table-dada2.qza" from ch_table_phylo_tree
    file "taxonomy.qza" from ch_taxonomy_phylo_tree
    file "graph.sh" from ch_graph_script
    file "filter_samples.py" from ch_filter_script

    output:
    path "phylo_trees/*" into ch_png_phylo_tree
    file "table-dada2.qza" into ch_table_lefse
    file "taxonomy.qza" into ch_tax_lefse
    

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import pandas as pd
    import numpy as np 
    import time

    metadata_table= pd.read_table(\"${metadata}\", sep='\t')
    metadata_table = metadata_table.drop([0,1])

    ioi_set = set(metadata_table[\"${ioi}\"])
    ioi = '${ioi}'

    subprocess.run(['mkdir phylo_trees'], shell=True)

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:

        # filters/splits the feature table based on the current ioi
        
        filter_command = "python3 filter_samples.py -m ${metadata} -i ${ioi} -c "+item
        result = subprocess.run([filter_command], shell=True)

        time.sleep(2)

        # adds taxonomic info needed for plotting
        collapse_command = 'qiime taxa collapse \
        --i-table '+item+'-filtered-table.qza \
        --o-collapsed-table collapse-'+item+'-table.qza \
        --p-level 7 \
        --i-taxonomy taxonomy.qza'

        result = subprocess.run([collapse_command], shell=True)

        # exports artifact so that the next step can collect it
        export_command='qiime tools export \
        --input-path collapse-'+item+'-table.qza \
        --output-path collapse-'+item+'-frequency/'
        
        result = subprocess.run([export_command], shell=True)

        # turns feature table into a human-reable format
        biom_command = 'biom convert -i collapse-'+item+\
        '-frequency/feature-table.biom -o otu-'+item+\
        '-table.tsv --to-tsv --header-key taxonomy'

        result = subprocess.run([biom_command], shell=True)

        # formatting the table so that it is in the correct order
        table = pd.read_table(\"otu-"+str(item)+"-table.tsv\", sep='\t', header=1)
        table = table.drop(columns=['taxonomy'])
        table = table.rename(columns={'#OTU ID':'taxonomy'})
        tax = table.pop('taxonomy')
        insertion_site = len(table.columns)
        table.insert(insertion_site, 'taxonomy', tax)
        table.insert(0, 'OTU_ID', np.arange(len(table)))
        table.to_csv('otu-'+str(item)+'-mod-table.tsv', sep='\t', index=False)

        # human readable table into compressed computer-readble format
        biom_format_command='biom convert -i otu-'+str(item)+ \
        '-mod-table.tsv -o otu-table-mod.biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy'

        result = subprocess.run([biom_format_command], shell=True)

        # Outputs the current ioi so that it can be annotatted in the graphlan image
        with open('current.txt', 'w') as file:
            file.write(item)

        # bash script call to handle the steps within a conda python 2.7.17 envionment
        generate_image_command = 'bash graph.sh'

        result = subprocess.run([generate_image_command], shell=True)

        # renaming otu tables so they have meaning
        rename_table = 'cp otu-table-mod.biom phylo_trees/otu-table-'+item+'-mod.biom'
        
        result = subprocess.run([rename_table],shell=True)

        # renaming the output of the graping bash script so that it has meaning
        rename_image = 'cp image_graph.png phylo_trees/image_'+item+'_graph.png'

        result = subprocess.run([rename_image], shell=True)

        # rename pdf quality image so that it has meaning
        rename_pdf_image = 'cp image_pdf_graph.png phylo_trees/image_'+item+'_pdf_g.png'

        result = subprocess.run([rename_pdf_image], shell=True)
        
    """

}

process LefseFormat {
    publishDir "${params.outdir}/lefse", mode: 'copy'

    //conda "${projectDir}/r_env.yml"
    //conda "r_env.yml"
    container "docker://lorentzb/r_latest"

    input:
    val ioi from ch_ioi_lefse
    file "table-dada2.qza" from ch_table_lefse
    file "rooted-tree.qza" from ch_tree_lefse
    file "taxonomy.qza" from ch_tax_lefse
    file metadata from ch_metadata_lefse
    file "qiime_to_lefse.R" from ch_lefse_format_script
    file "set.txt" from ch_r_wait
    file "init_and_refresh.R" from ch_r_init
    file "renv.lock" from ch_r_lock
    

    output:
    path "combos/*" into ch_paired_lefse_format
    file "table-dada2.qza" into ch_table_report
    file "rooted-tree.qza" into ch_tree_report
    file "taxonomy.qza" into ch_tax_report
    file "metadata.tsv" into ch_metadata_report


    script:
    """
    #!/usr/bin/env bash
    mkdir combos
    Rscript init_and_refresh.R
    cp ${metadata} "metadata.tsv"
    Rscript qiime_to_lefse.R ${ioi}
    mv lefse_formatted.txt combos/
    """
}

process LefseAnalysis{
    publishDir "${params.outdir}/lefse", mode: 'copy'

    //conda "${projectDir}/python2_env.yml"
    //conda "python2_env.yml"
    container "docker://lorentzb/py2_env"

    input:
    path "combos/*" from ch_paired_lefse_format
    file "lefse_analysis.sh" from ch_lefse_analysis_script
    file plot_clado from ch_clado_file
    file plot_res from ch_plot_res

    output:
    path "result/*" into ch_lefse_results

    script:
    """
    #!/usr/bin/env bash
    mkdir result
    bash lefse_analysis.sh
    """
}

process ExportSetup{
    publishDir "${params.outdir}", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    //conda "environment.yml"
    container "docker://lorentzb/automate_16_nf"

    input:
    file "stats-dada2.qzv" from ch_dada_stats_export
    file metadata from ch_metadata_finalize 

    output:
    file "dada2_stats.tsv" into ch_dada_stats_file
    file "metadata.tsv" into ch_metadata_renamed

    script:
    """
    #!/usr/bin/env bash

    qiime tools export \
    --input-path stats-dada2.qzv \
    --output-path stats-dada2

    cp stats-dada2/metadata.tsv dada2_stats.tsv

    cp ${metadata} metadata.tsv  
    """
}


process GenerateReport{
     publishDir "${baseDir}", mode: 'move'

    //conda "${projectDir}/r_env.yml"
    //conda "r_env.yml"
    container "docker://lorentzb/r_latest"
    label 'r'

    input:
    file "item_of_interest.csv" from ch_ioi_file_out
    file "table-dada2.qza" from ch_table_report
    file "rooted-tree.qza" from ch_tree_report
    file "taxonomy.qza" from ch_tax_report
    file metadata from ch_metadata_report
    path "phylo_trees/*" from ch_png_phylo_tree
    path "shannon/*" from ch_shannon_path
    path "simpson/*" from ch_simpson_path
    path "chao1/*" from ch_chao_path
    path "ace/*"  from ch_ace_path
    path "obs/*" from ch_obs_path
    path "faith_pd/*" from ch_faith_path
    path "core-metric-results/*" from ch_core_report
    path "alpha-rareplot/*" from ch_alpha_rare_viz
    path "unweighted-sig/*" from ch_u_unifrac_beta_path
    path "weighted-sig/*" from ch_w_unifrac_beta_path
    path "result/*" from ch_lefse_results
    file "report.Rmd" from ch_report_outline
    file "make_report.sh" from ch_report_bash_script
    file "order_item_of_interest.csv" from ch_format_ioi_order



    output:
    file "done.txt" into ch_done
       
    shell:
    '''
    #!/usr/bin/env bash
    echo "all files copied!"
    echo ":)" > done.txt
    cd !{baseDir}
    echo $PWD
    echo ' !{baseDir}/!{params.outdir} ' > out.txt
    '''
}


