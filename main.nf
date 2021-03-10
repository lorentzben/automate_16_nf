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
        --itemOfInterest [str]      Item of interest, group defining treatment vs control or longitudinal variable
        -name [str]                   Name for the analysis run, if not provided nextflow will generate one 
        --outdir [file]               The output directory where the results will be saved 
        --email [email]               Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
	    --email_on_fail [email]       Same as --email, except only send mail if the workflow is not successful

        

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
        .into{ ch_meta_feature_viz; ch_alpha_metadata }
}

Channel
    .fromPath("${baseDir}/plot_cladogram.py")
    .set{ ch_clado_file }

Channel
    .fromPath("${baseDir}/plot_res.py")
    .set{ ch_plot_res }

process SetupPy2CondaEnv{

    input:
    file plot_clado from ch_clado_file
    file plot_res from ch_plot_res
    //conda "${projectDir}/python2_env.yml"
    conda "python2_env.yml"

    shell:
    '''
    #!/usr/bin/env bash
    LEFSE=$(which lefse.py)
    LEFSE_DIR=${LEFSE::-8}
    cp plot_res.py $LEFSE_DIR
    cp plot_cladogram.py $LEFSE_DIR
    '''

}

process VerifyManifest{

    input:
    file manifest from ch_mani_veri
    path seqs_dir from ch_seqs_veri

    /*this is in place for local deployment, but the server does not give access to the dir for some reason
    The change is nessecary to do nextflow run -r main lorentzben/automate_16_nf
    */
    //conda "${projectDir}/environment.yml"
    conda "environment.yml"

    script:
    """
    #!/usr/bin/env python3
    import os 
    from pathlib import Path
    from pathlib import PurePath
    import pandas as pd 
    import csv 

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
    
    input: 
    file manifest from ch_single_pair

    output: 
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType

    //conda "${projectDir}/environment.yml"
    conda "environment.yml"

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import os 

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
    conda "environment.yml"

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
    module load  QIIME2/2020.11
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
    conda "environment.yml"

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
    conda "environment.yml"

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

        average_qual = np.round(mean_qual.mean(axis=1), 0)-1
        mean_qual_vals = np.array(mean_qual)[0]

        if int(average_qual) < 30:
            print(
                "The Average Quality of these sequences may be a concern would you like to continue?")
            exit(0)

        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[i] >= int(average_qual):
                left_cutoff = i+1
                break
        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
                right_cutoff = len(mean_qual_vals)-i
                break
        return(left_cutoff, right_cutoff)

    def find_rev_cutoffs(dataframe):
        mean_qual = dataframe[4:5]

        average_qual = np.round(mean_qual.mean(axis=1), 0)+2
        mean_qual_vals = np.array(mean_qual)[0]

        if int(average_qual) < 30:
            print(
                "The Average Quality of these sequences may be a concern would you like to continue?")
            exit(0)

        left_cutoff = 0
        right_cutoff = len(mean_qual_vals)-1

        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[i] >= int(average_qual):
                left_cutoff = i+1
                break

        for i in range(0, len(mean_qual_vals)):
            if mean_qual_vals[len(mean_qual_vals)-1-i] >= int(average_qual):
                right_cutoff = len(mean_qual_vals)-i
                break

        return(left_cutoff, right_cutoff)

    if seq_format == "single":
        print("determining left and right cutoffs based on qual score")

        input_file = str(wd)+"/demux_summary/forward-seven-number-summaries.csv"

        summary = pd.read_table(input_file, index_col=0, sep=',')
        left_cutoff, right_cutoff = find_cutoffs(summary)

        print("right cutoff: "+str(right_cutoff))
        print("left cutoff: " + str(left_cutoff))

        with open('cutoffs.csv', 'w', newline='') as csvfile:
            fieldnames = ['cutoff', 'value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            writer.writerow({'cutoff': 'right', 'value': right_cutoff})
            writer.writerow({'cutoff': 'left', 'value': left_cutoff})
            writer.writerow({'cutoff': 'filename', 'value': input_file})

        print(left_cutoff, right_cutoff)

    elif seq_format == "paired":
        print("determining forward and revese, left and right cutoffs based on qual score")
        forward_file = str(wd)+"/demux_summary/forward-seven-number-summaries.csv"
        fr_summary = pd.read_table(forward_file, index_col=0, sep=',')

        forward = find_cutoffs(fr_summary)

        reverse_file = str(wd)+"/demux_summary/reverse-seven-number-summaries.csv"
        rev_summary = pd.read_table(reverse_file, index_col=0, sep=',')

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
    conda "environment.yml"

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
    file "stats-dada2.qzv" into ch_dada2_stats_viz
    file "table.qzv" into ch_table_viz_obj
    file "rep-seqs.qzv" into ch_req_seq_vis_obj
    file "rep-seqs-dada2.qza" into ch_rep_seq_tree_gen
    file "table-dada2.qza" into ch_alpha_div_table

    //conda "${projectDir}/environment.yml"
    conda "environment.yml"

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
    conda "environment.yml"
    
    input:
    file "rep-seqs-dada2.qza" from ch_rep_seq_tree_gen

    output:
    file "aligned-rep-seqs.qza" into ch_aligned_rep_seqs
    file "masked-aligned-rep-seqs.qza" into ch_mask_align_rep_seq
    file "unrooted-tree.qza" into ch_unrooted_tree
    file "rooted-tree.qza" into ch_rooted_tree

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
    conda "environment.yml"
    
    input:
    file "table.qzv" from ch_table_viz_obj

    output:
    path "table_viz/*" into ch_table_viz_dir

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
    conda "environment.yml"

    input:
    path "table_viz/*" from ch_table_viz_dir

    output:
    file "sampling_depth.csv" into ch_sampling_depth_csv
    val samp_depth into ch_depth_val

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
    $samp_depth=sampling_depth
    """

}

process AlphaDiversityMeasure{
    publishDir "${params.outdir}/qiime", mode: 'copy'

    //conda "${projectDir}/environment.yml"
    conda "environment.yml"

    input:
    file "sampling_depth.csv" from ch_sampling_depth_csv
    file metadata from ch_alpha_metadata
    file "table-dada2.qza" from ch_alpha_div_table
    file "rooted-tree.qza" from ch_rooted_tree
    val samp_depth from ch_depth_val

    output:
    file "core-metric-results/*" into ch_core_div_res

    script:
    """
    #!/usr/bin/env bash
    echo ${samp_depth}
    qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table-dada2.qza \
    --p-sampling-depth ${samp_depth} \
    --m-metadata-file ${metadata} \
    --output-dir core-metric-results \
    --p-n-jobs-or-threads 'auto'
    """

}