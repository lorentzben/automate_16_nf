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
        --name [str]                  Name for the analysis run, if not provided nextflow will generate one
        --sampDepth [str]             Value of sampling depth derived from demux_summary.qzv
        --rareDepth [str]             Value of rarefaction depth derived from table.qzv
        --forward [str]               Value of the custom forward cutoff
        --rev [str]                   Value of the custom reverse cutoff
        --classify [path/to/folder]   Folder containing 515 and whole 16s classifier qza files.
        --classify_515 [path/to/file] Path to 515-806.qza classifier.
        --classify_full [path/to/file] Path to full 16s.qza classifier. 
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
        .into{ ch_meta_veri ; ch_meta_feature_viz; ch_alpha_metadata ; ch_metadata_rare_curve ; ch_metadata_alpha_sig ; 
        ch_metadata_beta_sig ; ch_metadata_phylo_tree ; ch_metadata_phylo_tree_run ; ch_metadata_lefse ; ch_metadata_finalize}
}


if(params.sampDepth){
    Channel
        .from(params.sampDepth)
        .set{ ch_user_sample_depth }
}
if (params.forward){
    Channel
        .from(params.forward)
        .set{ ch_user_forward }
}

if(!params.forward){
    Channel
        .from(0)
        .set{ ch_user_forward }
}

if(params.rev){
    Channel
        .from(params.rev)
        .set{ ch_user_rev }
}

if(!params.rev){
    Channel
        .from(0)
        .set{ ch_user_rev }
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

if(params.classify_515){
    fiveFile = file(params.classify_515).getName()
    Channel
        .fromPath(params.classify_515)
        .ifEmpty {exit 1, log.info "Cannot find the classifier!"}
        .set{ ch_515_classifier }
}

if(params.classify_full){
    fullFile = file(params.classify_full).getName()
    Channel
        .fromPath(params.classify_full)
        .ifEmpty {exit 1, log.info "Cannot find the classifier!"}
        .set{ ch_whole_classifier}
}


Channel
    .fromPath("${baseDir}/python_scripts/plot_cladogram.py")
    .ifEmpty {exit 1, log.info "Cannot find file plot_cladogram.py!"}
    .set{ ch_clado_file }

Channel
    .fromPath("${baseDir}/python_scripts/plot_res.py")
    .ifEmpty {exit 1, log.info "Cannot find file plot_res.py!"}
    .set{ ch_plot_res }



Channel
    .from(params.itemOfInterest)
    .ifEmpty {exit 1, log.info "Cannot find Item of interest"}
    .into{ ch_ioi_veri ; ch_ioi_beta_sig ; ch_ioi_phylo_tree ; ch_ioi_phylo_tree_run ; ch_ioi_lefse ;\
    ch_ioi_denoise_to_file ; ch_ioi_r01_csv ; ch_ioi_r02_csv ; ch_ioi_r03_csv ; ch_ioi_r04_csv ; ch_ioi_r05_csv ;
    ch_ioi_r06_csv ; ch_ioi_r07_csv ; ch_ioi_r08_csv; ch_ioi_r09_csv; ch_ioi_r10_csv; ch_ioi_r11_csv; ch_ioi_r12_csv; ch_ioi_r13_csv}

Channel
    .fromPath("${baseDir}/bash_scripts/graph.sh")
    .set{ ch_graph_script } 

Channel
    .fromPath("${baseDir}/r_scripts/qiime_to_lefse.R")
    .set { ch_lefse_format_script }

Channel
    .fromPath("${baseDir}/python_scripts/filter_samples.py")
    .set{ ch_filter_script }

Channel 
    .fromPath("${baseDir}/bash_scripts/lefse_analysis.sh")
    .set{ ch_lefse_analysis_script }

Channel
    .fromPath("${baseDir}/report_gen_files/01_report.Rmd")
    .set{ ch_01_report_file }

Channel
    .fromPath("${baseDir}/report_gen_files/02_report.Rmd")
    .set{ ch_02_report_file }

Channel
    .fromPath("${baseDir}/report_gen_files/03_report.Rmd")
    .set{ ch_03_report_file }

Channel
    .fromPath("${baseDir}/report_gen_files/04_report.Rmd")
    .set{ ch_04_report_file }

Channel
    .fromPath("${baseDir}/report_gen_files/05_report.Rmd")
    .set{ ch_05_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/06_report.Rmd")
    .set{ ch_06_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/07_report.Rmd")
    .set{ ch_07_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/08_report.Rmd")
    .set{ ch_08_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/09_report.Rmd")
    .set{ ch_09_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/10_report.Rmd")
    .set{ ch_10_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/11_report.Rmd")
    .set{ ch_11_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/12_report.Rmd")
    .set{ ch_12_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/13_report.Rmd")
    .set{ ch_13_report_file }

Channel 
    .fromPath("${baseDir}/report_gen_files/14_report.Rmd")
    .set{ ch_14_report_file }


process VerifyManifest{
    publishDir "${params.outdir}", mode: 'copy'
    input:
    file manifest from ch_mani_veri
    path seqs_dir from ch_seqs_veri
    file metadata from ch_meta_veri
    val ioi from ch_ioi_veri

    output:

    file "order_item_of_interest.csv" into ( ch_format_ioi_order, ch_oioi_r01_csv, ch_oioi_r02_csv, ch_oioi_r03_csv, ch_oioi_r04_csv,  ch_oioi_r05_csv ,
    ch_oioi_r06_csv, ch_oioi_r07_csv, ch_oioi_r08_csv, ch_oioi_r09_csv, ch_oioi_r10_csv, ch_oioi_r11_csv, ch_oioi_r12_csv, ch_oioi_r13_csv )

    label 'process_low'
    container "docker://lorentzb/automate_16_nf:2.0"

    script:
    """
    #!/usr/bin/env python3
    import os 
    from pathlib import Path
    from pathlib import PurePath
    import pandas as pd 
    import csv 

    cwd = os.getcwd()
    
    try:
        metadata_path = str('${metadata}')
        metadata_path = cwd+'/'+metadata_path
        read_metadata = pd.read_table(metadata_path, index_col=0, sep='\t')
    except FileNotFoundError:
        print("Incorrect Path for Metadata!")
        print("Path provided: "+metadata_path )
        exit(1)

    try:
        read_order = pd.read_table('${baseDir}/order_item_of_interest.csv', sep=',')
        pd.DataFrame.to_csv(read_order, 'order_item_of_interest.csv', index=False)
    except FileNotFoundError:
        iois = list(pd.Series.unique(read_metadata['${ioi}']))
        ioisdf = pd.DataFrame(iois[1:])
        ioisdf.columns = ['${ioi}']
        pd.DataFrame.to_csv(ioisdf, 'order_item_of_interest.csv', index=False)

    seq_dir = '${seqs_dir}'
    try:
        manifest_path = str('${manifest}')
        manifest_path = cwd+'/'+manifest_path
        read_manifest = pd.read_table(manifest_path, index_col=0, sep='\t')
    except FileNotFoundError:
        print("Incorrect Path for manifest!")
        print("Path provided: " +manifest_path)
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

    container "docker://lorentzb/automate_16_nf:2.0"
    
    input: 
    file manifest from ch_single_pair
    val ioi from ch_ioi_denoise_to_file

    output: 
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType
    file "item_of_interest.csv" into ch_ioi_file_out

    label 'process_low'
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
    
    container "docker://lorentzb/automate_16_nf:2.0"

    input: 
    file manifest from ch_make_qiime
    file manifest_format from manifest_type
    file data_type from dataType
    path seqs from ch_make_qiime_seq

    output: 
    file 'demux.qza' into ch_qiime_obj
    file manifest_format into ch_manifest_type
    
    label 'process_medium'
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
    file 'demux_summary.qzv' into ch_demux_export

    container "docker://lorentzb/automate_16_nf:2.0"

    label 'process_medium'
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
    
    publishDir "${params.outdir}/qiime", mode: 'copy'

    input:
    file 'manifest_format.txt' from ch_manifest_type
    file ('demux_summary/*') from ch_qiime_qual
    val forward_val from ch_user_forward 
    val reverse_val from ch_user_rev
    
    output: 
    file("cutoffs.csv") into (ch_cutoff_vals, ch_cutoff_r03)
    file("manifest_format.txt") into ch_manifest_type_denoise

    container "docker://lorentzb/automate_16_nf:2.0"

    label 'process_low'

    script:
    
    if(forward_val != 0)
        """
        #!/usr/bin/env python3
        import pandas as pd 
        from pathlib import Path
        import numpy as np 
        import csv 
        
        forward = [0,${forward_val}]
        reverse =[0, ${reverse_val}]
        with open('cutoffs.csv', 'w', newline='') as csvfile:
                fieldnames = ['cutoff', 'value']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

                writer.writeheader()
                writer.writerow({'cutoff': 'forward left', 'value': forward[0]})
                writer.writerow({'cutoff': 'forward right', 'value': forward[1]})
                writer.writerow({'cutoff': 'reverse left', 'value': reverse[0]})
                writer.writerow({'cutoff': 'reverse right', 'value': reverse[1]})
                writer.writerow({'cutoff': 'filename', 'value': "user_sub"})
                writer.writerow({'cutoff': 'filename', 'value': "user_sub"})
        """
    else 
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
    
    container "docker://lorentzb/automate_16_nf:2.0"

    label 'process_medium'

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
            --p-n-threads 0 \
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
            --p-n-threads 0 \
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
    file "table.qzv" into ch_table_viz_export
    file "rep-seqs.qzv" into ch_req_seq_vis_obj
    file "rep-seqs-dada2.qza" into ch_rep_seq_tree_gen
    file "table-dada2.qza" into ch_alpha_div_table

    container "docker://lorentzb/automate_16_nf:2.0"

    label 'process_low'

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

    container "docker://lorentzb/automate_16_nf:2.0"
    
    input:
    file "rep-seqs-dada2.qza" from ch_rep_seq_tree_gen

    output:
    file "aligned-rep-seqs.qza" into ch_aligned_rep_seqs
    file "masked-aligned-rep-seqs.qza" into ch_mask_align_rep_seq
    file "unrooted-tree.qza" into ch_unrooted_tree
    file "rooted-tree.qza" into (ch_rooted_tree, ch_rooted_tree_r01 , ch_root_tree_r03, ch_root_tree_r06, 
    ch_root_tree_r08, ch_root_tree_r11, ch_root_tree_r12)
    file "rep-seqs-dada2.qza" into ch_rep_seq_classify

    label 'process_medium'
    
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

    container "docker://lorentzb/automate_16_nf:2.0"
    
    input:
    file "table.qzv" from ch_table_viz_obj

    output:
    path "table_viz/*" into ch_table_viz_dir
    path "table_viz/*" into ch_table_viz_dir_rare

    label 'process_low'

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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    path "table_viz/*" from ch_table_viz_dir

    output:
    file "sampling_depth.csv" into (ch_sampling_depth_csv, ch_samp_depth_r03)
    file "samp_depth_simple.txt" into ch_depth

    label 'process_low'

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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file metadata from ch_alpha_metadata
    file "table-dada2.qza" from ch_alpha_div_table
    file "rooted-tree.qza" from ch_rooted_tree
    file "samp_depth_simple.txt" from ch_depth
    val user_depth from ch_user_sample_depth

    output:
    path "core-metric-results/*" into ch_core_beta_significance 
    path "core-metric-results/*" into ( ch_core_report , ch_rare_table_r01 , ch_core_metric_r03, ch_core_metric_r06, ch_core_metric_r08, 
    ch_core_metric_r09, ch_core_metric_r11, ch_core_metric_r12 )
    file "core-metric-results/rarefied_table.qza" into ch_phylo_tree_rare_table_run
    file "shannon.qza" into ch_shannon_qza
    file "simpson.qza" into ch_simpson_qza 
    file "chao1.qza" into ch_chao_qza
    file "ace.qza" into ch_ace_qza
    file "obs.qza" into ch_obs_qza
    file "faith_pd.qza" into ch_faith_qza
    file "table-dada2.qza" into ch_table_rare_curve
    file "rooted-tree.qza" into ch_tree_rare_curve
    

    label 'process_medium'

    shell:
    '''
    #!/usr/bin/env bash

    echo !{user_depth} > test_samp_depth.txt

    if [[ !{user_depth} == 0 ]]
    then
        SAMP_DEPTH=$(head samp_depth_simple.txt)
    else
        SAMP_DEPTH=!{user_depth}
    fi

    echo !{user_depth} > used_samp_depth.txt
  
    qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table-dada2.qza \
    --p-sampling-depth $SAMP_DEPTH \
    --m-metadata-file !{metadata} \
    --output-dir core-metric-results 

    qiime diversity alpha \
    --i-table core-metric-results/rarefied_table.qza \
    --p-metric shannon \
    --o-alpha-diversity shannon.qza

    qiime diversity alpha \
    --i-table core-metric-results/rarefied_table.qza \
    --p-metric simpson \
    --o-alpha-diversity simpson.qza 

    qiime diversity alpha \
    --i-table core-metric-results/rarefied_table.qza \
    --p-metric chao1 \
    --o-alpha-diversity chao1.qza

    qiime diversity alpha \
    --i-table core-metric-results/rarefied_table.qza \
    --p-metric ace \
    --o-alpha-diversity ace.qza

    qiime diversity alpha \
    --i-table core-metric-results/rarefied_table.qza \
    --p-metric observed_features \
    --o-alpha-diversity obs.qza 

    qiime diversity alpha-phylogenetic \
    --i-table core-metric-results/rarefied_table.qza \
    --i-phylogeny rooted-tree.qza \
    --p-metric faith_pd \
    --o-alpha-diversity faith_pd.qza 
    '''
}

process AssignTaxonomy{
    //TODO change out the classifier for the 515 only one
    //TODO check a different classifier out
    publishDir "${params.outdir}/qiime", mode: 'copy'

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file "rep-seqs-dada2.qza" from ch_rep_seq_classify
    file "16s-whole-seq-classifier.qza" from ch_whole_classifier
    file "515-806-classifier.qza" from ch_515_classifier

    output:
    file "taxonomy.qza" into ( ch_taxonomy_phylo_tree, ch_taxonomy_r01, ch_taxonomy_r03, 
    ch_taxonomy_r06, ch_taxonomy_r08, ch_taxonomy_r11, ch_taxonomy_r12)
    file "taxonomy.qza" into ch_taxonomy_phylo_tree_run
    file "taxonomy.qzv" into ch_classified_qzv
    
    label 'process_medium'

    script:
    """
    #!/usr/bin/env bash

    if [ ! -f "16s-whole-seq-classifier.qza" ]
    then 
        echo "Error, download the classifier from readme"
        exit 1
    fi
    if [ ! -f "515-806-classifier.qza" ]
    then 
        echo "Error, download the classifier from readme"
        exit 1
    fi

    qiime feature-classifier classify-sklearn \
    --i-classifier 515-806-classifier.qza \
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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    path "table_viz/*" from ch_table_viz_dir_rare

    output:
    file "rare_depth.txt" into ch_rare_curve_depth

    label 'process_low'

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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file "rare_depth.txt" from ch_rare_curve_depth
    file metadata from ch_metadata_rare_curve
    file "table-dada2.qza" from ch_table_rare_curve
    file "rooted-tree.qza" from ch_tree_rare_curve
    val user_rare_depth from ch_user_rarefaction_depth


    output:
    file "alpha-rarefaction.qzv" into ch_alpha_rare_obj
    path "alpha-rareplot/*" into (ch_alpha_rare_viz, ch_alpha_rare_r07)
    file "table-dada2.qza" into ch_table_phylo_tree_rare
    file "rooted-tree.qza" into ch_tree_lefse
    
    label 'process_medium'

    shell:
    '''
    #!/usr/bin/env bash

    if [[ !{user_rare_depth} == 0 ]]
    then
        DEPTH=$(head rare_depth.txt)
    else
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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file metadata from ch_metadata_alpha_sig
    file "shannon.qza" from ch_shannon_qza
    file "simpson.qza" from ch_simpson_qza
    file "chao1.qza" from ch_chao_qza
    file "ace.qza" from ch_ace_qza
    file "obs.qza" from ch_obs_qza
    file "faith_pd.qza" from ch_faith_qza

    output:
    path "shannon/*" into ( ch_shannon_path, ch_shannon_r04, ch_shannon_r05 )
    path "simpson/*" into ( ch_simpson_path, ch_simpson_r04, ch_simpson_r05 )
    path "chao1/*" into ( ch_chao_path, ch_chao_r04, ch_chao_r05 )
    path "ace/*" into ( ch_ace_path, ch_ace_r04, ch_ace_r05 ) 
    path "obs/*" into ( ch_obs_path, ch_obs_r04, ch_obs_r_05 ) 
    path "faith_pd/*" into ( ch_faith_path, ch_faith_r04, ch_faith_r05 ) 

    label 'process_medium'

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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    val ioi from ch_ioi_beta_sig
    file metadata from ch_metadata_beta_sig
    path "core-metric-results/*" from ch_core_beta_significance 

    output:
    path "unweighted-sig/*" into ( ch_u_unifrac_beta_path, ch_u_uni_r10 )
    path "weighted-sig/*" into ( ch_w_unifrac_beta_path, ch_w_uni_r10 )

    label 'process_medium'
    
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

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file metadata from ch_metadata_phylo_tree
    val ioi from ch_ioi_phylo_tree
    file "table-dada2.qza" from ch_table_phylo_tree_rare
    file "taxonomy.qza" from ch_taxonomy_phylo_tree
    file "filter_samples.py" from ch_filter_script

    output:
    
    file "table-dada2.qza" into ch_table_graphlan2
    file "taxonomy.qza" into ch_tax_lefse
    path "biom_tabs/*" into ch_biom_tabs

    label 'process_medium'

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
    subprocess.run(['mkdir biom_tabs'], shell=True)

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:
        item = str(item)
        # filters/splits the feature table based on the current ioi
        
        filter_command = "python3 filter_samples.py -m ${metadata} -i ${ioi} -c "+str(item)
        result = subprocess.run([filter_command], shell=True)

        time.sleep(2)

        # adds taxonomic info needed for plotting
        collapse_command = 'qiime taxa collapse \
        --i-table '+str(item)+'-filtered-table.qza \
        --o-collapsed-table collapse-'+str(item)+'-table.qza \
        --p-level 7 \
        --i-taxonomy taxonomy.qza'

        result = subprocess.run([collapse_command], shell=True)

        # exports artifact so that the next step can collect it
        export_command='qiime tools export \
        --input-path collapse-'+str(item)+'-table.qza \
        --output-path collapse-'+str(item)+'-frequency/'
        
        result = subprocess.run([export_command], shell=True)

        # turns feature table into a human-reable format
        biom_command = 'biom convert -i collapse-'+str(item)+\
        '-frequency/feature-table.biom -o otu-'+str(item)+\
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
        '-mod-table.tsv -o '+str(item)+'-otu-table-mod.biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy'

        result = subprocess.run([biom_format_command], shell=True)

        result = subprocess.run(['cp '+str(item)+'-otu-table-mod.biom biom_tabs'],shell=True)
    """

}

process runGraphlan{
    publishDir "${params.outdir}/graphlan", mode: 'copy'

    container "docker://lorentzb/py2_test:2.0"

    input:
    file metadata from ch_metadata_phylo_tree_run
    val ioi from ch_ioi_phylo_tree_run
    file "table-dada2.qza" from ch_table_graphlan2
    file "rarefied_table.qza" from ch_phylo_tree_rare_table_run
    file "taxonomy.qza" from ch_taxonomy_phylo_tree_run
    file "graph.sh" from ch_graph_script
    path "biom_tabs/*" from ch_biom_tabs
    
    output:
    path "phylo_trees/*" into (ch_png_phylo_tree,  ch_02_report_imgs) 
    file "table-dada2.qza" into ch_table_lefse_graphlan
    file "rarefied_table.qza" into ch_table_lefse
    
    label 'process_low'

    script:
    """
    #!/usr/bin/env python2
    import subprocess
    import csv
    import pandas as pd
    import numpy as np 
    import time
    import os

    metadata_table= pd.read_table(\"${metadata}\", sep='\t')
    metadata_table = metadata_table.drop([0,1])

    ioi_set = set(metadata_table[\"${ioi}\"])
    ioi = '${ioi}'

    os.system('cp biom_tabs/*-otu-table-mod.biom .')

    os.system('mkdir phylo_trees')

    # iterates over the items of interest to produce a circular phylogenetic tree per category e.g. CONTROL TREATMENT
    for item in ioi_set:
        item = str(item)
        # filters/splits the feature table based on the current ioi

        # Outputs the current ioi so that it can be annotated in the graphlan image
        with open('current.txt', 'w') as file:
            file.write(item)

        # bash script call to handle the steps within a conda python 2.7.17 envionment
        generate_image_command = 'bash graph.sh'
        result = os.system(generate_image_command)

    rename_image = 'cp *_image_graph.png phylo_trees/.'
    result = os.system(rename_image)

    rename_pdf_image = 'cp *_image_pdf_graph.png phylo_trees/.'
    result = os.system(rename_pdf_image)
    """

}

process LefseFormat {
    publishDir "${params.outdir}/lefse", mode: 'copy'

    container "docker://lorentzb/qiime2lefse:1.0"

    input:
    val ioi from ch_ioi_lefse
    file "table-dada2.qza" from ch_table_lefse_graphlan
    file "rarefied_table.qza" from ch_table_lefse
    file "rooted-tree.qza" from ch_tree_lefse
    file "taxonomy.qza" from ch_tax_lefse
    file metadata from ch_metadata_lefse
    file "qiime_to_lefse.R" from ch_lefse_format_script
      

    output:
    path "combos/*" into ch_paired_lefse_format
    file "table-dada2.qza" into ch_table_report_raw
    file "rarefied_table.qza" into ch_table_report_rare
    file "rooted-tree.qza" into ch_tree_report
    file "taxonomy.qza" into ch_tax_report
    file "metadata.tsv" into ( ch_metadata_report, ch_metadata_r01, ch_metadata_r02, ch_metadata_r03, 
    ch_metadata_r04, ch_metadata_r05, ch_metadata_r06, ch_metadata_r07, ch_metadata_r08, ch_metadata_r09, 
    ch_metadata_r10, ch_metadata_r11, ch_metadata_r12, ch_metadata_r13 )

    label 'process_medium'

    script:
    """
    #!/usr/bin/env bash
    mkdir combos
    cp ${metadata} "metadata.tsv"
    Rscript qiime_to_lefse.R ${ioi}
    mv lefse_formatted.txt combos/
    """
}

process LefseAnalysis{
    publishDir "${params.outdir}/lefse", mode: 'copy'

    container "docker://lorentzb/py2_env:2.0"

    input:
    path "combos/*" from ch_paired_lefse_format
    file "lefse_analysis.sh" from ch_lefse_analysis_script
    file plot_clado from ch_clado_file
    file plot_res from ch_plot_res

    output:
    path "result/*" into ( ch_lefse_results, ch_lefse_r13 ) 

    label 'process_medium'

    script:
    """
    #!/usr/bin/env bash
    mkdir result
    bash lefse_analysis.sh
    """
}

process ExportSetup{
    publishDir "${params.outdir}", mode: 'copy'

    container "docker://lorentzb/automate_16_nf:2.0"

    input:
    file "stats-dada2.qzv" from ch_dada_stats_export
    file metadata from ch_metadata_finalize 

    output:
    file "dada2_stats.tsv" into (ch_dada_stats_file, ch_dada_r03)
    file "metadata.tsv" into ch_metadata_renamed

    label 'process_medium'

    script:
    """
    #!/usr/bin/env bash

    qiime tools export \
    --input-path stats-dada2.qzv \
    --output-path stats-dada2

    cp stats-dada2/metadata.tsv ./dada2_stats.tsv

    cp ${metadata} ./metadata.tsv  
    """
}

process Report01 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_01:2.0"

    input:
    file "01_report.Rmd" from ch_01_report_file
    file "item_of_interest.csv" from ch_ioi_r01_csv
    file "order_item_of_interest.csv" from ch_oioi_r01_csv
    file "core-metric-results/*" from ch_rare_table_r01 
    file "rooted-tree.qza" from ch_rooted_tree_r01  
    file "taxonomy.qza" from ch_taxonomy_r01  
    file "metadata.tsv" from ch_metadata_r01

    output:
    path "01_report_*" into ch_01_reports
    path "Figures/*" into ch_01_figures

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    #Input: item_of_interest.csv order_item_of_interest.csv qiime/* metadata.tsv
    #Output: ../Figures/* 01_report_$dt.html 01_report_$dt.pdf

    echo "I am Here:"
    pwd
    ls
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('01_report.Rmd', output_file='$PWD/01_report_$dt.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"

    Rscript -e "rmarkdown::render('01_report.Rmd', output_file='$PWD/01_report_$dt.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
    '''


}

process Report02{
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_02:2.0"

    input:
    file "02_report.Rmd" from ch_02_report_file
    file "item_of_interest.csv" from ch_ioi_r02_csv
    file "order_item_of_interest.csv" from ch_oioi_r02_csv
    file "metadata.tsv" from ch_metadata_r02
    path "phylo_trees/*" from ch_02_report_imgs

    output:
    path "02_report_*" into ch_02_reports
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    #Rscript -e "rmarkdown::render('02_report.Rmd', output_file='$PWD/02_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''

}

process Report03{

    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_03:2.0"

    input:
    file "03_report.Rmd" from ch_03_report_file
    file "item_of_interest.csv" from ch_ioi_r03_csv
    file "order_item_of_interest.csv" from ch_oioi_r03_csv
    file "metadata.tsv" from ch_metadata_r03

    path "core-metric-results/*" from ch_core_metric_r03
    file "rooted-tree.qza" from ch_root_tree_r03
    file "taxonomy.qza" from ch_taxonomy_r03
    file "cutoffs.csv" from ch_cutoff_r03
    file "dada2_stats.tsv" from ch_dada_r03
    file "sampling_depth.csv" from ch_samp_depth_r03
    

    output:
    path "03_report_*" into ch_03_reports
    path "Figures/*" into ch_03_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('03_report.Rmd', output_file='$PWD/03_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('03_report.Rmd', output_file='$PWD/03_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''


}

process Report04{
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_04:2.0"

    input:
    file "04_report.Rmd" from ch_04_report_file
    file "item_of_interest.csv" from ch_ioi_r04_csv
    file "order_item_of_interest.csv" from ch_oioi_r04_csv
    file "metadata.tsv" from ch_metadata_r04

    path "shannon/*" from ch_shannon_r04
    path "simpson/*" from ch_simpson_r04
    path "chao1/*" from ch_chao_r04
    path "ace/*" from ch_ace_r04
    path "obs/*" from ch_obs_r04
    path "faith_pd/*" from ch_faith_r04
 
    output:
    path "04_report_*" into ch_04_reports
  

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    
    ls
    

    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('04_report.Rmd', output_file='$PWD/04_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('04_report.Rmd', output_file='$PWD/04_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''


}

process Report05{
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_05:2.0"

    input:
    file "05_report.Rmd" from ch_05_report_file
    file "item_of_interest.csv" from ch_ioi_r05_csv
    file "order_item_of_interest.csv" from ch_oioi_r05_csv
    file "metadata.tsv" from ch_metadata_r05

    path "obs/*" from ch_obs_r_05
    path "shannon/*" from ch_shannon_r05
    path "simpson/*" from ch_simpson_r05
    path "chao1/*" from ch_chao_r05
    path "ace/*" from ch_ace_r05
    path "faith_pd/*" from ch_faith_r05

    
    output:
    path "05_report_*" into ch_05_reports
    path "Figures/*" into ch_05_figures
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('05_report.Rmd', output_file='$PWD/05_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('05_report.Rmd', output_file='$PWD/05_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''

}

process Report06{
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_06:2.0"

    input:
    file "06_report.Rmd" from ch_06_report_file
    file "item_of_interest.csv" from ch_ioi_r06_csv
    file "order_item_of_interest.csv" from ch_oioi_r06_csv
    file "metadata.tsv" from ch_metadata_r06

    path "core-metric-results/*" from ch_core_metric_r06
    file "rooted-tree.qza" from ch_root_tree_r06
    file "taxonomy.qza" from ch_taxonomy_r06
        
    output:
    path "06_report_*" into ch_06_reports
    path "Figures/*" into ch_06_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('06_report.Rmd', output_file='$PWD/06_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('06_report.Rmd', output_file='$PWD/06_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''

}

process Report07{
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_07:2.0"

    input:
    file "07_report.Rmd" from ch_07_report_file
    file "item_of_interest.csv" from ch_ioi_r07_csv
    file "order_item_of_interest.csv" from ch_oioi_r07_csv
    file "metadata.tsv" from ch_metadata_r07

    path "alpha-rareplot/*" from ch_alpha_rare_r07

    
        
    output:
    path "07_report_*" into ch_07_reports
    path "Figures/*" into ch_07_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('07_report.Rmd', output_file='$PWD/07_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('07_report.Rmd', output_file='$PWD/07_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''


}

process Report08 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_08:2.0"

    input:
    file "08_report.Rmd" from ch_08_report_file
    file "item_of_interest.csv" from ch_ioi_r08_csv
    file "order_item_of_interest.csv" from ch_oioi_r08_csv
    file "metadata.tsv" from ch_metadata_r08

    path "core-metric-results/*" from ch_core_metric_r08
    file "rooted-tree.qza" from ch_root_tree_r08
    file "taxonomy.qza" from ch_taxonomy_r08
        
    output:
    path "08_report_*" into ch_08_reports
    path "Figures/*" into ch_08_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('08_report.Rmd', output_file='$PWD/08_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('08_report.Rmd', output_file='$PWD/08_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''
}

process Report09 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_09:2.0"

    input:
    file "09_report.Rmd" from ch_09_report_file
    file "item_of_interest.csv" from ch_ioi_r09_csv
    file "order_item_of_interest.csv" from ch_oioi_r09_csv
    file "metadata.tsv" from ch_metadata_r09

    path "core-metric-results/*" from ch_core_metric_r09
        
    output:
    path "09_report_*" into ch_09_reports
    path "Figures/*" into ch_09_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
   

    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('09_report.Rmd', output_file='$PWD/09_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('09_report.Rmd', output_file='$PWD/09_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''

}

process Report10 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_10:2.0"

    input:
    file "10_report.Rmd" from ch_10_report_file
    file "item_of_interest.csv" from ch_ioi_r10_csv
    file "order_item_of_interest.csv" from ch_oioi_r10_csv
    file "metadata.tsv" from ch_metadata_r10

    path "unweighted-sig/*" from ch_u_uni_r10 
    path "weighted-sig/*" from ch_w_uni_r10
        
    output:
    path "10_report_*" into ch_10_reports
    path "Figures/*" into ch_10_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls

    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('10_report.Rmd', output_file='$PWD/10_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('10_report.Rmd', output_file='$PWD/10_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''
}

process Report11 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_11:2.0"

    input:
    file "11_report.Rmd" from ch_11_report_file
    file "item_of_interest.csv" from ch_ioi_r11_csv
    file "order_item_of_interest.csv" from ch_oioi_r11_csv
    file "metadata.tsv" from ch_metadata_r11

    path "core-metric-results/*" from ch_core_metric_r11
    file "rooted-tree.qza" from ch_root_tree_r11
    file "taxonomy.qza" from ch_taxonomy_r11
        
    output:
    path "11_report_*" into ch_11_reports
    path "Figures/*" into ch_11_figures
    
    

    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
    
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('11_report.Rmd', output_file='$PWD/11_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('11_report.Rmd', output_file='$PWD/11_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''
}

process Report12 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_12:2.0"

    input:
    file "12_report.Rmd" from ch_12_report_file
    file "item_of_interest.csv" from ch_ioi_r12_csv
    file "order_item_of_interest.csv" from ch_oioi_r12_csv
    file "metadata.tsv" from ch_metadata_r12

    path "core-metric-results/*" from ch_core_metric_r12
    file "rooted-tree.qza" from ch_root_tree_r12
    file "taxonomy.qza" from ch_taxonomy_r12
        
    output:
    path "12_report_*" into ch_12_reports
    path "Figures/*" into ch_12_figures
    
    
    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls

    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('12_report.Rmd', output_file='$PWD/12_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('12_report.Rmd', output_file='$PWD/12_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''
}

process Report13 {
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_13:2.0"

    input:
    file "13_report.Rmd" from ch_13_report_file
    file "item_of_interest.csv" from ch_ioi_r13_csv
    file "order_item_of_interest.csv" from ch_oioi_r13_csv
    file "metadata.tsv" from ch_metadata_r13

    path "result/*" from ch_lefse_r13

    output:
    path "13_report_*" into ch_13_reports
    
    
    
    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('13_report.Rmd', output_file='$PWD/13_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    
    #Rscript -e "rmarkdown::render('13_report.Rmd', output_file='$PWD/13_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''


}

process Report14{
    publishDir "${params.outdir}/reports", mode: 'move'

    container "docker://lorentzb/r_14:2.0"

    input:
    file "14_report.Rmd" from ch_14_report_file

    output:
    path "14_report_*" into ch_14_reports
       
    label 'process_medium'
    script:
    '''
    #! /usr/bin/env bash

    ls
   
    dt=$(date '+%d-%m-%Y_%H.%M.%S');

    Rscript -e "rmarkdown::render('14_report.Rmd', output_file='$PWD/14_report_$dt.html', output_format='html_document', clean=TRUE,knit_root_dir='$PWD',intermediates_dir ='$PWD')"

    Rscript -e "rmarkdown::render('14_report.Rmd', output_file='$PWD/14_report_$dt.pdf', output_format='pdf_document', clean=TRUE,knit_root_dir='$PWD', intermediates_dir ='$PWD')"
    '''

}

/*
TODO Determine if we still want to publish this dir so we have the intermediate files for diagnositcs
process GenerateReport{
    publishDir "${baseDir}", mode: 'move'

    //conda "${projectDir}/r_env.yml"
    //conda "r_env.yml"
    container "docker://lorentzb/r_latest:2.0"
    //label 'r'

    input:
    file "item_of_interest.csv" from ch_ioi_file_out
    file "table-dada2.qza" from ch_table_report_raw
    file "rarefied_table.qza" from ch_table_report_rare
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
    //file "report.Rmd" from ch_report_outline
    //file "make_report.sh" from ch_report_bash_script
    file "order_item_of_interest.csv" from ch_format_ioi_order
    file "table.qzv" from ch_table_viz_export
    file 'demux_summary.qzv' from ch_demux_export



    output:
    file "done.txt" into ch_done

    label 'process_medium'

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
*/

