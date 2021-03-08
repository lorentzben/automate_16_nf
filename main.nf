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
        .into{ ch_placeholder }
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

    conda 'python2_env.yml'

    shell:
    '''
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

    conda 'environment.yml'

    script:
    """
    #!/usr/bin/env python3

    import os 
    from pathlib import Path
    from pathlib import PurePath
    import pandas as pd 
    import csv 

    try:
        read_manifest = pd.read_table('${manifest}', index_col=0, sep='\t')
    except FileNotFoundError:
        log.info("that manifest file does not exist")
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
                    log.info(filename + ' found')
                    found.append(filename)
                else:
                    if filename in gz_files:
                        log.info(filename + ' found')
                        found.append(filename)
                    else:
                        log.info(filename + ' missing')
                        missing.append(filename)
        except KeyError:
            log.info('single read project')

        # try except in the case that the user only has single end reads.
        try:
            for item in read_manifest['reverse-absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    log.info(filename + ' found')
                    found.append(filename)
                else:
                    if filename in gz_files:
                        log.info(filename + ' found')
                        found.append(filename)
                    else:
                        log.info(filename + ' missing')
                        missing.append(filename)

        except KeyError:
            log.info("looking for forward only reads")
    else:
        # this case is if there are only single reads and after which we can figure that the manifest file is wrong
        try:
            for item in read_manifest['absolute-filepath']:
                filename = os.path.split(item)[1]
                if filename in fastq_files:
                    log.info(filename + ' found')
                    found.append(filename)
                else:
                    if filename in gz_files:
                        log.info(filename + ' found')
                        found.append(filename)
                    else:
                        log.info(filename + ' missing')
                        missing.append(filename)
        except KeyError:
            log.info("headings in the manifest appear to be incorrect")
            exit(1)

    # if missing is not an empty list, i.e. a file listed in the manifest is not detected, it raises an error and
    # creates a list for the user
    if missing != []:

        log.info("files are missing, please see missing.csv to correct")

        with open('missing.csv', 'w', newline='') as csvfile:
            fieldnames = ['filename']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for filename in missing:
                writer.writerow({'filename': filename})

        exit(0)

    log.info("the manifest called: " + manifest +
                 " is valid and ready to go")
    """

}

process CheckSinglePaired { 
    
    input: 
    file manifest from ch_single_pair

    output: 
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType
    
    conda 'environment.yml'

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

    input: 
    file manifest from ch_make_qiime
    file manifest_format from manifest_type
    file data_type from dataType
    path seqs from ch_make_qiime_seq

    output: 
    file 'demux.qza' into qiime_obj

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

