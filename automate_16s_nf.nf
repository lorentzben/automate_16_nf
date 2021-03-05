#!/usr/bin/env nextflow

// This script takes a lot of inspiration from nf-core/ampliseq 

def helpMessage(){
    log.info nfcoreHeader()
	log.info"""
	Usage:

    This script should be called as follows:
    nextflow run automate_16_nf/analysis_16 --input seqs --metadata metadata.tsv --manifest manifest.tsv --item_of_interest "Day"

    Main arguments:
        --input [path/to/folder]      Folder containing demultiplexed fastq.gz files
        --metadata [path/to/file]     Path to metadata sheet in tsv format see EXAMPLE_METADATA.tsv
        --manifest [path/to/file]     Path to mapping file in tsv format see EXAMPLE_MAPPING.tsv 
        --item_of_interest [str]      Item of interest, group defining treatment vs control or longitudinal variable
        --email [email]               Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
	    --email_on_fail [email]       Same as --email, except only send mail if the workflow is not successful
        -name [str]                   Name for the analysis run, if not provided nextflow will generate one 
        --outdir [file]               The output directory where the results will be saved 
    """.stripIndent()
}

//show help message 
if (params.help){
    helpMessage()
    exit 0 
}


/*
 * Import input files
 */

if (params.metadata) {
	Channel.fromPath("${params.metadata}", checkIfExists: true)
		.into {ch_metadata_for_test ; ch_metadata_for_test_2}
if(!params.input){
    exit 1, "Option --input missing"
}

// check the user provided name 
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

//logger info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name'] = 'automate 16s nf'
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Input'] = params.manifest ?: params.input
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job" 
summary['Output dir'] = params.outdir
summary['Launch dir'] = params.workDir
summary['Script dir'] = params.projectDir
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

//Channel for manifest file

if(params.manifest) {
    tsvFile = file(params.manifest).getName()
    Channel
        .fromPath(params.manifest)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .into{ ch_single_pair ; ch_make_qiime }

}

//Channel for metadata file
if(params.metadata) {  
    Channel.fromPath("${params.metadata}", checkIfExists: true)
    .into { ch_feature_visualization }
}

// Determine Single or Paired analysis 

process check_single_paired { 
    
    input: 
    file manifest from ch_single_pair

    output: 
    value manifest_format into manifest_type
    value data_type into dataType

        script:
        """
        !#/usr/bin/python3
        import pandas as pd
        read_manifest = pd.read_table(${manifest}, index_col=0, sep='\t')

        if read_manifest.columns[0] == 'absolute-filepath':
            print("single end analysis")
            $manifest_format = "single"
        elif read_manifest.columns[0] == 'forward-absolute-filepath':
            print("paired end analsis")
            $manifest_format = 'paired'
        else:
            print(
                "cannot determine if paired or single end, check manifest file")
            exit(1)

        print($manifest_format)
        """
    
}