#!/usr/bin/env nextflow

params.manifest = "$baseDir/EXAMPLE_MANIFEST.tsv"

if(params.manifest) {
    tsvFile = file(params.manifest).getName()
    Channel
        .fromPath(params.manifest)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .into{ ch_single_pair ; ch_make_qiime }
}

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
        ${manifest_format} = "single"
    elif read_manifest.columns[0] == 'forward-absolute-filepath':
        print("paired end analsis")
        ${manifest_format} = 'paired'
    else:
        print(
            "cannot determine if paired or single end, check manifest file")
        exit(1)

    print(${manifest_format})
    """
    
}