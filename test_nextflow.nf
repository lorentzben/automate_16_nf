#!/usr/bin/env nextflow

params.manifest = "$baseDir/c5_litter_mapping.tsv"
params.sequence = "$baseDir/seqs"
params.outdir = "$baseDir/C5_litter"

if(params.manifest) {
    tsvFile = file(params.manifest).getName()
    Channel
        .fromPath(params.manifest)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .into{ ch_single_pair ; ch_make_qiime }
}

if(params.sequence){ 
    Channel
        .fromPath(params.sequence)
        .ifEmpty {exit 1, log.info "Cannot find path file ${sequences}"}
        .into{ ch_make_qiime_seq }

}

process check_single_paired { 
    
    input: 
    file manifest from ch_single_pair

    output: 
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType
    
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
result.subscribe { println it }

process generate_seq_object{

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
