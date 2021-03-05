#!/usr/bin/env nextflow

params.manifest = "$baseDir/c5_litter_mapping.tsv"

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
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType
    stdout into result

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
    input: 
    file manifest from ch_make_qiime
    file manifest_format from manifest_type
    file data_type from dataType
    

    output: 
    file 'demux.qza' into qiime_obj
    stdout into printer

    script:
    """
    module load  QIIME2/2020.11
    qiime tools import \
    --type ${data_type} \
    --input-path ${manifest} \
    --output-path demux.qza \
    --input-format ${manifest_format}
    """

}
printer.subscribe{ println it } 