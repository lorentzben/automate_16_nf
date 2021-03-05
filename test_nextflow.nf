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
    file 'manifest_format.txt' into manifest_type
    file 'data_type.txt' into dataType

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import os 
    read_manifest = pd.read_table(${manifest}, index_col=0, sep='\t')

    if read_manifest.columns[0] == 'absolute-filepath':
        print("single end analysis")
        format = "SingleEndFastqManifestPhred33V2"
        data = "SampleData[SequencesWithQuality]"
        with open("manifest_format.txt", "w") as file:
            file.write(format)

        with open("data_type.txt", "w") as d_file:
            d_file.write(data)

        
    elif read_manifest.columns[0] == 'forward-absolute-filepath':
        print("paired end analsis")
        format = "PairedEndFastqManifestPhred33V2"
        data = "SampleData[PairedEndSequencesWithQuality]"
        with open("manifest_format.txt", "w") as file:
            file.write(format)

         with open("data_type.txt", "w") as file:
            d_file.write(data)
    else:
        print(
            "cannot determine if paired or single end, check manifest file")
        exit(1)

    print(format + " " + data)
    """
    
}