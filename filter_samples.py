#!/usr/bin/env python3

import os
from os import path
import argparse
import subprocess
import pandas as pd
import numpy as np 

def filter_command(metadata_fi, item_of_int, current): 

    # filters/splits the feature table based on the current ioi
    filter_command = 'qiime feature-table filter-samples \
    --i-table table-dada2.qza \
    --m-metadata-file ' +metadata_fi+' \
    --p-where "\" '+ item_of_int +' \"=\"  '+ current +' "\"  \
    --o-filtered-table '+current+'-filtered-table.qza'

    result = subprocess.run([filter_command], shell=True)

def main(arg):
    filter_command(arg.metadata, arg.ioi, arg.current)



if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Handles the calling of the filter command for qiime because nextflow was messing with the string to call the command")
    parser.add_argument('-m', '--metadata', action='store', required=True,
                        help="name of the metadata file", dest='metadata')
    parser.add_argument('-i', '--ioi', action='store', required=True,
                        help="item of interest, column in the metadata file", dest='ioi')
    parser.add_argument('-c', '--current', action='store', required=True,
                        help="name of the current val in the ioi column being evaluated", dest='current')
    args = parser.parse_args()
    main(args)