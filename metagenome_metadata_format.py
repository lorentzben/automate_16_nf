#!/usr/bin/env python3

import os
from os import path
import argparse
import pandas as pd
import numpy as np 
import glob

def get_file_count(seq_dir):
    #print("returns a count of the number of sequence files in the directory provided by the user.")
    seq_list = os.listdir(seq_dir)
    seqs = []
    for item in seq_list:
        if "fastq.gz" in item:
            seqs.append(item)
        elif "fastq" in item:
            seqs.append(item)
    return(len(seqs))


def get_samp_name(seq_dir):
    #print("returns a list of assumed filenames by truncating the fastq files, should be half of the previous function.")
    seq_list = os.listdir(seq_dir)
    seqs = []
    for item in seq_list:
        if "fastq.gz" in item:
            seqs.append(item)
        elif "fastq" in item:
            seqs.append(item)
    seq_names = []
    for name in seqs:
        edited = (name.split('_S')[0]).replace('-','.')
        edited = edited.replace('_','.')
        edited = edited.replace('/','.')
        seq_names.append(edited)


    return(set(seq_names))
    


def main(arg):
    #print("you submitted: ")
    #print(arg.dir_name)
    #print(arg.bioproj)
    #print(arg.organism)
    #print(arg.host)
    #print(arg.isolation_source)
    #print(arg.samp_date)
    #print(arg.geo_loc)
    #print(arg.lat_long)

    seq_file_count = get_file_count(arg.dir_name)
    samp_names = get_samp_name(arg.dir_name)
    
    if arg.bioproj == None:
        arg.bioproj = ""
    frame = pd.DataFrame(columns=("sample_name","sample_title","bioproject_accession","organism","host","isolation_source","collection_date","geo_loc_name","lat_lon","ref_biomaterial","rel_to_oxygen","samp_collect_device","samp_mat_process","samp_size","source_material_id","description"))
    for samp_name in samp_names:
        values_to_add = {"sample_name":samp_name,"sample_title":"","bioproject_accession":arg.bioproj,"organism":arg.organism,"host":arg.host,"isolation_source":arg.isolation_source,"collection_date":arg.samp_date,"geo_loc_name":arg.geo_loc,"lat_lon":arg.lat_long,"ref_biomaterial":"","rel_to_oxygen":"","samp_collect_device":"","samp_mat_process":"","samp_size":"","source_material_id":"","description":""}
        row_to_add = pd.Series(values_to_add, name = 'samp_name')
        frame = frame.append(row_to_add)

    print("The formatted metadata file is saved to: Formatted.Metagenome.environmental.1.0.tsv")
    frame.to_csv("Formatted.Metagenome.environmental.1.0.tsv", index=False, sep='\t')
        




if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Assists in the Creation of the metadata file for SRA submission")
    parser.add_argument('-n', '--name', action='store', required=True,
                        help="name of the directory that sequences are stored in", dest='dir_name')
    parser.add_argument('-b', '--bioproject', action='store', required=False,
                        help="bioproject accession number that these seqs will be associated with", dest='bioproj')
    parser.add_argument('-o', '--organism', action='store', required=False,
                        help="type of organism sampled, fill if consistant across all samples", dest="organism")
    parser.add_argument('-s', '--host', action='store', required=False,
                        help="host for the organism isolation, either this or isolation_source req", dest="host")
    parser.add_argument('-i', '--isolation_source', action='store', required=False,
                        help="isolation source of samples, either this or host is req", dest="isolation_source")
    parser.add_argument('-d', '--date', action='store', required=False,
                        help="date that sampling occured, fill if consistant across all samples", dest="samp_date")
    parser.add_argument('-g', '--geolocation', action='store', required=False,
                        help="location of where the samples originated, fill if consistant across all samples", dest="geo_loc")
    parser.add_argument('-l', '--latlong', action='store', required=False,
                        help="latitude and longitude where the samples originated, fill if consistant across all samples", dest="lat_long")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    main(args)