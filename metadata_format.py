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

def get_file_names(seq_dir):
    seq_list = os.listdir(seq_dir)
    seqs = []
    for item in seq_list:
        if "fastq.gz" in item:
            seqs.append(item)
        elif "fastq" in item:
            seqs.append(item)
    seq_names = {}

    for name in seqs:
        edited = (name.split('_S')[0]).replace('-','.')
        edited = edited.replace('_','.')
        edited = edited.replace('/','.')
        seq_names.setdefault(edited,[]).append(name)
    return(seq_names)
    


def main(arg):
    #TODO add description section
    print("you submitted: ")
    print(arg.dir_name)
    print(arg.metadata)
    print(arg.title)
    print(arg.library_strategy)
    print(arg.library_source)
    print(arg.library_selection)
    print(arg.library_layout)
    print(arg.platform)
    print(arg.model)
    print(arg.file_format)

    num_files = get_file_count(arg.dir_name)
    samp_names = get_samp_name(arg.dir_name)
    file_names = get_file_names(arg.dir_name)
    print(str(num_files) +" \n"+str(samp_names)+" \n"+str(file_names))

    
    frame = pd.DataFrame(columns=("biosample_accession","library_ID","title","library_strategy","library_source","library_selection","library_layout","platform","instrument_model","design_description","filetype","filename","filename2","filename3","filename4","assembly","fasta_file"))
    for name in samp_names:

        fastq_1 = file_names[name][0]
        try:
            fastq_2 = file_names[name][1]
        except IndexError:
            fastq_2 = ""
        try:
            fastq_3 = file_names[name][2]
        except IndexError:
            fastq_3 = ""
        try:
            fastq_4 = file_names[name][3]
        except IndexError:
            fastq_4 = ""

        values_to_add = {"biosample_accession":"","library_ID":name,"title":arg.title,"library_strategy":arg.library_strategy,"library_source":arg.library_source ,"library_selection":arg.library_selection,"library_layout":arg.library_layout,"platform": arg.platform,"instrument_model": arg.model,"design_description":"","filetype":arg.file_format,"filename":fastq_1,"filename2":fastq_2,"filename3":fastq_3,"filename4":fastq_4,"assembly":"","fasta_file":""}
        row_to_add = pd.Series(values_to_add, name = 'samp_name')
        frame = frame.append(row_to_add)

    print("The formatted SRA metadata file is saved to: Formatted_SRA_metadata_acc.tsv")
    frame.to_csv("Formatted_SRA_metadata_acc.tsv", index=False, sep='\t')
    


if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Assists in the Creation of the metadata file for SRA submission")
    parser.add_argument('-n', '--name', action='store', required=True,
                        help="name of the directory that sequences are stored in", dest='dir_name')
    parser.add_argument('-m','--metadata', action='store', required=False,
                        help="name of the metadata file with sample names", dest='metadata')
    parser.add_argument('-t', '--title', action='store', required=False,
                        help="Short description: {methodology} of {organism}: {sample info}", dest='title')
    parser.add_argument('-l', '--library_strategy', action='store', required=False,
                        help="sample library type, dropdown", dest="library_strategy")
    parser.add_argument('-s', '--library_source', action='store', required=False,
                        help="strategy of project, dropdown", dest="library_source")
    parser.add_argument('-b', '--library_selection', action='store', required=False,
                        help="methodology of library selection", dest="library_selection")
    parser.add_argument('-y', '--library_layout', action='store', required=False,
                        help="date that sampling occured, fill if consistant across all samples", dest="library_layout")
    parser.add_argument('-p', '--platform', action='store', required=False,
                        help="platform that sequencing was performed on, dropdown", dest="platform")
    parser.add_argument('-i', '--instrument_model', action='store', required=False,
                        help="type of instrument that the platform was run on, dropdown", dest="model")
    parser.add_argument('-f', '--filetype', action='store', required=False,
                        help="file format that the sequences are stored, dropdown", dest="file_format")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    main(args)