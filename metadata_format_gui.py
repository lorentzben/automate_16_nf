#!/usr/bin/env python3

import PySimpleGUI as sg
import subprocess
import os
from os import path
import argparse
import pandas as pd
import numpy as np 
import glob

sg.theme('SystemDefault')

layout = [[sg.Text('SRA BioSample Metadata Editor (2)')],
          [sg.Text('Seq Dir', size=(20,1)), sg.FolderBrowse()],
          [sg.Text('Sample Title', size=(20,1)), sg.Input()],
          [sg.Text('Library Strategy')],
          [sg.Combo(['WGS','RNA-Seq','AMPLICON','OTHER'])],
          [sg.Text('Library Source')],
          [sg.Combo(['GENOMIC','TRANSCRIPTOMIC','METAGENOMIC','METATRANSCRIPTOMIC','SYNTHETIC','VIRAL RNA','GENOMIC SINGLE CELL', 'TRANSRIPTOMIC SINGLE CELL','OTHER'])],
          [sg.Text('Library Selection')],
          [sg.Combo(['RANDOM','PCR','RANDOM PCR','RT-PCR','HMPR','MF','CF-S','CF-M','CF-T','MDA','MSLL','cDNA','ChIP',',MNase','DNAase', 'Hybrid Selection','Reduced Representation', 'Restriction Digest','5-methylcytidine antibody','MBD2 protein methyl-CpG binding domain','CAGE','RACE','size fractionation','Padlock probes capture method', 'other','unspecified','cDNA_oligo_dT','cDNA_randomPriming','Inverse rRNA','Oligo-dT','PolyA','repeat fractionation'])],
          [sg.Text('Library Layout')],
          [sg.Combo(['single', 'paired'])],
          [sg.Text('Platform')],
          [sg.Combo(['ILLUMINA','ION_TORRENT','OXFORD_NANOPORE','PACBIO_SMRT'])],
          [sg.Text('Model')],
          [sg.Combo(['HiSeq X Five','HiSeq X Ten','Illumina Genome Analyzer','Illumina Genome Analyzer II','Illumina Genome Analyzer IIx','Illumina HiScanSQ','Illumina HiSeq 1000','Illumina HiSeq 1500','Illumina HiSeq 2000','Illumina HiSeq 2500','Illumina HiSeq 3000','Illumina HiSeq 4000','Illumina iSeq 100','Illumina NovaSeq 6000','Illumina MiniSeq','Illumina MiSeq','NextSeq 500','NextSeq 550', 'Ion Torrent PGM','Ion Torrent Proton','Ion Torrent S5 XL','Ion Torrent S5','GridION','MinION','PromethION','PacBio RS','PacBio RS II','PacBio Sequel','PacBio Sequel II'])],
          [sg.Text('Filetype')],
          [sg.Combo(['bam','srf','sff','fastq','454_native','Helicos_native','SOLiD_native','PacBio_HDF5','CompleteGenomics_native','OxfordNanopore_native'])],
          [sg.Text('Outfile Name', size=(20,1)), sg.Input()],
          [sg.Submit(), sg.Cancel()]]

window = sg.Window('Automate 16s nf', layout)
event, values = window.read()
window.close()

print("you chose: " +str(values))

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
    


#TODO add description section
dir_name = values['Browse']
title = values[0]
library_strategy = values[1]
library_source = values[2]
library_selection = values[3]
library_layout = values[4]
platform = values[5]
model = values[6]
file_format = values[7]
outfile_name = values[8]

num_files = get_file_count(dir_name)
samp_names = get_samp_name(dir_name)
file_names = get_file_names(dir_name)
#print(str(num_files) +" \n"+str(samp_names)+" \n"+str(file_names))

    
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

    values_to_add = {"biosample_accession":"","library_ID":name,"title":title,"library_strategy":library_strategy,"library_source":library_source ,"library_selection":library_selection,"library_layout":library_layout,"platform": platform,"instrument_model": model,"design_description":"","filetype":file_format,"filename":fastq_1,"filename2":fastq_2,"filename3":fastq_3,"filename4":fastq_4,"assembly":"","fasta_file":""}
    row_to_add = pd.Series(values_to_add, name = 'samp_name')
    frame = frame.append(row_to_add)

if outfile_name == "":
    print("The formatted SRA metadata file is saved to: Formatted_SRA_metadata_acc.tsv")
    frame.to_csv("Formatted_SRA_metadata_acc.tsv", index=False, sep='\t')
else:
    print("The formatted SRA metadata file is saved to: "+outfile_name+".tsv")
    frame.to_csv(outfile_name+".tsv",index=False,sep='\t')
