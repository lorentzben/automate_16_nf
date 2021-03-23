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

layout = [[sg.Text('SRA Bioproject Metadata Editor (1)')],
          [sg.Text('Seq Dir', size=(20,1)), sg.FolderBrowse()],
          [sg.Text('Bioproj Accession', size=(20,1)), sg.Input()],
          [sg.Text('Organism', size=(20,1)), sg.Input()],
          [sg.Text('Host*', size=(20,1)), sg.Input()],
          [sg.Text('Isolation Source*', size=(20,1)), sg.Input()],
          [sg.Text('Sample Date', size=(20,1)), sg.Input()],
          [sg.Text('Geolocation', size=(20,1)), sg.Input()],
          [sg.Text('LatLong', size=(20,1)), sg.Input()],
          [sg.Text('Outfile Name', size=(20,1)), sg.Input()],
          [sg.Submit(), sg.Cancel()]]

def get_file_count(seq_dir):
    #print("returns a count of the number of sequence files in the directory provided by the user.")
    try:
        seq_list = os.listdir(seq_dir)
    except FileNotFound:
        print("could not find dir provided, or no dir provided")
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
    


#print(arg.dir_name)
#print(arg.bioproj)
#print(arg.organism)
#print(arg.host)
#print(arg.isolation_source)
#print(arg.samp_date)
#print(arg.geo_loc)
#print(arg.lat_long)
window = sg.Window('Automate 16s nf', layout)
event, values = window.read()
window.close()

print("you chose: " +str(values))
if values['Browse'] != None:
    dir_name = values['Browse']
bioproj = values[0]
organism = values[1]
host = values[2]
isolation_source = values[3]
samp_date = values[4]
geo_loc = values[5]
lat_long = values[6]
outfile_name = values[7]
try:
    seq_file_count = get_file_count(dir_name)
    samp_names = get_samp_name(dir_name)
except NameError:
    print('whoopsie')
    samp_names = []
    
if bioproj == None:
    bioproj = ""
frame = pd.DataFrame(columns=("sample_name","sample_title","bioproject_accession","organism","host","isolation_source","collection_date","geo_loc_name","lat_lon","ref_biomaterial","rel_to_oxygen","samp_collect_device","samp_mat_process","samp_size","source_material_id","description"))
for samp_name in samp_names:
    values_to_add = {"sample_name":samp_name,"sample_title":"","bioproject_accession":bioproj,"organism":organism,"host":host,"isolation_source":isolation_source,"collection_date":samp_date,"geo_loc_name":geo_loc,"lat_lon":lat_long,"ref_biomaterial":"","rel_to_oxygen":"","samp_collect_device":"","samp_mat_process":"","samp_size":"","source_material_id":"","description":""}
    row_to_add = pd.Series(values_to_add, name = 'samp_name')
    frame = frame.append(row_to_add)

if outfile_name == None:
    exit(1)
elif outfile_name == '':
    print("The formatted metadata file is saved to: Formatted.Metagenome.environmental.1.0.tsv")
    frame.to_csv("Formatted.Metagenome.environmental.1.0.tsv", index=False, sep='\t')
else:
    try:
        print("The formatted metadata file is saved to: "+outfile_name+".tsv")
        frame.to_csv(outfile_name+".tsv", index=False, sep='\t')
    except TypeError:
        print("whoopsie")
