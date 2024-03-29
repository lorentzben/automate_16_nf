#!/usr/bin/env python3

import PySimpleGUI as sg
import subprocess
import os

#TODO add a box for the ordered ioi and then a call to rename it to where it should be.

sg.theme('System Default')

layout = [[sg.Text('Welcome to the Automate 16s Submitter')],
          [sg.Text('Work Directory', size=(20,1)), sg.Input(), sg.FolderBrowse()],
          [sg.Text('Sequence Directory', size=(20,1)), sg.Input(), sg.FolderBrowse()],
          [sg.Text('Metadata', size=(20,1)), sg.Input(), sg.FileBrowse()],
          [sg.Text('Manifest', size=(20,1)), sg.Input(), sg.FileBrowse()],
          [sg.Text('Item of Interest', size=(20,1)), sg.Input()],
          [sg.Text('Output Directory', size=(20,1)), sg.Input(), sg.FolderBrowse()],
          [sg.Checkbox('Resume',default=False)],
          [sg.Submit(), sg.Cancel()]]

window = sg.Window('Automate 16s nf', layout)

event, values = window.read()
window.close()

print("you chose: " +str( values))



os.chdir(str(values[0]))
if values[6] :
    result = subprocess.run(['nextflow run ~/Documents/automate_16_nf/main.nf --input '+values[1]+" --metadata " +values[2]+" --manifest "+values[3]+" --itemOfInterest "+values[4]+" --outdir "+values[5]+ "--resume"], shell=True)
else:
    result = subprocess.run(['nextflow run ~/Documents/automate_16_nf/main.nf --input '+values[1]+" --metadata " +values[2]+" --manifest "+values[3]+" --itemOfInterest "+values[4]+" --outdir "+values[5]], shell=True)

#TODO run this command in a docker/singulatiry lorentzb/r_latest 

if result.returncode != 0:
    print("error")
else:
    subprocess.run(['bash report.sh'], shell=True)
