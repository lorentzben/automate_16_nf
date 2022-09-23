#!/usr/bin/env bash

#eval "$(conda shell.bash hook)"
#conda activate python2
    
LEFSE=$(which run_lefse.py)
LEFSE_DIR=${LEFSE::-8}
cp plot_res.py $LEFSE_DIR
cp plot_cladogram.py $LEFSE_DIR

cd combos
lefse_files=`ls *lefse_formatted.txt`
for eachfile in $lefse_files
do
    format_input.py $eachfile lefse_formatted.in -c 2 -u 1
    run_lefse.py lefse_formatted.in lefse_result.res 
    plot_res.py lefse_result.res "${eachfile::-4}"_res.png --dpi 100 
    plot_res.py lefse_result.res "${eachfile::-4}"_pdf_r.png --dpi 72
    plot_cladogram.py lefse_result.res "${eachfile::-4}"_cladogram.png --format png --dpi 200
    plot_cladogram.py lefse_result.res "${eachfile::-4}"_pdf_qual_clad.png --format png --dpi 72
    cp lefse_result.res "${eachfile::-4}"_result.res
    cp "${eachfile::-4}"_res.png ../result
    cp "${eachfile::-4}"_pdf_r.png ../result
    cp "${eachfile::-4}"_cladogram.png ../result
    cp "${eachfile::-4}"_pdf_qual_clad.png ../result
done
