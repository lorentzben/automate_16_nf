#!/usr/bin/env bash

lefse_files=`ls combos/*lefse_formatted.txt`
for eachfile in $lefse_files
do
    format_input.py $eachfile lefse_formatted.in -c 2 -u 1
    run_lefse.py lefse_formatted.in lefse_result.res 
    plot_res.py lefse_result.res result/"${eachfile::-4}"_res.png --dpi 100 
    plot_res.py lefse_result.res result/"${eachfile::-4}"_pdf_r.png --dpi 72
    plot_cladogram.py lefse_result.res result/"${eachfile::-4}"_cladogram.png --format png --dpi 200
    plot_cladogram.py lefse_result.res result/"${eachfile::-4}"_pdf_qual_clad.png --format png --dpi 72
    cp lefse_result.res result/"${eachfile::-4}"_result.res
done
