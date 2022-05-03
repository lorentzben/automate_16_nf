#!/usr/bin/env bash

cd $(cat out.txt)

#report section 1 

if [[ ! -e 01_report ]]; then
    mkdir 01_report
fi

cp ../01_report.Rmd 01_report
cp ../01_make_report.sh 01_report
cd 01_report
singularity run docker://lorentzb/r_01 bash 01_make_report.sh 
cd ..


#report section 2

if [[ ! -e 02_report ]]; then
    mkdir 02_report
fi
cp ../02_report.Rmd 02_report
cp ../02_make_report.sh 02_report
cd 02_report
singularity run docker://lorentzb/r_02 bash 02_make_report.sh
cd ..

#report section 3
if [[ ! -e 03_report ]]; then
    mkdir 03_report
fi
cp ../03_report.Rmd 03_report
cp ../03_make_report.sh 03_report
cd 03_report
singularity run docker://lorentzb/r_03 bash 03_make_report.sh
cd ..

#report section 4

if [[ ! -e 04_report ]]; then
    mkdir 04_report
fi
cp ../04_report.Rmd 04_report
cp ../04_make_report.sh 04_report
cd 04_report
singularity run docker://lorentzb/r_04 bash 04_make_report.sh
cd ..

#report section 5

if [[ ! -e 05_report ]]; then
    mkdir 05_report
fi
cp ../05_report.Rmd 05_report
cp ../05_make_report.sh 05_report
cd 05_report
singularity run docker://lorentzb/r_05 bash 05_make_report.sh
cd ..

#report section 6

if [[ ! -e 06_report ]]; then
    mkdir 06_report
fi
cp ../06_report.Rmd 06_report
cp ../06_make_report.sh 06_report
cd 06_report
singularity run docker://lorentzb/r_06 bash 06_make_report.sh
cd ..

#report section 7 

if [[ ! -e 07_report ]]; then
    mkdir 07_report
fi
cp ../07_report.Rmd 07_report
cp ../07_make_report.sh 07_report
cd 07_report
singularity run docker://lorentzb/r_07 bash 07_make_report.sh
cd ..

#report section 8

if [[ ! -e 08_report ]]; then
    mkdir 08_report
fi
cp ../08_report.Rmd 08_report
cp ../08_make_report.sh 08_report
cd 08_report
singularity run docker://lorentzb/r_08 bash 08_make_report.sh
cd ..

#report section 9

if [[ ! -e 09_report ]]; then
    mkdir 09_report
fi
cp ../09_report.Rmd 09_report
cp ../09_make_report.sh 09_report
cd 09_report
singularity run docker://lorentzb/r_09 bash 09_make_report.sh
cd ..

#report section 10

if [[ ! -e 10_report ]]; then
    mkdir 10_report
fi
cp ../10_report.Rmd 10_report
cp ../10_make_report.sh 10_report
cd 10_report
singularity run docker://lorentzb/r_10 bash 10_make_report.sh
cd ..

#report section 11

if [[ ! -e 11_report ]]; then
    mkdir 11_report
fi
cp ../11_report.Rmd 11_report
cp ../11_make_report.sh 11_report
cd 11_report
singularity run docker://lorentzb/r_11 bash 11_make_report.sh
cd ..

#report section 12

if [[ ! -e 12_report ]]; then
    mkdir 12_report
fi
cp ../12_report.Rmd 12_report
cp ../12_make_report.sh 12_report
cd 12_report
singularity run docker://lorentzb/r_12 bash 12_make_report.sh
cd ..

#report section 13

if [[ ! -e 13_report ]]; then
    mkdir 13_report
fi
cp ../13_report.Rmd 13_report
cp ../13_make_report.sh 13_report
cd 13_report
singularity run docker://lorentzb/r_13 bash 13_make_report.sh
cd ..

#report section 14

if [[ ! -e 14_report ]]; then
    mkdir 14_report
fi
cp ../14_report.Rmd 14_report
cp ../14_make_report.sh 14_report
cd 14_report
singularity run docker://lorentzb/r_14 bash 14_make_report.sh
cd ..

#consolidate html files
if [[ ! -e final_html ]]; then
    mkdir final_html
fi
cp */*.html final_html

#consolidate pdf files

if [[ ! -e final_pdf ]]; then
    mkdir final_pdf
fi
cp */*.pdf final_pdf

