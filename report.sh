#!/usr/bin/env
conda activate work/r_env*
cp report.Rmd $(cat out.txt)
cp make_report.sh $(cat out.txt)
cd $(cat out.txt)
echo $PWD
bash make_report.sh