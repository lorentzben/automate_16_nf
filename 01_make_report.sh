#! /usr/bin/env bash

echo "I am Here:"
pwd
ls

cp -rf /renv_dev/renv .
cp -rf /renv_dev/renv.lock . 

cp ../item_of_interest.csv .
cp ../order_item_of_interest.csv .

dt=$(date '+%d-%m-%Y_%H.%M.%S');

Rscript -e "rmarkdown::render('01_report.Rmd', output_file='report_$dt.html', output_format='html_document', clean=TRUE)"

Rscript -e "rmarkdown::render('01_report.Rmd', output_file='report_$dt.pdf', output_format='pdf_document', clean=TRUE)"

