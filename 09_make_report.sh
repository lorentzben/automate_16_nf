#! /usr/bin/env bash

echo "I am Here:"
pwd
ls

cp -rf /renv_dev/renv .
cp -rf /renv_dev/renv.lock . 

Rscript -e "renv::init()"
Rscript -e "renv::install('rmarkdown')"

cp ../item_of_interest.csv .
cp ../order_item_of_interest.csv .

dt=$(date '+%d-%m-%Y_%H.%M.%S');

Rscript -e "rmarkdown::render('09_report.Rmd', output_file='09_report_$dt.html', output_format='html_document', clean=TRUE)"

Rscript -e "rmarkdown::render('09_report.Rmd', output_file='09_report_$dt.pdf', output_format='pdf_document', clean=TRUE)"

