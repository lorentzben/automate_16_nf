#! /usr/bin/env bash

dt=$(date '+%d-%m-%Y_%H.%M.%S');

Rscript -e "rmarkdown::render('report.Rmd', output_file='report_$dt.html', output_format='html_document', clean=TRUE)"

Rscript -e "rmarkdown::render('report.Rmd', output_file='report_$dt.pdf', output_format='pdf_document', clean=TRUE)"

