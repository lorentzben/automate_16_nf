#!/usr/bin/env bash
cp init_and_refresh.R $(cat out.txt)
cp report.Rmd $(cat out.txt)
cp make_report.sh $(cat out.txt)
cp order_item_of_interest.csv $(cat out.txt)
cd $(cat out.txt)
echo $PWD
Rscript init_and_refresh.R
bash make_report.sh