#!/usr/bin/env bash

cp -rf /renv_dev/renv .
cp -rf /renv_dev/renv.lock .
Rscript -e "renv::restore(library='./renv/library/R-4.1/x86_64-pc-linux-gnu/', lockfile='./renv.lock')"
echo 'done' > set.txt