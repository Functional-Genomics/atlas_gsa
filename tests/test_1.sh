#!/bin/bash

# 2 experiments

# collect all analytics from one folder
ls -1 test_1/*-analytics.tsv > tsv_files_1
# collect and summarise the data
../scripts/gsa_prepare_data.R -c 4 -i tsv_files_1 -o test1.po

# run
../scripts/gsa_run.R --db test1.po --out TEST --gs "ENSMUSG00000048852 ENSMUSG00000053931  ENSMUSG00000062991 ENSMUSG00000097928"
