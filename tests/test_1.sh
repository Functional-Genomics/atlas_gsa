#!/bin/bash

# 2 experiments

ls -1 test_1/*-analytics.tsv > tsv_files_1
sed "s/-analytics.tsv/-configuration.xml/" tsv_files_1 > xml_files_1
../scripts/prepare_data.R -c 4 -i tsv_files_1  -x xml_files_1 -o test1.po
