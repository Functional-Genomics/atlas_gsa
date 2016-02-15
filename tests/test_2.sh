#!/bin/bash

# Human experiments
source test_shared.sh
species=human
logfile=$species.log
rm -f $logfile

# Collect all analytics from one folder
ls -1 test2/*-analytics.tsv.gz > tsv_files_2
sed "s/-analytics.tsv.gz/-configuration.xml/" tsv_files_2 > xml_files_2
# collect and summarise the data
echo "Indexing..."
if [ ! -e v2_$species.po ]; then
    run_wrapper "$species\tIndex\t." $logfile ../scripts/gsa_prepare_data.R -c 4 -i tsv_files_2  -x xml_files_2 -o $species.po
    # ~15m to index
    echo "Indexing...done."
else
    echo "Found v2_human.po file...skipping indexing"
fi


############################################################
# run
gset0="ENSMUSG00000048852 ENSMUSG00000053931  ENSMUSG00000062991 ENSMUSG00000097928"
gset5="ENSG00000282031 ENSG00000282458 ENSG00000282883 ENSG00000282907 ENSG00000282988"
gset10="ENSG00000280401 ENSG00000280407 ENSG00000280417 ENSG00000280543 ENSG00000280693 ENSG00000280798 ENSG00000280809 ENSG00000280832 ENSG00000281005 ENSG00000281026"
gset100="ENSG00000260837 ENSG00000260852 ENSG00000260862 ENSG00000260907 ENSG00000260916 ENSG00000260941 ENSG00000260966 ENSG00000260966 ENSG00000260971 ENSG00000260992 ENSG00000261040 ENSG00000261051 ENSG00000261061 ENSG00000261087 ENSG00000261097 ENSG00000261098 ENSG00000261116 ENSG00000261136 ENSG00000261150 ENSG00000261150 ENSG00000261150 ENSG00000261167 ENSG00000261183 ENSG00000261189 ENSG00000261213 ENSG00000261236 ENSG00000261242 ENSG00000261254 ENSG00000261268 ENSG00000261295 ENSG00000261314 ENSG00000261324 ENSG00000261326 ENSG00000261329 ENSG00000261335 ENSG00000261338 ENSG00000261342 ENSG00000261359 ENSG00000261371 ENSG00000261371 ENSG00000261373 ENSG00000261386 ENSG00000261441 ENSG00000261474 ENSG00000261485 ENSG00000261488 ENSG00000261488 ENSG00000261490 ENSG00000261526 ENSG00000261553 ENSG00000261609 ENSG00000261609 ENSG00000261613 ENSG00000261652 ENSG00000261737 ENSG00000261742 ENSG00000261757 ENSG00000261780 ENSG00000261786 ENSG00000261795 ENSG00000261799 ENSG00000261801 ENSG00000261824 ENSG00000261824 ENSG00000261845 ENSG00000261845 ENSG00000261970 ENSG00000261971 ENSG00000262089 ENSG00000262152 ENSG00000262228 ENSG00000262251 ENSG00000262406 ENSG00000262580 ENSG00000262655 ENSG00000262919 ENSG00000263002 ENSG00000263002 ENSG00000263013 ENSG00000263063 ENSG00000263146 ENSG00000263276 ENSG00000263276 ENSG00000263335 ENSG00000263391 ENSG00000263412 ENSG00000263465 ENSG00000263465 ENSG00000263528 ENSG00000263528 ENSG00000263624 ENSG00000263731 ENSG00000263731 ENSG00000263823 ENSG00000263874 ENSG00000263900 ENSG00000263900 ENSG00000263900 ENSG00000263900 ENSG00000263900 ENSG00000264070 ENSG00000264112 ENSG00000264176 ENSG00000264207 ENSG00000264207"



for gs in gset0 gset5 gset10 gset100; do
    echo -n $gs
    run_wrapper "$species\t$gs\tv2" $logfile ../scripts/gsa_run_v2.R --db v2_human.po --out res_v2_$gs.tsv --gs "'${!gs}'"
    run_wrapper "$species\t$gs\tv1" $logfile ../scripts/gsa_run.R --db human.po --out res_$gs.tsv --gs "'${!gs}'"
done


cat $logfile
