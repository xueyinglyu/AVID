#!/bin/bash
current_dir=$(pwd)
python_executable=$current_dir/3rdParty/bin/python

$python_executable AVID.py \
-1 testdata/test_1.fastq \
-2 testdata/test_2.fastq \
-d test \
-r testdata/DQ089769.fasta \
-p testdata/DQ089769 \
-H testdata/chr5.fa \
-a testdata/DQ089769.bed \
-l 10 -q 10 -@ 1 -R 100 -I 10000 -s test 

