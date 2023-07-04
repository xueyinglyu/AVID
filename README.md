# AVID – Accurate Viral Integration Detector

#### AVID is a sensitive and accurate tools for viral integration detection by next generation sequencing data.

### Installation
##### conda install -c lyuxueying avid

### Test data
#### python AVID.py -1 ./testdata/test_1.fastq -2 ./testdata/test_2.fastq -d ./testdata/ -s test -r ./testdata/DQ089769.fasta -l 10 -q 10 -t 1 -@ 1 -v hg38 -p ./testdata/DQ089769 -H ./testdata/chr5.fa -a ./testdata/DQ089769.bed -R 100 -I 250 -T DNA
