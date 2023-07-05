# AVID – Accurate Viral Integration Detector

#### AVID is a sensitive and accurate tools for viral integration detection by next generation sequencing data.

## 1. Installation
##### conda install -c lyuxueying avid

## 2. Run test data
#### sh run_test.sh
## 3. Build index
### 3.1 For the host genome
#### e.g. bwa index -a bwtsw -p hg38.fa hg38.fa

                                          
### 3.2 For the virus genome
#### e.g. bwa index -a bwtsw -p HBV.fa HBV.fa
#### e.g. makeblastdb -in HBV.fa -dbtype nucl -parse_seqids -out test_data


                                        
## 4. Run AVID

### python AVID.py -1 test_1.fastq -2 test_2.fastq -d testdata -s test -r HBV.fa -l 10 -q 10 -t 1 -@ 1 -v hg38 -p HBV -H hg38.fa -a HBV.bed -R 100 -I 250 -T DNA
## 5. Parameters

