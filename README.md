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

#### python AVID.py -1 test_1.fastq -2 test_2.fastq -d testdata -s test -r HBV.fa -l 10 -q 10 -t 1 -@ 1 -v hg38 -p HBV -H hg38.fa -a HBV.bed -R 100 -I 250 -T DNA
## 5. Parameters
#### -1,--fastq1:
#### <br />fastq/fastq.gz for read1
#### -2,--fastq2:
#### <br />fastq/fastq.gz for read2
#### -d,--directory:
#### <br />output directory
#### -s,--sample:
#### sample ID
#### -r,--reference:
#### virus genome index by BWA
#### -p,--pathogen:
#### virus genome index by BLAST
#### -H,--human:
#### host genome index by BWA
#### -l,--len: 
#### minimum length of soft-clip reads for integration detection ()
#### -q,--quality:
#### minimum quality for BWA mapping
#### -t,--threshold:
#### minimum supported reads including soft-clip reads and chimeric reads for integration identification
#### -@,--threads:
#### threads of BWA
#### -v,--version:
#### annovar annotation version(hg19,hg38)
#### -a,--annotation:
#### virus annotation file
#### -R,--Readlen: 
#### sequencing read length (e.g. 150)
#### -I,--insert_size:
#### DNA library fragment size



 
