# AVID – Accurate Viral Integration Detector

#### AVID is a sensitive and accurate tools for viral integration detection by next generation sequencing data.

### 1. Installation
##### conda install -c lyuxueying avid

### 2. Run test data
#### sh run_test.sh
### 3. Build index
#### 3.1 For the host genome
#### e.g. bwa index -a bwtsw -p hg38.fa hg38.fa

#### 3.2 For the virus genome
#### e.g. bwa index -a bwtsw -p HBV.fa HBV.fa
#### e.g. makeblastdb -in HBV.fa -dbtype nucl -parse_seqids -out <dirname>

#### 4. Run AVID

### Parameters

