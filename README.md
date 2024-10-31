# AVID – Accurate Viral Integration Detector

### AVID is a sensitive and accurate tool for viral integration detection by next generation sequencing data.

# 1. Installation

## 1.1 Conda

### 1.1.1 Add channel:
#### conda config --add channels lyuxueying
### 1.1.2 Install AVID:
#### conda create -n avid
#### conda activate avid
#### conda install lyuxueying::avid
### 1.1.3 Export library:
#### export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:[your conda path]/bin/AVID/"
#### (e.g. export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/lyuxueying/anaconda3/envs/avid/bin/AVID/")


## 1.2 Github
#### 1.2.1 git clone https://github.com/xueyinglyu/AVID.git
#### 1.2.2 export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:[your dir]/AVID/"
#### (e.g. export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/lyuxueying/test/AVID/")
#### [your dir]: Directory where you download AVID.



# 2. Test AVID
## AVID package includes a dataset, necessary reference, and index for testing. Users can run the command below directly to test the installation without file preparation.
#### 
#### [your dir]: your directory where download AVID
#### [outdir]: output directory
####
#### python [your dir]/AVID/AVID.py 
#### -1 [your dir]/AVID/testdata/test_1.fastq 
#### -2 [your dir]/AVID/testdata/test_2.fastq 
#### -d [outdir] -r [your dir]/AVID/testdata/DQ089769.fasta 
#### -p [your dir]/AVID/testdata/DQ089769 
#### -H [your dir]/AVID/testdata/chr5.fa 
#### -a [your dir]/AVID/testdata/DQ089769.bed 
#### -l 10 -q 10 -@ 1 -R 100 -`I` 10000 -s test 
####
#### (Final output file: test.AVID.final.txt)

#### Example:
#### python /home/lyuxueying/anaconda3/envs/avid/bin/AVID/AVID.py 
#### -1 /home/lyuxueying/anaconda3/envs/avid/bin/AVID/testdata/test_1.fastq 
#### -2 /home/lyuxueying/anaconda3/envs/avid/bin/AVID/testdata/test_2.fastq 
#### -d /home/lyuxueying/test/avid/ 
#### -r /home/lyuxueying/anaconda3/envs/avid/bin/AVID/testdata/DQ089769.fasta 
#### -p /home/lyuxueying/anaconda3/envs/avid/bin/AVID/testdata/DQ089769 
#### -H /home/lyuxueying/anaconda3/envs/avid/bin/AVID/testdata/chr5.fa 
#### -a /home/lyuxueying/anaconda3/envs/avid/bin/AVID/testdata/DQ089769.bed 
#### -l 10 -q 10 -@ 1 -R 100 -`I` 10000 -s test 



# 3. Build index

## Before running AVID, users must build the index for the host and virus genomes, respectively.
## 3.1 For the host genome
#### BWA: bwa index -a bwtsw -p hg38.fa hg38.fa

                                          
## 3.2 For the virus genome
### Virus genome needs both BWA and BLAST index. 
#### BWA: bwa index -a bwtsw -p HBV.fa HBV.fa
#### BLAST: makeblastdb -in HBV.fa -dbtype nucl -parse_seqids -out HBV

# 4. Annotation files
## 4.1 Human genome
### Annotation file for the human genome could be found in annotation directory in AVID package or downloaded from annovar database(e.g. hg38_refGene.txt).
## 4.2 Virus genome
### Annotation file for virus genome is bed format.
| HBV | 0 | 155 | gene | + | PreS2 |
| -------- | -------- | -------- | -------- | -------- | -------- |
| HBV | 155 | 835 | gene | + | S |
| HBV | 1374 | 1838 | gene | + | X |


# 5. Run AVID
### Once index files and annotation files are prepared, users can run AVID.

#### python [your dir]/AVID/AVID.py 
#### -1 [your dir]/AVID/testdata/test_1.fastq 
#### -2 [your dir]/AVID/testdata/test_2.fastq 
#### -d [outdir] -s test -r [your dir]/AVID/testdata/DQ089769.fasta 
#### -p [your dir]/AVID/testdata/DQ089769 
#### -H [your dir]/AVID/testdata/chr5.fa 
#### -a [your dir]/AVID/testdata/DQ089769.bed 
#### -l 10 -q 10 -t 1 -@ 1 -v hg38 -R 100 -`I` 250 -T DNA

# 6. Parameters
#### -1,--fastq1:
#### &emsp;fastq/fastq.gz for read1
#### -2,--fastq2:
#### &emsp;fastq/fastq.gz for read2
#### -d,--directory:
#### &emsp;Output directory
#### -s,--sample:
#### &emsp;Sample ID
#### -r,--reference:
#### &emsp;Virus genome index by BWA
#### -p,--pathogen:
#### &emsp;Virus genome index by BLAST
#### -H,--human:
#### &emsp;Host genome index by BWA
#### -l,--len: 
#### &emsp;Minimum length(bp) of soft-clip reads for integration detection (default: 10)
#### -q,--quality:
#### &emsp;Minimum quality for BWA mapping
#### -t,--threshold:
#### &emsp;Minimum supported reads including soft-clip reads and chimeric reads for integration identification (default: 3)
#### -@,--threads:
#### &emsp;Threads of BWA
#### -v,--version:
#### &emsp;Annovar annotation version[hg19 or hg38, default: hg38]
#### -a,--annotation:
#### &emsp;virus annotation file, bed format
#### -R,--Readlen: 
#### &emsp;Sequencing read length (e.g. 150)
#### -`I`,--insert_size:
#### &emsp;DNA library fragment size (e.g. 100000)
#### -T,--data_type:
#### &emsp;Aligner difference, DNA (BWA) or RNA (STAR), default:DNA
#### -F,--chr-incompatibility:
#### &emsp;Activate the filtering of chromosomal incompatibility (False or True, default: True)
#### -S,--filter-insertsize:
#### &emsp;Activate the filtering of DNA fragment insert size (False or True, default: True)

# 7. Output format
# test.AVID.final.txt (test data result)

| chr_human | breakpoint_human | chr_virus | breakpoint_virus | breaktype | virus_direction | chr_human_chimeric | start_human_chimeric | end_human_chimeric | chr_virus_chimeric | start_virus_chimeric | end_virus_chimeric | support_reads_softclip | support_reads_chimeric | total_reads | repeat_mark | region_human | gene_human | gene_virus |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |-------- | -------- | -------- |-------- |-------- | -------- | -------- |
|chr5	|2347	|DQ089769.1	|1746	|upstream	|negative	|0	|0	|0	|0	|0	|0	|9	|0	|9	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89821)	|X|
|chr5	|2465	|DQ089769.1	|1815	|upstream	|negative	|0	|0	|0	|0	|0	|0	|7	|0	|7	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89703)	|X|
|chr5	|2669	|DQ089769.1	|1818	|upstream	|negative	|0	|0	|0	|0	|0	|0	|9	|0	|9	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89499)	|X|![image](https://github.com/xueyinglyu/AVID/assets/78334297/f939bcc2-78ae-40bc-847a-199208398c65)

#### chr_human: 
#### &emsp;Chromosome of breakpoint in human genome
#### breakpoint_human:
#### &emsp;Breakpoint location in human genome
#### chr_virus:
#### &emsp;Chromosome of breakpoint in virus genome
#### breakpoint_virus:
#### &emsp;Breakpoint location in virus genome
#### breaktype:
#### &emsp;breakpoint type, upstream or downstrem of the breakpoint. There are two breakpoints of a HBV insertion event. Upstream and downstream are used to distinguish them.
#### virus_direction:
#### &emsp;Strand of virus(positive or negative) inserted into human genome(refer to human positive strand) identified by soft-clip reads
#### chr_human_chimeric:
#### &emsp;Chromosome of integration region in human genome identified by clustering chimeric reads 
#### start_human_chimeric:
#### &emsp;Start site of the cluster of chimeric reads in human genome
#### end_human_chimeric:
#### &emsp;End site of the cluster of chimeric reads in human genome
#### start_virus_chimeric:
#### &emsp;Start site of the cluster of chimeric reads in virus genome
#### end_virus_chimeric:
#### &emsp;End site of the cluster of chimeric reads in virus genome
#### support_reads_softclip:
#### &emsp;The number of soft-clip reads supporting the breakpoint
#### support_reads_chimeric:
#### &emsp;The number of clustered chimeric reads supporting the breakpoint
#### total_reads:
#### &emsp;Total number of soft-clip reads and chimeric reads supporting the breakpoint
#### repeat_region:
#### &emsp;Whether the breakpoint is supported by reads with multi-site mapping (Yes or No). 
#### region_human:
#### &emsp;Breakpoint in the region of human genome according to annovar annotation result
#### gene_human:
#### &emsp;Breakpoint to the nearest gene in human genome according to annovar annotation result
#### gene_virus:
#### &emsp;Breakpoint to the nearest gene in virus genome according to virus annotation file providing by users


# 8. Visualization
### 8.1 Breakpoint distribution in chromosomes in human genome
![Viral integration breakpoint distribution](https://github.com/xueyinglyu/AVID/blob/main/circos.png)
### 8.2 Breakpoint distribution in virus genome
![Breakpoints in viral genome](https://github.com/xueyinglyu/AVID/blob/main/virus_hist.png)

# 9. Promblem shooting




 
