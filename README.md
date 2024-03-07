# AVID – Accurate Viral Integration Detector

#### AVID is a sensitive and accurate tool for viral integration detection by next generation sequencing data.

## 1. Installation

#### git clone https://github.com/xueyinglyu/AVID.git
#### export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:[your dir]/AVID/"
#### [your dir]: your directory where download AVID





## 2. Test AVID

#### [your dir]: your directory where download AVID
#### [outdir]: output directory
####
#### python [your dir]/AVID/AVID.py 
#### -1 [your dir]/AVID/testdata/test_1.fastq 
#### -2 [your dir]/AVID/testdata/test_2.fastq 
#### -d [outdir] -s test -r [your dir]/AVID/testdata/DQ089769.fasta 
#### -p [your dir]/AVID/testdata/DQ089769 
#### -H [your dir]/AVID/testdata/chr5.fa 
#### -a [your dir]/AVID/testdata/DQ089769.bed 
#### -l 10 -q 10 -t 1 -@ 1 -v hg38 -R 100 -`I` 250 -T DNA



## 3. Build index
### 3.1 For the host genome
#### e.g. bwa index -a bwtsw -p hg38.fa hg38.fa

                                          
### 3.2 For the virus genome
#### e.g. bwa index -a bwtsw -p HBV.fa HBV.fa
#### e.g. makeblastdb -in HBV.fa -dbtype nucl -parse_seqids -out test_data



																				
## 4. Run AVID
### Please refer to 2.Test AVID step

#### python [your dir]/AVID/AVID.py 
#### -1 [your dir]/AVID/testdata/test_1.fastq 
#### -2 [your dir]/AVID/testdata/test_2.fastq 
#### -d [outdir] -s test -r [your dir]/AVID/testdata/DQ089769.fasta 
#### -p [your dir]/AVID/testdata/DQ089769 
#### -H [your dir]/AVID/testdata/chr5.fa 
#### -a [your dir]/AVID/testdata/DQ089769.bed 
#### -l 10 -q 10 -t 1 -@ 1 -v hg38 -R 100 -`I` 250 -T DNA

## 5. Parameters
#### -1,--fastq1:
#### &emsp;fastq/fastq.gz for read1
#### -2,--fastq2:
#### &emsp;fastq/fastq.gz for read2
#### -d,--directory:
#### &emsp;output directory
#### -s,--sample:
#### &emsp;sample ID
#### -r,--reference:
#### &emsp;virus genome index by BWA
#### -p,--pathogen:
#### &emsp;virus genome index by BLAST
#### -H,--human:
#### &emsp;host genome index by BWA
#### -l,--len: 
#### &emsp;minimum length of soft-clip reads for integration detection ()
#### -q,--quality:
#### &emsp;minimum quality for BWA mapping
#### -t,--threshold:
#### &emsp;minimum supported reads including soft-clip reads and chimeric reads for integration identification
#### -@,--threads:
#### &emsp;threads of BWA
#### -v,--version:
#### &emsp;annovar annotation version(hg19,hg38)
#### -a,--annotation:
#### &emsp;virus annotation file
#### -R,--Readlen: 
#### &emsp;sequencing read length (e.g. 150)
#### -`I`,--insert_size:
#### &emsp;DNA library fragment size
#### -T,--data_type:
#### &emsp;aligner difference, DNA (BWA) or RNA (STAR)

## 6. Output format
## sample.virusclip2.final.txt (test data result)

| chr_human | breakpoint_human | chr_virus | breakpoint_virus | breaktype | virus_direction | chr_human_chimeric | start_human_chimeric | end_human_chimeric | chr_virus_chimeric | start_virus_chimeric | end_virus_chimeric | support_reads_softclip | support_reads_chimeric | total_reads | repeat_mark | region_human | gene_human | gene_virus |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |-------- | -------- | -------- |-------- |-------- | -------- | -------- |
|chr5	|2347	|DQ089769.1	|1746	|upstream	|negative	|0	|0	|0	|0	|0	|0	|9	|0	|9	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89821)	|X|
|chr5	|2465	|DQ089769.1	|1815	|upstream	|negative	|0	|0	|0	|0	|0	|0	|7	|0	|7	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89703)	|X|
|chr5	|2669	|DQ089769.1	|1818	|upstream	|negative	|0	|0	|0	|0	|0	|0	|9	|0	|9	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89499)	|X|![image](https://github.com/xueyinglyu/AVID/assets/78334297/f939bcc2-78ae-40bc-847a-199208398c65)

#### chr_human: 
#### &emsp;chromosome of breakpoint in human genome
#### breakpoint_human:
#### &emsp;the breakpoint location in human genome
#### chr_virus:
#### &emsp;chromosome of breakpoint in virus genome
#### breakpoint_virus:
#### &emsp;the breakpoint location in virus genome
#### breaktype:
#### &emsp;breakpoint type, upstream or downstrem of the breakpoint
#### virus_direction:
#### &emsp;the strand of virus(positive or negative) inserted into human genome(refer to human positive strand) identified by soft-clip reads
#### chr_human_chimeric:
#### &emsp;the chromosome of integration region in human genome identified by clustering chimeric reads 
#### start_human_chimeric:
#### &emsp;the start site of the cluster of chimeric reads in human genome
#### end_human_chimeric:
#### &emsp;the end site of the cluster of chimeric reads in human genome
#### start_virus_chimeric:
#### &emsp;the start site of the cluster of chimeric reads in virus genome
#### end_virus_chimeric:
#### &emsp;the end site of the cluster of chimeric reads in virus genome
#### support_reads_softclip:
#### &emsp;the number of soft-clip reads supporting the breakpoint
#### support_reads_chimeric:
#### &emsp;the number of chimeric reads supporting the breakpoint
#### total_reads:
#### &emsp;the total number of soft-clip reads and chimeric reads supporting the breakpoint
#### repeat_region:
#### &emsp;whether the breakpoint is supported by reads with multi-site mapping (Yes or No). 
#### region_human:
#### &emsp;breakpoint in the region of human genome according to annovar annotation result
#### gene_human:
#### &emsp;breakpoint to the nearest gene in human genome according to annovar annotation result
#### gene_virus:
#### &emsp;breakpoint to the nearest gene in virus genome according to virus annotation file providing by users


## 7. Visualization
![Alt Text](https://github.com/xueyinglyu/AVID/tree/main/testdata/Screenshot 2024-03-07 at 18.56.30.png)




 
