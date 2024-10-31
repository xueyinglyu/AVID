# AVID – Accurate Viral Integration Detector

### AVID is a sensitive and accurate tool for viral integration detection by next-generation sequencing data.

# 1. Installation


#### conda config --add channels lyuxueying
#### conda create -n avid
#### conda activate avid
#### conda install lyuxueying::avid



# 2. Test AVID
### AVID package includes a dataset, necessary references, and indexs for testing. Users can directly run the command below to test the installation without any file preparation.
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

### Before running AVID, users must build the index for the host and virus genomes, respectively.
## 3.1 For the host genome
#### BWA: bwa index -a bwtsw -p hg38.fa hg38.fa

                                          
## 3.2 For the virus genome
### Virus genome needs both BWA and BLAST index. 
#### BWA: bwa index -a bwtsw -p HBV.fa HBV.fa
#### BLAST: makeblastdb -in HBV.fa -dbtype nucl -parse_seqids -out HBV

# 4. Annotation files
## 4.1 Human genome
### The annotation file for the human genome can be found in the annotation directory of the AVID package on GitHub or downloaded from the ANNOVAR database (e.g., hg38_refGene.txt).
## 4.2 Virus genome
### Annotation file for virus genome is bed format.
| HBV | 0 | 155 | gene | + | PreS2 |
| -------- | -------- | -------- | -------- | -------- | -------- |
| HBV | 155 | 835 | gene | + | S |
| HBV | 1374 | 1838 | gene | + | X |


# 5. Run AVID
### After preparing the index and annotation files, users can run AVID.

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
## test.AVID.final.txt (test data result)

| chr_human | breakpoint_human | chr_virus | breakpoint_virus | breaktype | virus_direction | chr_human_chimeric | start_human_chimeric | end_human_chimeric | chr_virus_chimeric | start_virus_chimeric | end_virus_chimeric | support_reads_softclip | support_reads_chimeric | total_reads | repeat_mark | region_human | gene_human | gene_virus |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |-------- | -------- | -------- |-------- |-------- | -------- | -------- |
|chr5	|2347	|DQ089769.1	|1746	|upstream	|negative	|0	|0	|0	|0	|0	|0	|9	|0	|9	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89821)	|X|
|chr5	|2465	|DQ089769.1	|1815	|upstream	|negative	|0	|0	|0	|0	|0	|0	|7	|0	|7	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89703)	|X|
|chr5	|2669	|DQ089769.1	|1818	|upstream	|negative	|0	|0	|0	|0	|0	|0	|9	|0	|9	|No	|intergenic	|NONE(dist=NONE),PLEKHG4B(dist=89499)	|X|![image](https://github.com/xueyinglyu/AVID/assets/78334297/f939bcc2-78ae-40bc-847a-199208398c65)

#### chr_human: Chromosome of breakpoint in human genome
#### breakpoint_human: Breakpoint location in human genome
#### chr_virus: Chromosome of breakpoint in virus genome
#### breakpoint_virus: Breakpoint location in virus genome
#### breaktype: Breakpoint type, upstream or downstrem of the breakpoint. There are two breakpoints of a HBV insertion event. Upstream and downstream are used to distinguish them.
#### virus_direction: Strand of virus(positive or negative) inserted into human genome(refer to human positive strand) identified by soft-clip reads
#### chr_human_chimeric: Chromosome of integration region in human genome identified by clustering chimeric reads 
#### start_human_chimeric: Start site of the cluster of chimeric reads in human genome
#### end_human_chimeric: End site of the cluster of chimeric reads in human genome
#### start_virus_chimeric: Start site of the cluster of chimeric reads in virus genome
#### end_virus_chimeric: End site of the cluster of chimeric reads in virus genome
#### support_reads_softclip: The number of soft-clip reads supporting the breakpoint
#### support_reads_chimeric: The number of clustered chimeric reads supporting the breakpoint
#### total_reads: Total number of soft-clip reads and chimeric reads supporting the breakpoint
#### repeat_region: Whether the breakpoint is supported by reads with multi-site mapping (Yes or No). 
#### region_human: Breakpoint in the region of human genome according to annovar annotation result
#### gene_human: Breakpoint to the nearest gene in human genome according to annovar annotation result
#### gene_virus: Breakpoint to the nearest gene in virus genome according to virus annotation file providing by users


# 8. Visualization
### 8.1 Breakpoint distribution in chromosomes in human genome
![Viral integration breakpoint distribution](https://github.com/xueyinglyu/AVID/blob/main/circos.png)
### 8.2 Breakpoint distribution in virus genome
![Breakpoints in viral genome](https://github.com/xueyinglyu/AVID/blob/main/virus_hist.png)

# 9. Tips for AVID users

### 9.1 How to set the parameters?
#### Parameters are typically configured based on users' needs. If users wish to implement strict parameters to minimize false positive events, they can adjust -l, -q, and -t. Increasing the length of the soft-clip part (-l), mapping quality (-q), and the minimum supporting reads (-t) will help achieve this.

### 9.2 How to choose the solid viral integration events from the result?
#### The number of supporting soft-clipped reads (support_reads_softclip) serves as valid evidence of viral integration, and users can increase this number if necessary. Integration events within repeat regions can be filtered out by selecting 'No' in the repeat_region column. There is a high likelihood of false positive events when repetitive breakpoints in the viral genome occur across multiple instances.

### 9.3 Why some breakpoints in human genome are quite close to each other?
#### AVID can identify two breakpoints of an integration event, defined by the upstream and downstream events.













































 
