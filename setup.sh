#!/bin/bash

dest=$PWD
echo $dest
#mkdir -p $dest/3rdParty
conda create -p $dest/3rdParty -c bioconda -c conda-forge -y bwa star blast samtools python
$dest/3rdParty/bin/pip install numpy pandas pysam
conda create -p $dest/R -c bioconda -c conda-forge -y r r-circlize r-base r-ggplot2 r-plyr
 



