#!/bin/bash

dest=$PWD
echo $dest
#mkdir -p $dest/3rdParty
conda create -p $dest/3rdParty -c bioconda -c conda-forge -y python=3.13 bwa star blast samtools
$dest/3rdParty/bin/pip install numpy==2.2.6 pandas pysam
conda create -p $dest/R -c bioconda -c conda-forge -y r r-circlize r-ggplot2 r-plyr
 



