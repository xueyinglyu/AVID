import os,sys
import argparse
from optparse import OptionParser
import pysam
import re

usage=""
parser=OptionParser(usage=usage)
parser.add_option("-i","--infile",dest="infile",help="Input files, summary out file")
parser.add_option('-o',"--outfile",dest="outfile",help="Output file")

(options,args)=parser.parse_args()
infile=options.infile
outfile=options.outfile

if not os.path.exists(infile):
	print("input file does not exist (from subtype result)......\n")

IN=open(infile,'r')
OUT=open(outfile,'w')
title=IN.readline()

for line in IN:
    line=line.strip()
    arr=line.split("\t")
    mapping_type=arr[9]
    full_seq=arr[11]
    unmap_seq=arr[12]
    if mapping_type=="F":
        OUT.write(">%s\n%s\n"%(arr[0],full_seq))
    if mapping_type=="P":
        test1=full_seq[0:len(unmap_seq)]
        test2=full_seq[len(full_seq)-len(unmap_seq):len(full_seq)]
        if test1==unmap_seq:
            mapping_seq=full_seq[len(unmap_seq):len(full_seq)]
            OUT.write(">%s\n%s\n"%(arr[0],mapping_seq))
        if test2==unmap_seq:
            mapping_seq=full_seq[0:len(full_seq)-len(unmap_seq)]
            OUT.write(">%s\n%s\n"%(arr[0],mapping_seq))
        if (test1!=unmap_seq) & (test2!=unmap_seq):
            next

IN.close()
OUT.close()
























