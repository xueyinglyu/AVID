import os,sys
import argparse
from optparse import OptionParser
import numpy as np 
import pandas as pd 

usage=""
parser=OptionParser(usage=usage)
parser.add_option("-i","--infile",dest="infile",help="Input files, sam file")
parser.add_option('-d',"--outdir",dest="outdir",help="Output dir")
parser.add_option('-s',"--sample",dest="sample",help="output file prefix")
parser.add_option('-l',"--slen",dest="slen",help="minimun softclip length")

(options,args)=parser.parse_args()
infile=options.infile
outdir=options.outdir
sample=options.sample
slen=int(float(options.slen))

def Tofasta(datframe,outfile):
	tmpframe=datframe
	tmpframe.columns=["read_name","seq"]
	OUT=open(outfile,'w')
	for index,row in datframe.iterrows():
		OUT.write(">%s\n%s\n"%(row["read_name"],row["seq"]))
	OUT.close()


data_result=pd.read_table(infile,engine='python')
data_result["softclip_len"]=data_result["softclip_len"].astype(int)
unmapped_data=data_result[data_result["mapping_type"]=="U"]
partial_data=data_result[data_result["mapping_type"]=="P"]
#partial_data=partial_data[partial_data["softclip_len"]>=8]
#partial_data_bwa=partial_data[partial_data["softclip_len"]>6]
#partial_data_blastn=partial_data[partial_data["softclip_len"]<=6]
full_data=data_result[data_result["mapping_type"]=="F"]


unmapp_outfile=outdir+"/"+sample+"_Unmapped.fasta"
full_outfile=outdir+"/"+sample+"_Full.fasta"

Tofasta(unmapped_data[["read_name","softclip_seq"]],unmapp_outfile)

full_data.to_csv("%s/%s_Full.out"%(outdir,sample),sep='\t',index=False)
unmapped_data.to_csv("%s/%s_Unmapped.out"%(outdir,sample),sep='\t',index=False)
partial_data.to_csv("%s/%s_Partial.out"%(outdir,sample),sep='\t',index=False)

if sample=="HBV":
	partial_data=partial_data[partial_data["softclip_len"]>=slen]
	partial_data_bwa=partial_data[partial_data["softclip_len"]>26]
	partial_data_blastn=partial_data[partial_data["softclip_len"]<=26]
	partial_bwa_outfile=outdir+"/"+sample+"_Partial_BWA.fasta"
	partial_blastn_outfile=outdir+"/"+sample+"_Partial_BLASTN.fasta"
	Tofasta(partial_data_bwa[["read_name","softclip_seq"]],partial_bwa_outfile)
	Tofasta(partial_data_blastn[["read_name","softclip_seq"]],partial_blastn_outfile)
else:
	partial_data=partial_data[partial_data["softclip_len"]>8]
	partial_outfile=outdir+"/"+sample+"_Partial.fasta"
	Tofasta(partial_data[["read_name","softclip_seq"]],partial_outfile)


























