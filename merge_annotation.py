import os,sys
import argparse
from optparse import OptionParser
import pysam
import numpy as np
import pandas as pd 

usage=""
parser=OptionParser(usage=usage)
parser.add_option("-H","--human",dest="human",help="annovar human result variant_function")
parser.add_option("-v","--virus",dest="virus",help="virus annotation bed file, chr start end region strand gene")
parser.add_option("-i","--infile",dest="infile",help="find_breakpoint.py result all_type_result.out")
parser.add_option("-o","--output",dest="output",help="output file")
(options,args)=parser.parse_args()

human=options.human
virus=options.virus
infile=options.infile
output=options.output

infile_data=pd.read_table(infile)

if len(infile_data.index)>0:
	virus_annotation=pd.read_table(virus,names=("chr","start","end","region","strand","gene"))
	human_annovar=pd.read_table(human,names=("region_human","gene_human","chr","start","end","ref","variant"))

	human_annovar=human_annovar.assign(key=['0']*len(human_annovar.index))
	human_annovar["key"]=human_annovar["chr"]+human_annovar["start"].astype(str)
	human_annovar=human_annovar[["key","region_human","gene_human"]]


	infile_data=infile_data.assign(key=['0']*len(infile_data.index))
	infile_data_region=infile_data.copy()
	infile_data_region=infile_data_region[infile_data_region["breakpoint_human"]==0]
	infile_data_region["key"]=infile_data_region["chr_human_chimeric"]+infile_data_region["start_human_chimeric"].astype(str)
	infile_data_breakpoint=infile_data.copy()
	infile_data_breakpoint=infile_data_breakpoint[infile_data_breakpoint["breakpoint_human"]>0]
	infile_data_breakpoint["key"]=infile_data_breakpoint["chr_human"]+infile_data_breakpoint["breakpoint_human"].astype(str)
	#infile_data["key"]=infile_data["chr_human"]+infile_data["breakpoint_human"].astype(str)
	infile_data=pd.concat([infile_data_breakpoint,infile_data_region])
	merge_data2=infile_data.merge(human_annovar,on=["key"])

	#region_data=infile_data[~(infile_data["key"].isin(human_annovar["key"]))]
	#region_data=region_data.assign(region_human=['0']*len(region_data.index),gene_human=['0']*len(region_data.index))
	#merge_data2=pd.concat([merge_data2,region_data])
	merge_data2=merge_data2.assign(region_virus=['0']*len(merge_data2.index),gene_virus=['0']*len(merge_data2.index))

	for index,row in merge_data2.iterrows():
		for i,r in virus_annotation.iterrows():
			if int(row["breakpoint_virus"]) <= int(r["end"]) and int(row["breakpoint_virus"]) >= int(r["start"]):
				merge_data2.loc[index,"region_virus"]=r["region"]
				merge_data2.loc[index,"gene_virus"]=r["gene"]
				break
	merge_data2=merge_data2.drop(columns=["key"])
	merge_data2.to_csv(output,index=False,sep='\t')
else:
	merge_data2=pd.DataFrame()
	merge_data2.to_csv(output,index=False,sep='\t')



