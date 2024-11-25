import os,sys
import argparse
from optparse import OptionParser
import numpy as np 
import pandas as pd
import re 
sys.setrecursionlimit(10000000)

usage=""
parser=OptionParser(usage=usage)
parser.add_option('-i',"--infile",dest="infile",help="input file")
parser.add_option('-o',"--outfile",dest="outfile",help="output file")
parser.add_option('-t',"--threshold",dest="threshold",help="threshold")
(options,args)=parser.parse_args()
infile=options.infile
outfile=options.outfile
threshold=int(options.threshold)

def primary_merge(dat):
	data_repeat_id=dat.groupby(["repeat_region"]).size().reset_index(name="repeat_count")
	data_repeat_id=data_repeat_id[data_repeat_id["repeat_count"]>1]
	data_repeat_id.index=range(len(data_repeat_id.index))
	data_repeat_id=pd.unique(data_repeat_id["repeat_region"])
	for i in data_repeat_id:
		tmp_dat=dat[dat["repeat_region"]==i]
		chr_human=";".join(tmp_dat["chr_human"].astype("str").tolist())
		breakpoint_human=";".join(tmp_dat["breakpoint_human"].astype("str").tolist())
		chr_virus=";".join(tmp_dat["chr_virus"].astype("str").tolist())
		breakpoint_virus=";".join(tmp_dat["breakpoint_virus"].astype("str").tolist())
		breaktype=";".join(tmp_dat["breaktype"].astype("str").tolist())
		virus_direction=";".join(tmp_dat["virus_direction"].astype("str").tolist())
		chr_human_chimeric=";".join(tmp_dat["chr_human_chimeric"].astype("str").tolist())
		start_human_chimeric=";".join(tmp_dat["start_human_chimeric"].astype("str").tolist())
		end_human_chimeric=";".join(tmp_dat["end_human_chimeric"].astype("str").tolist())
		chr_virus_chimeric=";".join(tmp_dat["chr_virus_chimeric"].astype("str").tolist())
		start_virus_chimeric=";".join(tmp_dat["start_virus_chimeric"].astype("str").tolist())
		end_virus_chimeric=";".join(tmp_dat["end_virus_chimeric"].astype("str").tolist())
		repeat_region="".join(tmp_dat["repeat_region"].astype("str").tolist())
		region_human=";".join(tmp_dat["region_human"].astype("str").tolist())
		gene_human=";".join(tmp_dat["gene_human"].astype("str").tolist())
		region_virus=";".join(tmp_dat["region_virus"].astype("str").tolist())
		repeat_mark=";".join(tmp_dat["repeat_mark"].astype("str").tolist())
		gene_virus=";".join(tmp_dat["gene_virus"].astype("str").tolist())
		tmp_repeat=[]
		repeat_list=repeat_region.split(";")
		repeat_count=0
		no_count=0
		for k in repeat_list:
			if ("repeat" in k) and (k!="") and (k!="No") and (k not in tmp_repeat):
				repeat_count=repeat_count+1
				tmp_repeat.append(k)
			if k=="No":
				no_count=no_count+1
		support_reads_softclip=repeat_count+no_count
		support_reads_chimeric=0
		total_reads=support_reads_softclip
		temp_dataframe=pd.DataFrame(data=dict(zip(dat.columns,[chr_human,breakpoint_human,chr_virus,breakpoint_virus,breaktype,virus_direction,chr_human_chimeric,start_human_chimeric,end_human_chimeric,chr_virus_chimeric,start_virus_chimeric,end_virus_chimeric,support_reads_softclip,support_reads_chimeric,total_reads,repeat_region,repeat_mark,region_human,gene_human,region_virus,gene_virus])),index=[0])
		dat=dat[dat["repeat_region"]!=i]
		dat=pd.concat([temp_dataframe,dat],ignore_index=True)
	return dat




def merge_two_line(dat1,dat2,data_repeat):
	
	arr1=dat1["repeat_region"][0].split(";")
	arr1.pop()
	arr2=dat2["repeat_region"][0].split(";")
	arr2.pop()
	No_count=0
	repeat_count=0
	tmp_list=[]
	#print(arr1)
	#print(arr2)
	for k in arr1:
		if k not in tmp_list and (k !="") and ("repeat" in k):
			tmp_list.append(k)
		if k=="No":
			No_count=No_count+1
			tmp_list.append(k)
	for k in arr2:
		if k not in tmp_list and (k !="") and ("repeat" in k):
			tmp_list.append(k)
		if k=="No":
			No_count=No_count+1
			tmp_list.append(k)
	chr_human=";".join([dat1["chr_human"][0],dat2["chr_human"][0]])
	breakpoint_human=";".join([dat1["breakpoint_human"][0],dat2["breakpoint_human"][0]])
	chr_virus=";".join([dat1["chr_virus"][0],dat2["chr_virus"][0]])
	breakpoint_virus=";".join([dat1["breakpoint_virus"][0],dat2["breakpoint_virus"][0]])
	breaktype=";".join([dat1["breaktype"][0],dat2["breaktype"][0]])
	virus_direction=";".join([dat1["virus_direction"][0],dat2["virus_direction"][0]])
	chr_human_chimeric=";".join([dat1["chr_human_chimeric"][0],dat2["chr_human_chimeric"][0]])
	start_human_chimeric=";".join([dat1["start_human_chimeric"][0],dat2["start_human_chimeric"][0]])
	end_human_chimeric=";".join([dat1["end_human_chimeric"][0],dat2["end_human_chimeric"][0]])
	chr_virus_chimeric=";".join([dat1["chr_virus_chimeric"][0],dat2["chr_virus_chimeric"][0]])
	start_virus_chimeric=";".join([dat1["start_virus_chimeric"][0],dat2["start_virus_chimeric"][0]])
	end_virus_chimeric=";".join([dat1["end_virus_chimeric"][0],dat2["end_virus_chimeric"][0]])
	support_reads_softclip=len(tmp_list)
	support_reads_chimeric=0
	total_reads=support_reads_softclip
	repeat_region=";".join(tmp_list)
	region_human=";".join([dat1["region_human"][0],dat2["region_human"][0]])
	gene_human=";".join([dat1["gene_human"][0],dat2["gene_human"][0]])
	region_virus=";".join([dat1["region_virus"][0],dat2["region_virus"][0]])
	#print(dat1["repeat_mark"][0],dat2["repeat_mark"][0])
	repeat_mark=";".join([dat1["repeat_mark"][0],dat2["repeat_mark"][0]])
	# print("dat1:")
	# print(dat1[["chr_human"]])
	# print(dat1[["breakpoint_human"]])
	# print(dat1[["repeat_region"]])
	# print("dat2:")
	# print(dat2[["chr_human"]])
	# print(dat2[["breakpoint_human"]])
	# print(dat2[["repeat_region"]])
	gene_virus=";".join([str(dat1["gene_virus"][0]),str(dat2["gene_virus"][0])])
	data_repeat=data_repeat[data_repeat["repeat_region"]!=dat1["repeat_region"][0]]
	data_repeat=data_repeat[data_repeat["repeat_region"]!=dat2["repeat_region"][0]]
	temp_dataframe=pd.DataFrame(data=dict(zip(data_repeat.columns,[chr_human,breakpoint_human,chr_virus,breakpoint_virus,breaktype,virus_direction,chr_human_chimeric,start_human_chimeric,end_human_chimeric,chr_virus_chimeric,start_virus_chimeric,end_virus_chimeric,support_reads_softclip,support_reads_chimeric,total_reads,repeat_region,repeat_mark,region_human,gene_human,region_virus,gene_virus])),index=[0])
	# print("merge:")
	# print(temp_dataframe[["chr_human"]])
	# print(temp_dataframe[["breakpoint_human"]])
	# print(temp_dataframe[["repeat_region"]])
	data_repeat=pd.concat([data_repeat,temp_dataframe],ignore_index=True)
	return data_repeat


def findoverlap(data_repeat):
	mark=0
	data_repeat.loc[data_repeat.index,"breakpoint_human"]=data_repeat.loc[data_repeat.index,"breakpoint_human"].apply(str)
	data_repeat.loc[data_repeat.index,"breakpoint_virus"]=data_repeat.loc[data_repeat.index,"breakpoint_virus"].apply(str)
	data_repeat.loc[data_repeat.index,"chr_human_chimeric"]=data_repeat.loc[data_repeat.index,"chr_human_chimeric"].apply(str)
	data_repeat.loc[data_repeat.index,"start_human_chimeric"]=data_repeat.loc[data_repeat.index,"start_human_chimeric"].apply(str)
	data_repeat.loc[data_repeat.index,"end_human_chimeric"]=data_repeat.loc[data_repeat.index,"end_human_chimeric"].apply(str)
	data_repeat.loc[data_repeat.index,"chr_virus_chimeric"]=data_repeat.loc[data_repeat.index,"chr_virus_chimeric"].apply(str)
	data_repeat.loc[data_repeat.index,"start_virus_chimeric"]=data_repeat.loc[data_repeat.index,"start_virus_chimeric"].apply(str)
	data_repeat.loc[data_repeat.index,"end_virus_chimeric"]=data_repeat.loc[data_repeat.index,"end_virus_chimeric"].apply(str)
	data_repeat_id=data_repeat["repeat_region"]
	for index,row in data_repeat.iterrows():
		kk_index=data_repeat.loc[index,"repeat_region"]
		kk_index_list=kk_index.split(";")
		kk_index_list.pop()
		for k in data_repeat_id:
			for kk in kk_index_list:
				if (kk+";" in k) and (data_repeat.loc[index,"repeat_region"]!=k) and (kk!="No"):
					tmp1=data_repeat[data_repeat["repeat_region"]==data_repeat.loc[index,"repeat_region"]].copy()
					tmp1.index=range(len(tmp1.index))
					#print(tmp1[["breakpoint_human","repeat_region"]])
					#print(kk)
					#print(k)
					tmp2=data_repeat[data_repeat["repeat_region"]==k].copy()
					tmp2.index=range(len(tmp2.index))
					#print(tmp2[["breakpoint_human","repeat_region"]])
					# print("tmp1:")
					# print(tmp1[["chr_human"]])
					# print(tmp1[["breakpoint_human"]])
					# print(tmp1[["repeat_region"]])
					# print("tmp2:")
					# print(tmp2[["chr_human"]])
					# print(tmp2[["breakpoint_human"]])
					# print(tmp2[["repeat_region"]])
					data_repeat=merge_two_line(tmp1,tmp2,data_repeat)
					mark=1
					return findoverlap(data_repeat)
	if mark==0:
		return data_repeat



in_data=pd.read_table(infile)

if len(in_data.index)>0:


	data_uniq=in_data.copy()
	data_repeat=in_data.copy()
	data_uniq=data_uniq[(data_uniq["repeat_mark"]=="No")]
	data_uniq=data_uniq.drop(["repeat_region"], axis=1)

	data_repeat=data_repeat[~(data_repeat["repeat_mark"]=="No")]


	if len(data_repeat)>0:
		data_repeat["repeat_region"]=data_repeat["repeat_region"]+";"
		data_repeat=primary_merge(data_repeat)
		data_repeat_out=findoverlap(data_repeat)
		data_repeat_out["repeat_mark"]=['Yes']*len(data_repeat_out.index)
		#data_repeat_out[data_repeat_out["chr_human"].map(lambda x:(";" in x))]["repeat_mark"]=['Possible']*sum(data_repeat_out["chr_human"].map(lambda x:(";" in x)))
		data_repeat_out=pd.concat([data_uniq,data_repeat_out],ignore_index=True)
		data_repeat_out=data_repeat_out[data_repeat_out["total_reads"]>=threshold]
		data_repeat_out=data_repeat_out.drop(["repeat_region","region_virus"], axis=1)

		data_repeat_out.loc[~(data_repeat_out['chr_human'].astype(str).str.contains(";")),"repeat_mark"]="No"
		data_repeat_out.to_csv(outfile,index=False,sep='\t')
	else:
		data_uniq=data_uniq.drop(["region_virus"], axis=1)
		data_uniq.to_csv(outfile,index=False,sep='\t')	

else:
	data_uniq=pd.DataFrame()
	data_uniq.to_csv(outfile,index=False,sep='\t')







