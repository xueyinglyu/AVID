import os,sys
import argparse
from optparse import OptionParser
import pysam
import re

usage=""
parser=OptionParser(usage=usage)
parser.add_option('-q',"--quality",dest="quality",help="mapping quality")
parser.add_option("-i","--infile",dest="infile",help="Input files, sam file")
parser.add_option('-o',"--outfile",dest="outfile",help="Output file")

(options,args)=parser.parse_args()
quality=options.quality
infile=options.infile
outfile=options.outfile

if not os.path.exists(infile):
	print("input file does not exist......\n")

def define_pair(flag,read_id):
	mark="None"
	if flag & 64 :
		mark="read1"
	if flag & 128 :
		mark="read2"
	mark=read_id+"_"+mark
	return mark

def define_strand(flag):
	if flag & 16:
		mark="minus"
	else:
		mark="plus"
	return mark

def define_map(cigar,flag):
	if flag & 4:
		mark="U"
		return mark
	else:
		if "S" in cigar:
			mark="P"
			return mark
		else:
			mark="F"
			return mark

def softclip_one(cigar_info):
	soft_clip_start=0
	soft_clip_end=0
	match_base=0
	count=0
	mark="S"
	for k in cigar_info:
		count=count+1
		if k[1]=='M' and count>1:
			match_base=int(k[0])
		if k[1]=="M" and count==1:
			soft_clip_start=soft_clip_start+int(k[0])
			soft_clip_end=soft_clip_end+int(k[0])
			match_base=int(k[0])
			mark="M"
		if k[1]=="S":
			soft_clip_end=soft_clip_end+int(k[0])
		if k[1]=="I":
			soft_clip_start=soft_clip_start+int(k[0])
			soft_clip_end=soft_clip_end+int(k[0])
	return (soft_clip_start,soft_clip_end,match_base,mark)


def softclip_twice(cigar_info):
	soft_clip_start=0
	soft_clip_end=0
	soft_count=0
	match_base=0
	mark="S"
	count=0
	exist=0
	for k in cigar_info:
		count=count+1
		if k[1]=='M':
			soft_clip_start=soft_clip_start+int(k[0])
			soft_clip_end=soft_clip_end+int(k[0])
			match_base=int(k[0])
		if k[1]=="M" and count==1:
			mark="M"
		if k[1]=="S" and soft_count==0:
			soft_clip_end=soft_clip_end+int(k[0])
			softlen=int(k[0])
			if softlen > 6:
				exist=exist+1
			soft_clip_first_end=soft_clip_end
			soft_clip_first_start=soft_clip_start
			soft_clip_start=soft_clip_start+int(k[0])
			soft_count+=1
			continue
		if k[1]=="S" and soft_count>0:
			if int(k[0])>6:
				exist=exist+1
			if int(k[0])>=softlen:
				soft_clip_end=soft_clip_end+int(k[0])
				return(soft_clip_start,soft_clip_end,match_base,mark,exist)
			else:
				return(soft_clip_first_start,soft_clip_first_end,match_base,mark,exist)

		if k[1]=="I":
			soft_clip_start=soft_clip_start+int(k[0])
			soft_clip_end=soft_clip_end+int(k[0])




def Complement_Reverse(seq):
	trantab=str.maketrans('ATGCatgcnN','TACGtacgnN')
	seq_c=seq.translate(trantab)
	seq_c_r=seq_c[::-1]
	return(seq_c_r)


def softclip_region(start,end,strand):
	if start<=5 and strand=="plus":
		soft_human="forward"
	if start<=5 and strand =="minus":
		soft_human="backward"
	if start>5 and strand =="plus":
		soft_human="backward"
	if start>5 and strand =="minus":
		soft_human="forward"
	return soft_human


def deal_line(line,mark,write_h,quality):
	
	line=line.strip()
	arr=line.split("\t")
	read_id=arr[0]
	read_flag=int(arr[1])
	read_chr=arr[2]
	read_start=int(arr[3])
	read_mqual=int(arr[4])
	if read_mqual<=int(quality) and mark!="Yes":
		return "None"
	read_cigar=arr[5]
	mate_chr=arr[6]
	mate_start=int(arr[7])
	insert_size=int(arr[8])
	sequence=arr[9]
	multi_mapping_site=mark
	mapping_type=define_map(read_cigar,read_flag)
	read_strand=define_strand(read_flag)
	if "_read" in read_id:
		read_name=read_id
		read_id=read_id[0:len(read_id)-6]
	else:
		read_name=define_pair(read_flag,read_id)
	cigar_info=re.findall(r'(\d+)(\w)',read_cigar)
	softclip_count=re.findall(r'(\d+)S',read_cigar)
	if len(softclip_count)==1:
		(softclip_start,softclip_end,match_base,mark)=softclip_one(cigar_info)
	if len(softclip_count)>1:
		(softclip_start,softclip_end,match_base,mark,exist)=softclip_twice(cigar_info)
		if exist>1:
			return "None"
	if len(softclip_count)==0:
		softclip_start=0
		softclip_end=len(sequence)
		match_base=softclip_end
		mark="M"
	softclip_seq=sequence[int(softclip_start):int(softclip_end)]
	# if mark=="M" and len(softclip_count) >0:
	# 	softclip_site=read_start+match_base
	# else:
	# 	softclip_site=read_start
	softclip_site=0

	read_end=read_start+int(match_base)
	#print("%s\t%d\t%s\t%d\t%d\t%d"%(read_name,match_base,cigar_info,len(softclip_count),int(softclip_start),int(softclip_end)))
	if read_strand == "minus" and (not(mapping_type == "U")):
		softclip_seq=Complement_Reverse(softclip_seq)
		sequence=Complement_Reverse(sequence)
	softclip_human_region=softclip_region(softclip_start,softclip_end,read_strand)
	if mapping_type == "U":
		read_chr="*"
		read_start=-1
		read_end=-1
		read_strand="*"
		softclip_site=-1
		if read_flag & 16:
			softclip_seq=Complement_Reverse(sequence)
		else:
			softclip_seq=sequence
		softclip_human_region="*"
	softclip_len=len(softclip_seq)
	write_h.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(read_name,read_id,read_chr,read_strand,str(read_start),str(read_end),str(softclip_site),softclip_human_region,softclip_len,mapping_type,read_cigar,sequence,softclip_seq,multi_mapping_site))




OUT=open(outfile,'w')
OUT.write("read_name\tread_id\tchr\tstrand\tstart\tend\tsoftclip_site\tsoftclip_type\tsoftclip_len\tmapping_type\tcigar\tread_seq\tsoftclip_seq\tmulti_mapping\n")
IN=open(infile,'r')
for line in IN:
	if line.startswith('@'):
		continue

	if "chrUn_" in line:
		continue
	arr2=line.split("\t")
	multi_mapping_site1=re.findall('XA:Z:(chr.*);',line)
	if len(multi_mapping_site1)==0:
		multi_mapping_site1=re.findall('XA:Z:(chr.*);',line)
	multi_mapping_site2=re.findall('SA:Z:(chr.*);	',line)
	if len(multi_mapping_site2)==0:
		multi_mapping_site2=re.findall('SA:Z:(chr.*);',line)
	if len(multi_mapping_site1)>0:
		deal_line(line,"Yes",OUT,quality)
		site_list=multi_mapping_site1[0].split(";")
		line2=line.strip()
		arr=line2.split("\t")
		for k in site_list:
			tt=k.split(",")
			k_chr=tt[0]
			if "_" in k_chr:
				continue
			k_start=tt[1]
			k_start=tt[1][1:]
			k_cigar=tt[2]
			k_strand=tt[1][0:1]
			if (k_strand=="+"):
				arr[1]=str(0)
			if (k_strand=="+") & (int(arr[1]) & 64):
				arr[1]=str(64)
			if (k_strand=="+") & (~(int(arr[1]) & 64)):
				arr[1]=str(128)
			

			if (k_strand=="-"):
				arr[1]=str(16)
			if (k_strand=="-") & (int(arr[1]) & 64):
				arr[1]=str(80)
			if (k_strand=="-") & (~(int(arr[1]) & 64)):
				arr[1]=str(144)
			
			arr[2]=k_chr
			arr[3]=k_start
			arr[5]=k_cigar
			line2="\t".join(arr)
			deal_line(line2,"Yes",OUT,quality)
		continue

	if len(multi_mapping_site2)>0:
		deal_line(line,"Yes",OUT,quality)
		site_list=multi_mapping_site2[0].split(";")
		line2=line.strip()
		arr=line2.split("\t")
		for k in site_list:
			tt=k.split(",")
			k_chr=tt[0]
			k_start=tt[1]
			k_strand=tt[2]
			k_cigar=tt[3]
			if "_" in k_chr:
				continue
			if (k_strand=="+") & (int(arr[1]) & 64):
				arr[1]=str(64)
			if (k_strand=="+") & (~(int(arr[1]) & 64)):
				arr[1]=str(128)
			if (k_strand=="-") & (int(arr[1]) & 64):
				arr[1]=str(80)
			if (k_strand=="-") & (~(int(arr[1]) & 64)):
				arr[1]=str(144)
			arr[2]=k_chr
			arr[3]=k_start
			arr[5]=k_cigar
			line2="\t".join(arr)
			deal_line(line2,"Yes",OUT,quality)
		continue

	else:
		line2=line.strip()
		arr=line2.split("\t")
		if "_" in arr[2] and "chr" in arr[2]:
			continue
		deal_line(line,"No",OUT,quality)	







	
IN.close()
OUT.close()


