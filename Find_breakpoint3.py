import os,sys
import argparse
from optparse import OptionParser
import numpy as np 
import pandas as pd
import re 
sys.setrecursionlimit(100000)

usage=""
parser=OptionParser(usage=usage)
parser.add_option('-d',"--outdir",dest="outdir",help="Output dir")
parser.add_option('-s',"--insert_size",dest="insert_size",help="insert size")
parser.add_option('-t',"--threshold",dest="threshold",help="threshold")
parser.add_option('-l',"--len",dest="read_len",help="read_len")
parser.add_option('-v',"--vlen",dest="virus_len",help="virus_len")
parser.add_option("-S","--filter-insertsize",dest="filter_insert",help="Activate the filtering of DNA fragment insert size (False or True)")
parser.add_option("-F","--chr-incompatibility",dest="chr_incompatibility",help="Activate the filtering of chromosomal incompatibility (False or True)")

(options,args)=parser.parse_args()
outdir=options.outdir
insert_size=options.insert_size
insert_size=int(float(insert_size))
threshold=int(float(options.threshold))
read_len=int(float(options.read_len))
virus_len=int(float(options.virus_len))
filter_size=options.filter_insert
chr_incompatibility=options.chr_incompatibility



def find_mode(x):
	a=pd.Series.mode(x)[0]
	return a

############################################### define integration strand , site and breakpoint ###############################################


def define_breakpoint(df):
	df['breakpoint_human']=[0]*len(df.index)
	df['breaktype']=['0']*len(df.index)
	df.index=range(len(df.index))
	df['softclip_type_virus']=df['softclip_type_virus'].astype("str")
	for index,row in df.iterrows():

		if row['softclip_type_virus']=="forward":
			df.loc[index,"softclip_type_human"]="backward"
			row['softclip_type_human']="backward"
		if row['softclip_type_virus']=="backward":
			df.loc[index,"softclip_type_human"]="forward"
			row['softclip_type_human']="forward"

		if row['softclip_type_virus']=='0' and row['softclip_type_human']=="forward":
			row['softclip_type_virus']="backward"
			df.loc[index,"softclip_type_virus"]="backward"
		if row['softclip_type_virus']=='0' and row['softclip_type_human']=="backward":
			row['softclip_type_virus']="forward"
			df.loc[index,"softclip_type_virus"]="forward"

		if row['softclip_type_human']=="forward" and row['strand_human']=="plus":
			df.loc[index,"breaktype"]="downstream"
			if df.loc[index,"breakpoint_human"]==0:
				df.loc[index,"breakpoint_human"]=row["start_human"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="minus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="plus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]

		if row['softclip_type_human']=="backward" and row['strand_human']=="minus":
			df.loc[index,"breaktype"]="downstream"
			if df.loc[index,"breakpoint_human"]==0:
				df.loc[index,"breakpoint_human"]=row["start_human"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="plus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="minus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]
					
		if row['softclip_type_human']=="forward" and row['strand_human']=="minus":
			df.loc[index,"breaktype"]="upstream"	
			if df.loc[index,"breakpoint_human"]==0:
				df.loc[index,"breakpoint_human"]=row["end_human"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="plus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="minus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]

		if row['softclip_type_human']=="backward" and row['strand_human']=="plus":
			df.loc[index,"breaktype"]="upstream"	
			if df.loc[index,"breakpoint_human"]==0:
				df.loc[index,"breakpoint_human"]=row["end_human"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="plus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
			if row["softclip_site_virus"]==0 and row['strand_virus']=="minus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]

	df=df[df['breakpoint_human']>0]
	return df




def define_breakpoint_type2(df):
	df['breakpoint_human']=[0]*len(df.index)
	df['breaktype']=['0']*len(df.index)
	df.index=range(len(df.index))
	df['softclip_type_virus']=df['softclip_type_virus'].astype("str")
	for index,row in df.iterrows():
		if row['softclip_type_virus']=="forward":
			df.loc[index,"softclip_type_human"]="backward"
			row['softclip_type_human']="backward"
		if row['softclip_type_virus']=="backward":
			df.loc[index,"softclip_type_human"]="forward"
			row['softclip_type_human']="forward"

		if row['strand_human']=="plus" and df.loc[index,"softclip_type_human"]=="forward":
			df.loc[index,"breaktype"]="downstream"
			df.loc[index,"breakpoint_human"]=row["start_human"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="plus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="minus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
		if row['strand_human']=="plus" and df.loc[index,"softclip_type_human"]=="backward":
			df.loc[index,"breaktype"]="upstream"
			df.loc[index,"breakpoint_human"]=row["end_human"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="plus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="minus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]



		if row['strand_human']=="minus" and df.loc[index,"softclip_type_human"]=="forward":
			df.loc[index,"breaktype"]="uptream"
			df.loc[index,"breakpoint_human"]=row["end_human"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="plus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="minus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
		if row['strand_human']=="minus" and df.loc[index,"softclip_type_human"]=="backward":
			df.loc[index,"breaktype"]="downstream"
			df.loc[index,"breakpoint_human"]=row["start_human"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="plus":
				df.loc[index,"softclip_site_virus"]=row["start_virus"]
			if row["softclip_site_virus"]==0 and row["strand_virus"]=="minus":
				df.loc[index,"softclip_site_virus"]=row["end_virus"]



	df=df[df['breakpoint_human']>0]
	return df



######################################################### cluster chimeric distance ################################


def calculate_distance(df,insert_size,result,read_len,virus_len,depth=0):

	depth=depth+1
	read_len=int(read_len)
	start_human=df["start_human"]
	start_virus=df["start_virus"]
	diff_value=[abs(y - x) for x,y in zip(start_human,start_human[1:])]
	diff_value=np.array(diff_value)
	mark=diff_value<int(1.5*insert_size)
	if sum(mark)<len(mark):
		reads=np.where(mark == False)[0][0]
		reads=reads+1
	if sum(mark)==len(mark):
		reads=len(mark)+1
	
	human_list=start_human[0:reads].values.tolist()
	virus_list=start_virus[0:reads].values.tolist()
	virus_list2=virus_list.copy()
	virus_list2=sorted(virus_list2)


	human_max=max(df.iloc[0:reads]["end_human"])
	human_min=min(df.iloc[0:reads]["start_human"])

	virus_max=max(df[0:reads]["end_virus"])
	virus_min=min(df[0:reads]["start_virus"])
	if (virus_min<read_len):
		virus_min=read_len
	if(virus_max> int(virus_len-read_len)):
		virus_max=int(virus_len-read_len)


	if reads >1:
		diff_value_virus=[abs(y - x) for x,y in zip(virus_list2,virus_list2[1:])]
		if max(diff_value_virus) > int(1*insert_size):
			diff_value_virus=[abs(y - x) for x,y in zip(virus_list,virus_list[1:])]
			diff_value_virus=np.array(diff_value_virus)
			diff_value_virus_dat=pd.DataFrame()
			diff_value_virus_dat["distance"]=diff_value_virus
			mark_virus=diff_value_virus_dat['distance']<int(1.5*insert_size)
			if sum(mark_virus)<len(mark_virus):
				reads_virus=np.where(mark_virus == False)[0][0]
				reads_virus=reads_virus
			if sum(mark_virus)==len(mark_virus):
				reads_virus=len(mark_virus)
			mark_virus=np.array(mark_virus)
			mark_virus=np.insert(mark_virus,reads_virus,"True").tolist()
			tmp=pd.DataFrame(columns=df.columns)
			tmp=df.iloc[0:reads_virus+1,:]
			tmp=tmp.assign(tmp=mark_virus[0:int(reads_virus+1)])
			tmp=tmp[tmp["tmp"]==True]
			tmp=tmp.drop("tmp", axis=1)
			human_max=max(tmp["end_human"])
			human_min=min(tmp["start_human"])

			virus_max=max(tmp["end_virus"])
			virus_min=min(tmp["start_virus"])

			if (virus_min<read_len):
				virus_min=read_len
			if(virus_max> int(virus_len-read_len)):
				virus_max=int(virus_len-read_len)
			if len(tmp.index)>0:
				dd={"chr_human":tmp.iloc[0]["chr_human"],"start_human":human_min-read_len,"end_human":human_max+read_len,"chr_virus":tmp.iloc[0]["chr_virus"],"start_virus":virus_min-read_len,"end_virus":virus_max+read_len,"cluster":depth,"supported":reads_virus}
				dd=pd.DataFrame(data=dd,index=[0])
				result=pd.concat([result,dd],ignore_index=True)
			df=df.iloc[reads_virus+1:,:]
			
			if len(df.index)==0:
				return result
			else:
				return calculate_distance(df,insert_size,result,read_len,virus_len,depth)
		else:
			dd={"chr_human":df.iloc[0]["chr_human"],"start_human":human_min-read_len,"end_human":human_max+read_len,"chr_virus":df.iloc[0]["chr_virus"],"start_virus":virus_min-read_len,"end_virus":virus_max+read_len,"cluster":depth,"supported":reads}
			dd=pd.DataFrame(data=dd,index=[0])
			result=pd.concat([result,dd],ignore_index=True)
			df=df.iloc[reads:,:]
			if len(df.index)==0:
				return result
			else:
				return calculate_distance(df,insert_size,result,read_len,virus_len,depth)

	else:
		if len(start_human)>0:
			dd={"chr_human":df.iloc[0]["chr_human"],"start_human":human_min-read_len,"end_human":human_max+read_len,"chr_virus":df.iloc[0]["chr_virus"],"start_virus":virus_min-read_len,"end_virus":virus_max+read_len,"cluster":depth,"supported":reads}
			dd=pd.DataFrame(data=dd,index=[0])
			result=pd.concat([result,dd],ignore_index=True)
		start_human=start_human[reads:]
		df=df.iloc[reads:,:]
		if len(start_human)==0:
			return result
		else:
			return calculate_distance(df,insert_size,result,read_len,virus_len,depth)

################################ cluster chimeric reads ################################


def Cluster(df,insert_size,read_len,virus_len):

	chrlist=pd.unique(df["chr_human"])
	chrlist=chrlist[~pd.Series(chrlist).str.contains(pat='chrUn|_'),]
	result_cluster=pd.DataFrame(columns=["chr_human","start_human","end_human","chr_virus","start_virus","end_virus","cluster","supported"])
	for chri in chrlist:
		datatmp=df[df["chr_human"]==chri]
		datatmp=datatmp.sort_values(by=["start_human"])
		result_tmp=pd.DataFrame(columns=["chr_human","start_human","end_human","chr_virus","start_virus","end_virus","cluster","supported"])
		result_distance=calculate_distance(datatmp,insert_size,result_tmp,read_len,virus_len)
		result_cluster=pd.concat([result_cluster,result_distance],ignore_index=True)
	return result_cluster


def get_repeat(dat):
	if ((sum(dat=="Yes"))/len(dat))>0.5:
		return "Yes"
	if (sum(dat=="No"))==len(dat):
		return "No"
	else:
		return "Possible"

################################ calculate all types distance ################################



def calculate_distance_result(datatmp,result,cutoff,depth=0):
	depth=depth+1
	datatmp.index=range(len(datatmp.index))
	site=datatmp["breakpoint_human"]
	softclip_reads=datatmp["support_reads_softclip"]
	virus_site=datatmp["breakpoint_virus"]
	repeat_mark=datatmp["multi_mapping_human"]
	repeat_region=datatmp["repeat_region"]
	diff_value=[y - x for x,y in zip(site,site[1:])]
	diff_value=np.array(diff_value)
	#print(diff_value)
	mark=diff_value<=10
	if sum(mark)<len(mark):
		reads=np.where(mark == False)[0][0]
		reads=reads+1
	if sum(mark)==len(mark):
		reads=len(mark)+1
	#print(reads)
	tmp_human=datatmp.iloc[0:reads,]
	tmp_human2=datatmp.iloc[reads:,]
	tmp_human2.index=range(len(tmp_human2.index))
	virus_site2=sorted(virus_site[0:reads])
	if reads > 1:
		if max([abs(y - x) for x,y in zip(virus_site2[0:reads],virus_site2[1:reads])])>10:    ################## need to cluster virus reads
			diff_value_virus=[abs(y - x) for x,y in zip(virus_site2[0:reads],virus_site2[1:reads])]
			diff_value_virus=np.array(diff_value_virus)
			mark_virus=diff_value_virus<=10
			if sum(mark_virus)<len(mark_virus):
				reads_virus=np.where(mark_virus == False)[0][0]
				reads_virus=reads_virus+1
			if sum(mark_virus)==len(mark_virus):
				reads_virus=len(mark_virus)+1
			tmp_virus=tmp_human.loc[tmp_human['breakpoint_virus'].isin(virus_site2[0:reads_virus])]
			tmp_virus.index = range(len(tmp_virus.index))
			if tmp_virus["support_reads_softclip"].astype(int).values.sum()>=cutoff:
				dd={'chr_human':tmp_virus.loc[0,"chr_human"],'breakpoint_human':int(tmp_virus["breakpoint_human"].astype("int").values.mean()),'chr_virus':tmp_virus.loc[0,"chr_virus"],'breakpoint_virus':int(tmp_virus["breakpoint_virus"].astype("int").values.mean()),'breaktype':tmp_virus.loc[0,"breaktype"],'virus_direction':tmp_virus.loc[0,"virus_direction"], 'support_reads_softclip':tmp_virus["support_reads_softclip"].astype("int").values.sum(),'repeat_region':";".join(tmp_virus["repeat_region"]),'repeat_mark':get_repeat(np.array(tmp_virus["multi_mapping_human"]))}
				dd=pd.DataFrame(data=dd,index=[0])
				result=pd.concat([result,dd],ignore_index=True)
				tmp_human=tmp_human.drop(tmp_human.index[tmp_human['breakpoint_virus'].isin(virus_site2[0:reads_virus])])
				datatmp=pd.concat([tmp_human,tmp_human2],ignore_index=True)
				datatmp.index = range(len(datatmp.index))
				if len(datatmp.index)>1:
					return calculate_distance_result(datatmp,result,depth,cutoff)
				else:
					return result
			else:
				tmp_human=tmp_human.drop(tmp_human.index[tmp_human['breakpoint_virus'].isin(virus_site2[0:reads_virus])])
				datatmp=pd.concat([tmp_human,tmp_human2],ignore_index=True)
				datatmp.index = range(len(datatmp.index))
			if datatmp.empty:
				return result
			else:
				return calculate_distance_result(datatmp,result,depth,cutoff)
		else:
			dd={'chr_human':datatmp.loc[0,"chr_human"],'breakpoint_human':int(site[0:reads].values.mean()),'chr_virus':datatmp.loc[0,"chr_virus"],'breakpoint_virus':int(virus_site[0:reads].values.mean()),'breaktype':datatmp.loc[0,"breaktype"],'virus_direction':datatmp.loc[0,"virus_direction"], 'support_reads_softclip':softclip_reads[0:reads].values.sum(),'repeat_region':";".join(repeat_region[0:reads]),'repeat_mark':get_repeat(repeat_mark[0:reads])}
			dd=pd.DataFrame(data=dd,index=[0])
			result=pd.concat([result,dd],ignore_index=True)
			datatmp=datatmp.iloc[reads:,]
			datatmp.index = range(len(datatmp.index))
			if len(datatmp.index)==0:
				return result
			else:
				return calculate_distance_result(datatmp,result,depth,cutoff)

	if reads==1 and len(datatmp.index)>1:
		dd={'chr_human':datatmp.loc[0,"chr_human"],'breakpoint_human':int(site[0:reads].values.mean()),'chr_virus':datatmp.loc[0,"chr_virus"],'breakpoint_virus':int(virus_site[0:reads].values.mean()),'breaktype':datatmp.loc[0,"breaktype"],'virus_direction':datatmp.loc[0,"virus_direction"], 'support_reads_softclip':softclip_reads[0:reads].values.sum(),'repeat_region':";".join(repeat_region[0:reads]),'repeat_mark':get_repeat(repeat_mark[0:reads])}
		dd=pd.DataFrame(data=dd,index=[0])
		result=pd.concat([result,dd],ignore_index=True)
		datatmp=datatmp.iloc[reads:,]
		datatmp.index = range(len(datatmp.index))
		if len(datatmp.index)==0:
				return result
		else:
			return calculate_distance_result(datatmp,result,depth,cutoff)

	if reads==1 and len(datatmp.index)==1:
		dd={'chr_human':datatmp.loc[0,"chr_human"],'breakpoint_human':int(site[0:reads].values.mean()),'chr_virus':datatmp.loc[0,"chr_virus"],'breakpoint_virus':int(virus_site[0:reads].values.mean()),'breaktype':datatmp.loc[0,"breaktype"],'virus_direction':datatmp.loc[0,"virus_direction"], 'support_reads_softclip':softclip_reads[0:reads].values.sum(),'repeat_region':";".join(repeat_region[0:reads]),'repeat_mark':get_repeat(repeat_mark[0:reads])}
		dd=pd.DataFrame(data=dd,index=[0])
		result=pd.concat([result,dd],ignore_index=True)
		return result

	datatmp=datatmp.iloc[reads:,]
	datatmp.index = range(len(datatmp.index))
	if len(datatmp.index)==0:
		return result
	else:
		dd=pd.DataFrame(data=dd,index=[0])
		result=pd.concat([result,dd],ignore_index=True)
		return calculate_distance_result(datatmp,result,depth,cutoff)



################################ cluster all types result ################################

def Cluster_result(key_uniq,dataframe,cutoff):   
	result=pd.DataFrame(columns=dataframe.columns)
	for k in key_uniq:
		datatmp=dataframe[dataframe["key"]==k]
		datatmp=datatmp.sort_values(by=["chr_human","breakpoint_human","breakpoint_virus"])
		result_cluster_up=pd.DataFrame(columns=dataframe.columns)
		result_cluster_down=pd.DataFrame(columns=dataframe.columns)
		datatmp_up=datatmp[(datatmp["breaktype"]=="upstream")]
		datatmp_down=datatmp[(datatmp["breaktype"]=="downstream")]
		distance_result_up=calculate_distance_result(datatmp_up,result_cluster_up,cutoff)
		distance_result_down=calculate_distance_result(datatmp_down,result_cluster_down,cutoff)
		distance_result=pd.concat([distance_result_up,distance_result_down],ignore_index=True)
		result=pd.concat([result,distance_result],ignore_index=True)
	return result



################################################################  ###############################################################################################################


HBV_Full_data=pd.read_table("%s/HBV_Full.out"%(outdir))
HBV_Partial_data=pd.read_table("%s/HBV_Partial.out"%(outdir))
HBV_Unmapped_data=pd.read_table("%s/HBV_Unmapped.out"%(outdir))

Partial_HBV_Full_human_data=pd.read_table("%s/Partial_HBV_Full.out"%(outdir))    ##### reads partially mapped to HBV and softclip sequence fully mapped to human genome
Partial_HBV_Partial_human_data=pd.read_table("%s/Partial_HBV_Partial.out"%(outdir)) ##### reads partially mapped to HBV and softclip sequence partially mapped to human genome
Partial_HBV_Unmapped_human_data=pd.read_table("%s/Partial_HBV_Unmapped.out"%(outdir)) ##### reads partially mapped to HBV and softclip sequence does not mapped to human genome


################################################################ deal with repeat regions (extract multiple sites) ###############################################################################################################


Unmap_HBV_Full_human_data=pd.read_table("%s/Unmap_HBV_Full.out"%(outdir))
two_end_full=Unmap_HBV_Full_human_data[Unmap_HBV_Full_human_data["multi_mapping"]=="No"].groupby(["read_id"]).size().reset_index(name="read_count")
two_end_full=two_end_full[two_end_full["read_count"]>1]["read_id"]
Unmap_HBV_Full_human_data=Unmap_HBV_Full_human_data[~(Unmap_HBV_Full_human_data["read_id"].isin(two_end_full))]    ######## extract uniq mapped reads
Unmap_HBV_Full_human_data.index=range(len(Unmap_HBV_Full_human_data.index))

Unmap_HBV_Partial_human_data=pd.read_table("%s/Unmap_HBV_Partial.out"%(outdir))
Unmap_HBV_Unmapped_human_data=pd.read_table("%s/Unmap_HBV_Unmapped.out"%(outdir))



############################################################## merge human part of HBV partial result ###############################################################

HBV_Partial_mapping_data=pd.concat([Partial_HBV_Full_human_data,Partial_HBV_Partial_human_data])
HBV_Partial_unmapping_data=HBV_Partial_data[~(HBV_Partial_data["read_name"].isin(HBV_Partial_mapping_data["read_name"]))]

############################################################# type1A (one read softclip read, another read fully map to human genome) #############################################################

type1A=HBV_Partial_mapping_data.merge(HBV_Partial_data,on="read_name",suffixes=('_human', '_virus'))
type1A["read_id"]=type1A["read_id_virus"]
type1A=type1A.merge(Unmap_HBV_Full_human_data,on="read_id")
type1A=type1A[~(type1A["read_name_x"]==type1A["read_name_y"])]
if(chr_incompatibility!="False"):
	type1A=type1A[type1A["chr_human"]==type1A["chr"]]
if(filter_size!="False"):
	type1A=type1A[abs(type1A["start_human"]-type1A["start"])<int(1.5*insert_size)]   #### sites mapping to human genome are close within 1.5 insert size
type1A=type1A[(abs(type1A["start_human"]-type1A["end_human"]+type1A["start_virus"]-type1A["end_virus"])<read_len+int(0.1*read_len))&(abs(type1A["start_human"]-type1A["end_human"]+type1A["start_virus"]-type1A["end_virus"])>int(0.9*read_len))] #########  total length (virus+human) no more than 1.2 read len
type1A=type1A[~(type1A["strand_human"]==type1A["strand"])]
type1A_result=define_breakpoint(type1A)
type1A_result.to_csv("%s/type1A.out"%(outdir),sep='\t',index=False)

############################################################# type1B (one read softclip read （human bwa, virus bwa）, another read partially map to human genome) #############################################################

#####check reads (cannot align to virus genome by bwa) whehter align to virus genome by blastn
if os.path.isfile("%s/Unmap_HBV_Partial_HBV_blastn.out"%(outdir)):
	Unmap_HBV_Partial_HBV_blastn=pd.read_table("%s/Unmap_HBV_Partial_HBV_blastn.out"%(outdir),header=None,names=("softclip_name","softclip_chr","softclip_ident","softclip_len","softclip_mis","softclip_gap","softclip_s","softclip_e","softclip_start","softclip_end","softclip_evalue","softclip_score","softclip_strand"))
else:
	Unmap_HBV_Partial_HBV_blastn=pd.DataFrame(columns=["softclip_name","softclip_chr","softclip_ident","softclip_len","softclip_mis","softclip_gap","softclip_s","softclip_e","softclip_start","softclip_end","softclip_evalue","softclip_score","softclip_strand"])
Unmap_HBV_Partial_HBV_blastn=Unmap_HBV_Partial_HBV_blastn.drop_duplicates(subset=["softclip_name"],keep="first")
Unmap_HBV_Partial_HBV_blastn["key"]=Unmap_HBV_Partial_HBV_blastn["softclip_name"].map(lambda x:x[0:len(x)-6])
Unmap_HBV_Partial_HBV_blastn=Unmap_HBV_Partial_HBV_blastn[Unmap_HBV_Partial_HBV_blastn["softclip_ident"]>80]
Unmap_HBV_Partial_HBV_blastn=Unmap_HBV_Partial_HBV_blastn.rename(columns={"softclip_name":"read_name","key":"read_id","softclip_chr":"chr","softclip_start":"start","softclip_end":"end","softclip_strand":"strand"})
Unmap_HBV_Partial_HBV_blastn=Unmap_HBV_Partial_HBV_blastn.assign(softclip_site=[0]*len(Unmap_HBV_Partial_HBV_blastn.index),softclip_type=['0']*len(Unmap_HBV_Partial_HBV_blastn.index),cigar=['0']*len(Unmap_HBV_Partial_HBV_blastn.index),read_seq=['0']*len(Unmap_HBV_Partial_HBV_blastn.index),softclip_seq=['0']*len(Unmap_HBV_Partial_HBV_blastn.index),mapping_type=['F']*len(Unmap_HBV_Partial_HBV_blastn.index),multi_mapping=['No']*len(Unmap_HBV_Partial_HBV_blastn.index))
Unmap_HBV_Partial_HBV_blastn=Unmap_HBV_Partial_HBV_blastn[["read_name","read_id","chr","strand","start","end","softclip_site","softclip_type","softclip_len","mapping_type","cigar","read_seq","softclip_seq","multi_mapping"]]
Unmap_HBV_Partial_HBV_blastn["chr"]=Unmap_HBV_Partial_HBV_blastn["chr"].map(lambda x:x.split("|")[1])
Unmap_HBV_Partial_HBV_blastn2=Unmap_HBV_Partial_HBV_blastn.copy()
Unmap_HBV_Partial_HBV_blastn.loc[Unmap_HBV_Partial_HBV_blastn["strand"]=="minus","start"]=Unmap_HBV_Partial_HBV_blastn2.loc[Unmap_HBV_Partial_HBV_blastn2["strand"]=="minus","end"]
Unmap_HBV_Partial_HBV_blastn.loc[Unmap_HBV_Partial_HBV_blastn["strand"]=="minus","end"]=Unmap_HBV_Partial_HBV_blastn2.loc[Unmap_HBV_Partial_HBV_blastn2["strand"]=="minus","start"]


type1B=HBV_Partial_mapping_data.merge(HBV_Partial_data,on="read_name",suffixes=('_human', '_virus'))
type1B["read_id"]=type1B["read_id_virus"]

#################### remove reads aligned to human but exist in blastn result ##########
Unmap_HBV_Partial_HBV_unmap=Unmap_HBV_Partial_human_data[~(Unmap_HBV_Partial_human_data["read_id"].isin(Unmap_HBV_Partial_HBV_blastn["read_id"]))]

type1B=type1B.merge(Unmap_HBV_Partial_HBV_unmap,on="read_id")  

type1B=type1B[type1B["read_id"].isin(Unmap_HBV_Partial_HBV_unmap)]

type1B=type1B[~(type1B["read_name_x"]==type1B["read_name_y"])]

type1B=type1B[type1B["chr_human"]==type1B["chr"]]
if(filter_size!="False"):
	type1B=type1B[abs(type1B["start_human"]-type1B["start"])<int(1.5*insert_size)]
type1B=type1B[(abs(type1B["start_human"]-type1B["end_human"]+type1B["start_virus"]-type1B["end_virus"])<read_len+int(0.1*read_len))&(abs(type1B["start_human"]-type1B["end_human"]+type1B["start_virus"]-type1B["end_virus"])>int(0.9*read_len))]

type1B=type1B[~(type1B["strand_human"]==type1B["strand"])]

type1B_result=define_breakpoint(type1B)
type1B_result.to_csv("%s/type1B.out"%(outdir),sep='\t',index=False)

################################################################# type1C (one read softclip read (human by bwa, virus by blastn), another read fully map to human ) ###########################################################

type1C=Unmap_HBV_Partial_human_data.merge(Unmap_HBV_Partial_HBV_blastn,on="read_name",suffixes=('_human', '_virus'))
type1C["read_id"]=type1C["read_id_virus"]
type1C=type1C.merge(Unmap_HBV_Full_human_data,on="read_id")
type1C=type1C[~(type1C["read_name_x"]==type1C["read_name_y"])]
type1C=type1C[type1C["chr_human"]==type1C["chr"]]
if(filter_size!="False"):
	type1C=type1C[abs(type1C["start_human"]-type1C["start"])<int(1.5*insert_size)]
type1C=type1C[(abs(type1C["start_human"]-type1C["end_human"]+type1C["start_virus"]-type1C["end_virus"])<read_len+int(0.1*read_len))&(abs(type1C["start_human"]-type1C["end_human"]+type1C["start_virus"]-type1C["end_virus"])>int(0.9*read_len))]
type1C=type1C[~(type1C["strand_human"]==type1C["strand"])]
type1C_result=define_breakpoint(type1C)
type1C_result.to_csv("%s/type1C.out"%(outdir),sep='\t',index=False)

############################################################### type1D (one softclip read (human by bwa, virus by blastn), another read partially map to human genome ################################################################

type1D=Unmap_HBV_Partial_human_data.merge(Unmap_HBV_Partial_HBV_blastn,on="read_name",suffixes=('_human', '_virus'))
type1D["read_id"]=type1D["read_id_virus"]
type1D=type1D.merge(Unmap_HBV_Partial_human_data,on="read_id")
type1D=type1D[~(type1D["read_name_x"]==type1D["read_name_y"])]
type1D=type1D[type1D["chr_human"]==type1D["chr"]]
if(filter_size!="False"):
	type1D=type1D[abs(type1D["start_human"]-type1D["start"])<int(1.5*insert_size)]
softclip_read_count=type1D.groupby(["read_id"]).size().reset_index(name="read_count")
softclip_read_count_type1=softclip_read_count[softclip_read_count["read_count"]==1]["read_id"]  
softclip_read_count_type4=softclip_read_count[softclip_read_count["read_count"]>1]["read_id"]  ######### remove type4 from type1D, type4 double softclip
type4B=type1D[type1D["read_id"].isin(softclip_read_count_type4)]
type1D=type1D[type1D["read_id"].isin(softclip_read_count_type1)]
type1D=type1D[(abs(type1D["start_human"]-type1D["end_human"]+type1D["start_virus"]-type1D["end_virus"])<read_len+int(0.1*read_len))&(abs(type1D["start_human"]-type1D["end_human"]+type1D["start_virus"]-type1D["end_virus"])>int(0.9*read_len))]
type1D=type1D[~(type1D["strand_human"]==type1D["strand"])]
type1D_result=define_breakpoint(type1D)
type1D_result.to_csv("%s/type1D.out"%(outdir),sep='\t',index=False)

type1_result=pd.concat([type1A_result,type1B_result,type1C_result,type1D_result])
type1_result=type1_result.assign(type=["type1"]*len(type1_result.index))


########################################################## type4 #############################################################################


########################################################## type4B, double softclip reads, human by bwa, virus by blastn #############################################################################

type4B=type4B[(abs(type4B["start_human"]-type4B["end_human"]+type4B["start_virus"]-type4B["end_virus"])<read_len+int(0.1*read_len))&(abs(type4B["start_human"]-type4B["end_human"]+type4B["start_virus"]-type4B["end_virus"])>int(0.9*read_len))]
type4B=type4B[~(type4B["strand_human"]==type4B["strand"])]
type4B_result=define_breakpoint(type4B)
type4B_result.to_csv("%s/type4B.out"%(outdir),sep='\t',index=False)


########################################################## type4B, double softclip reads, human by bwa, virus by bwa #############################################################################


type4A_2=HBV_Partial_mapping_data.merge(HBV_Partial_data,on="read_name",suffixes=('_human', '_virus'))
type4A_2_raw=type4A_2.copy()
type4A_2_raw=type4A_2_raw[['read_name','read_id_human', 'chr_human', 'strand_human','start_human', 'end_human', 'softclip_site_human','softclip_type_human', 'softclip_len_human', 'mapping_type_human','cigar_human', 'read_seq_human', 'softclip_seq_human','multi_mapping_human']]
type4A_2_raw=type4A_2_raw.rename(columns={'read_id_human':'read_id', 'chr_human':'chr', 'strand_human':'strand','start_human':'start', 'end_human':'end', 'softclip_site_human':'softclip_site','softclip_type_human':'softclip_type', 'softclip_len_human':'softclip_len', 'mapping_type_human':'mapping_type','cigar_human':'cigar', 'read_seq_human':'read_seq', 'softclip_seq_human':'softclip_seq','multi_mapping_human':'multi_mapping'})
type4A_2["read_id"]=type4A_2["read_id_virus"]
type4A_2=type4A_2.merge(type4A_2_raw,on="read_id")
type4A_2=type4A_2[~(type4A_2["read_name_x"]==type4A_2["read_name_y"])]


if (chr_incompatibility=="True" and filter_size=="True"):
	type4A_2=type4A_2[type4A_2["chr_human"]==type4A_2["chr"]]
	type4A_2=type4A_2[abs(type4A_2["start_human"]-type4A_2["start"])<int(1.5*insert_size)]
	type4A=type4A_2.copy()
	type4_softclip_read_count=type4A.groupby(["read_id"]).size().reset_index(name="read_count")
	type4_softclip_read_count=type4_softclip_read_count[type4_softclip_read_count["read_count"]>1]["read_id"]
	type4A=type4A[type4A["read_id"].isin(type4_softclip_read_count)]
	type4A=type4A[(abs(type4A["start_human"]-type4A["end_human"]+type4A["start_virus"]-type4A["end_virus"])<read_len+int(0.1*read_len))&(abs(type4A["start_human"]-type4A["end_human"]+type4A["start_virus"]-type4A["end_virus"])>int(0.9*read_len))]
	type4A=type4A[~(type4A["strand_human"]==type4A["strand"])]
	type4A_result=define_breakpoint(type4A)

elif(filter_size=="False" and chr_incompatibility=="True"):
	type4A_2=type4A_2[type4A_2["chr_human"]==type4A_2["chr"]]
	type4A=type4A_2.copy()
	type4_softclip_read_count=type4A.groupby(["read_id"]).size().reset_index(name="read_count")
	type4_softclip_read_count=type4_softclip_read_count[type4_softclip_read_count["read_count"]>1]["read_id"]
	type4A=type4A[type4A["read_id"].isin(type4_softclip_read_count)]
	type4A=type4A[(abs(type4A["start_human"]-type4A["end_human"]+type4A["start_virus"]-type4A["end_virus"])<read_len+int(0.1*read_len))&(abs(type4A["start_human"]-type4A["end_human"]+type4A["start_virus"]-type4A["end_virus"])>int(0.9*read_len))]
	type4A=type4A[~(type4A["strand_human"]==type4A["strand"])]
	type4A_result=define_breakpoint(type4A)

elif(chr_incompatibility=="False"):
	type4A=type4A_2.copy()
	type4_softclip_read_count=type4A.groupby(["read_id"]).size().reset_index(name="read_count")
	type4_softclip_read_count=type4_softclip_read_count[type4_softclip_read_count["read_count"]>1]["read_id"]
	type4A=type4A[type4A["read_id"].isin(type4_softclip_read_count)]
	type4A=type4A[(abs(type4A["start_human"]-type4A["end_human"]+type4A["start_virus"]-type4A["end_virus"])<read_len+int(0.1*read_len))&(abs(type4A["start_human"]-type4A["end_human"]+type4A["start_virus"]-type4A["end_virus"])>int(0.9*read_len))]
	type4A_result=define_breakpoint(type4A)






type4A_result.to_csv("%s/type4A.out"%(outdir),sep='\t',index=False)
type4_result=pd.concat([type4A_result,type4B_result])
type4_result=type4_result.assign(type=["type4"]*len(type4_result.index))

##########################################################type3 cluster#####################################################




type3B=Unmap_HBV_Full_human_data.merge(HBV_Full_data,on="read_id",suffixes=('_human', '_virus'))
type3B=type3B[~(type3B["read_name_human"]==type3B["read_name_virus"])]
type3B=type3B.drop_duplicates(subset=["read_name_human"],keep=False)  ############## remove all multi-mapped reads for chimeric pairsfxxx
type3B=type3B.drop_duplicates(subset=["chr_human","start_human","end_human","chr_virus","start_virus","end_virus"],keep="first")
type3B=type3B[((type3B["multi_mapping_human"]=="No")&(type3B["multi_mapping_virus"]=="No"))]
type3_result=Cluster(type3B,insert_size,read_len,virus_len)
type3_result=type3_result.sort_values(by=["chr_human","start_human"])
type3_result.to_csv("%s/type3.out"%(outdir),sep='\t',index=False)


############################################################ type2 #####################################################

######################### type2B (one read softclip read (human by bwa, virus by blastn), another read fully map to virus  #########################
type2B=Unmap_HBV_Partial_human_data.merge(Unmap_HBV_Partial_HBV_blastn,on="read_name",suffixes=('_human', '_virus'))
type2B["read_id"]=type2B["read_id_virus"]
type2B=type2B.merge(HBV_Full_data,on="read_id")
type2B=type2B[~(type2B["read_name_x"]==type2B["read_name_y"])]
type2B=type2B[(abs(type2B["start_human"]-type2B["end_human"]+type2B["start_virus"]-type2B["end_virus"])<read_len+int(0.1*read_len))&(abs(type2B["start_human"]-type2B["end_human"]+type2B["start_virus"]-type2B["end_virus"])>int(0.9*read_len))]

type2B_result=define_breakpoint_type2(type2B)
type2B_result.to_csv("%s/type2B.out"%(outdir),sep='\t',index=False)

######################### type2A (one read softclip read (human by bwa, virus by blastn), another read partially map to virus  #########################

type2A=Unmap_HBV_Partial_human_data.merge(Unmap_HBV_Partial_HBV_blastn,on="read_name",suffixes=('_human', '_virus'))
type2A["read_id"]=type2A["read_id_virus"]
type2A=type2A.merge(HBV_Partial_unmapping_data,on="read_id")
type2A=type2A[~(type2A["read_name_x"]==type2A["read_name_y"])]
type2A=type2A[(abs(type2A["start_human"]-type2A["end_human"]+type2A["start_virus"]-type2A["end_virus"])<read_len+int(0.1*read_len))&(abs(type2A["start_human"]-type2A["end_human"]+type2A["start_virus"]-type2A["end_virus"])>int(0.9*read_len))]

type2A_result=define_breakpoint_type2(type2A)
type2A_result.to_csv("%s/type2A.out"%(outdir),sep='\t',index=False)

######################### type2A (one read softclip read (human by bwa, virus by bwa), another read Fully map to virus  #########################


type2C=HBV_Partial_mapping_data.merge(HBV_Partial_data,on="read_name",suffixes=('_human', '_virus'))
type2C["read_id"]=type2C["read_id_virus"]
type2C=type2C.merge(HBV_Full_data,on="read_id")
type2C=type2C[~(type2C["read_name_x"]==type2C["read_name_y"])]
type2C=type2C[(abs(type2C["start_human"]-type2C["end_human"]+type2C["start_virus"]-type2C["end_virus"])<read_len+int(0.1*read_len))&(abs(type2C["start_human"]-type2C["end_human"]+type2C["start_virus"]-type2C["end_virus"])>int(0.9*read_len))]
type2C.to_csv("%s/type2C.txt"%(outdir),sep='\t',index=False)
type2C_result=define_breakpoint_type2(type2C)
type2C_result.to_csv("%s/type2C.out"%(outdir),sep='\t',index=False)

type2_result=pd.concat([type2B_result,type2A_result,type2C_result])
type2_result=type2_result.assign(type=["type2"]*len(type2_result.index))


########################################################## integrate result ##################################################

type2_result=type2_result[~type2_result["read_id_human"].isin(type1_result["read_id_human"])]

type4_result=type4_result[~type4_result["read_id_human"].isin(type1_result["read_id_human"])]


all_type_result=pd.concat([type1_result,type2_result,type4_result])

########################### remove PCR duplicates ###########################
all_type_result=all_type_result.drop_duplicates(subset=["chr_human","start_human","end_human","chr_virus","start_virus","end_virus","chr","start","end"],keep="first")

########################### remove chimerics in BWA ###########################
all_type_result=all_type_result.drop_duplicates(subset=["read_name_x","chr_human","start_human","end_human","chr_virus","start_virus","end_virus"],keep="first")

all_type_result_nodup=all_type_result[(all_type_result["multi_mapping_human"]=="No")]
all_type_result_nodup=all_type_result_nodup.assign(repeat_region=["No"]*len(all_type_result_nodup.index))
all_type_result_dup=all_type_result[(all_type_result["multi_mapping_human"]=="Yes")]
all_type_result_dup=all_type_result_dup.sort_values(by=["read_name_x"])
all_type_result_dup_id=all_type_result_dup.drop_duplicates(subset=["read_name_x"],keep="first")
all_type_result_dup_count=all_type_result_dup.groupby("read_name_x").size().reset_index(name="read_count")
all_type_result_dup_count["id"]=range(0,len(all_type_result_dup_count.index))
all_type_result_dup_count["id"]=all_type_result_dup_count["id"]+1
all_type_result_dup["repeat_region"]=[item for item, count in zip(all_type_result_dup_count["id"],all_type_result_dup_count["read_count"]) for i in range(count)]
all_type_result_dup["repeat_region"]="repeat"+all_type_result_dup["repeat_region"].astype(str)

all_type_result=pd.concat([all_type_result_nodup,all_type_result_dup],ignore_index=True)

##############################################################################

all_type_result["softclip_site_virus"]=all_type_result["softclip_site_virus"].astype("int")
all_type_result["breakpoint_human"]=all_type_result["breakpoint_human"].astype("int")
all_type_result.to_csv("%s/all_breakpoints.out"%(outdir),sep='\t',index=False)
if len(all_type_result.index)>0:
	all_type_result["virus_direction"]=all_type_result.apply(lambda x: "positive" if x["strand_human"]==x["strand_virus"] else "negative",axis=1)
	all_type_result_support_read=all_type_result.groupby(["chr_human","virus_direction","breakpoint_human",'chr_virus','softclip_site_virus','breaktype','multi_mapping_human']).size().reset_index(name="support_reads_softclip")
	all_type_result_multi_site=all_type_result.groupby(["chr_human","virus_direction","breakpoint_human",'chr_virus','softclip_site_virus','breaktype','multi_mapping_human'],as_index=False)['repeat_region'].apply(';'.join)
	all_type_result_support_read=all_type_result_support_read.merge(all_type_result_multi_site,on=["chr_human","virus_direction","breakpoint_human",'chr_virus','softclip_site_virus','breaktype','multi_mapping_human'])
	all_type_result_support_read.rename(columns={'softclip_site_virus':'breakpoint_virus'},inplace=True)
	all_type_result_support_read.to_csv("%s/all_breakpoints_2.out"%(outdir),sep='\t',index=False)
else:
	all_type_result_support_read=pd.DataFrame(columns=all_type_result.columns)

################################################################ count reads for breakpoint #################################################################



################################# if breakpoints and cluster exist #################################

if (len(all_type_result_support_read.index)>0) and (len(type3_result.index)>0):
	all_type_result=all_type_result_support_read
	all_type_result["breakpoint_virus"]=all_type_result["breakpoint_virus"].astype("int")


############################################################ merge softclip and chimeric result ############################################################
	all_type_result["key"]=all_type_result[["chr_human","chr_virus","breaktype","virus_direction"]].apply(lambda x: '_'.join(x), axis=1)  ############# consider breaktype, virus_direction, when cluster breakpoint#########################
	key_uniq=pd.unique(all_type_result["key"])
	all_type_result.to_csv("%s/all_breakpoints_3.out"%(outdir),sep='\t',index=False)
	all_type_result_cluster=Cluster_result(key_uniq,all_type_result,threshold)
	all_type_result_cluster=all_type_result_cluster.assign(chr_human_chimeric=['0']*len(all_type_result_cluster.index),start_human_chimeric=[0]*len(all_type_result_cluster.index),end_human_chimeric=[0]*len(all_type_result_cluster.index),chr_virus_chimeric=["0"]*len(all_type_result_cluster.index),start_virus_chimeric=[0]*len(all_type_result_cluster.index),end_virus_chimeric=[0]*len(all_type_result_cluster.index),support_reads_chimeric=[0]*len(all_type_result_cluster.index),total_reads=[0]*len(all_type_result_cluster.index))
	all_type_result_cluster.index=range(len(all_type_result_cluster.index))
	chimeric_result=pd.DataFrame(columns=all_type_result_cluster.columns)
	for index_c,row_c in type3_result.iterrows():
		mark=0
		for index_s,row_s in all_type_result_cluster.iterrows():
			# if(filter_size!="False"):
			# 	if (row_s["chr_human"]==row_c["chr_human"]) and (abs(row_s["breakpoint_human"]-row_c["start_human"])<int(1.5*insert_size) or abs(row_s["breakpoint_human"]-row_c["end_human"])<int(1.5*insert_size)) and (abs(row_s["breakpoint_virus"]-row_c["start_virus"])<int(1.5*insert_size) or abs(row_s["breakpoint_virus"]-row_c["end_virus"])<int(1.5*insert_size)):
			# 		all_type_result_cluster.loc[index_s,"chr_human_chimeric"]=row_c["chr_human"]
			# 		all_type_result_cluster.loc[index_s,"start_human_chimeric"]=row_c["start_human"]
			# 		all_type_result_cluster.loc[index_s,"end_human_chimeric"]=row_c["end_human"]
			# 		all_type_result_cluster.loc[index_s,"chr_virus_chimeric"]=row_c["chr_virus"]
			# 		all_type_result_cluster.loc[index_s,"start_virus_chimeric"]=row_c["start_virus"]
			# 		all_type_result_cluster.loc[index_s,"end_virus_chimeric"]=row_c["end_virus"]
			# 		all_type_result_cluster.loc[index_s,"support_reads_chimeric"]=int(row_c["supported"])
			# 		mark=1
			# if(filter_size=="False"):
			# 	if (row_s["chr_human"]==row_c["chr_human"]):
			# 		all_type_result_cluster.loc[index_s,"chr_human_chimeric"]=row_c["chr_human"]
			# 		all_type_result_cluster.loc[index_s,"start_human_chimeric"]=row_c["start_human"]
			# 		all_type_result_cluster.loc[index_s,"end_human_chimeric"]=row_c["end_human"]
			# 		all_type_result_cluster.loc[index_s,"chr_virus_chimeric"]=row_c["chr_virus"]
			# 		all_type_result_cluster.loc[index_s,"start_virus_chimeric"]=row_c["start_virus"]
			# 		all_type_result_cluster.loc[index_s,"end_virus_chimeric"]=row_c["end_virus"]
			# 		all_type_result_cluster.loc[index_s,"support_reads_chimeric"]=int(row_c["supported"])
			# 		mark=1
			if (row_s["chr_human"]==row_c["chr_human"]) and (abs(row_s["breakpoint_human"]-row_c["start_human"])<int(1.5*insert_size) or abs(row_s["breakpoint_human"]-row_c["end_human"])<int(1.5*insert_size)) and (abs(row_s["breakpoint_virus"]-row_c["start_virus"])<int(1.5*insert_size) or abs(row_s["breakpoint_virus"]-row_c["end_virus"])<int(1.5*insert_size)):
					all_type_result_cluster.loc[index_s,"chr_human_chimeric"]=row_c["chr_human"]
					all_type_result_cluster.loc[index_s,"start_human_chimeric"]=row_c["start_human"]
					all_type_result_cluster.loc[index_s,"end_human_chimeric"]=row_c["end_human"]
					all_type_result_cluster.loc[index_s,"chr_virus_chimeric"]=row_c["chr_virus"]
					all_type_result_cluster.loc[index_s,"start_virus_chimeric"]=row_c["start_virus"]
					all_type_result_cluster.loc[index_s,"end_virus_chimeric"]=row_c["end_virus"]
					all_type_result_cluster.loc[index_s,"support_reads_chimeric"]=int(row_c["supported"])
					mark=1


		if mark==0:
			row2=all_type_result_cluster.iloc[[index_s]].copy()
			row2["chr_human"]='0'
			row2["breakpoint_human"]=0
			row2["chr_virus"]='0'
			row2["breakpoint_virus"]=0
			row2["support_reads_softclip"]=0
			row2["breaktype"]='0'
			row2["chr_human_chimeric"]=row_c["chr_human"]
			row2["start_human_chimeric"]=int(row_c["start_human"])
			row2["end_human_chimeric"]=int(row_c["end_human"])
			row2["chr_virus_chimeric"]=row_c["chr_virus"]
			row2["start_virus_chimeric"]=int(row_c["start_virus"])
			row2["end_virus_chimeric"]=int(row_c["end_virus"])
			row2["support_reads_chimeric"]=int(row_c["supported"])
			row2["total_reads"]=int(row_c["supported"])
			row2["virus_direction"]='0'
			row2["repeat_region"]='No'
			row2["repeat_mark"]='No'
			chimeric_result=pd.concat([chimeric_result,row2])


	
	all_type_result_cluster["total_reads"]=all_type_result_cluster.apply(lambda x:x["support_reads_softclip"]+x["support_reads_chimeric"],axis=1)
	all_type_result_cluster=all_type_result_cluster.sort_values(by=['chr_human','breakpoint_human'])
	all_type_result_cluster=all_type_result_cluster.drop(columns=["key"])
	all_type_result_cluster=pd.concat([all_type_result_cluster,chimeric_result],ignore_index=True)
	all_type_result_cluster=all_type_result_cluster[["chr_human","breakpoint_human","chr_virus","breakpoint_virus","breaktype","virus_direction","chr_human_chimeric","start_human_chimeric","end_human_chimeric","chr_virus_chimeric","start_virus_chimeric","end_virus_chimeric","support_reads_softclip","support_reads_chimeric","total_reads","repeat_region","repeat_mark"]]

################################# if breakpoints does not exist and clusters exist #################################

if (len(all_type_result.index)==0) and (len(type3_result.index)>0):
	type3_result.index=range(len(type3_result.index))
	all_type_result_cluster=pd.DataFrame(0,index=np.arange(len(type3_result.index)),columns=["chr_human","breakpoint_human","chr_virus","breakpoint_virus","breaktype","virus_direction","chr_human_chimeric","start_human_chimeric","end_human_chimeric","chr_virus_chimeric","start_virus_chimeric","end_virus_chimeric","support_reads_softclip","support_reads_chimeric","total_reads","repeat_region"])
	for index_c,row_c in all_type_result_cluster.iterrows():
		all_type_result_cluster.loc[index_c,"chr_human_chimeric"]=type3_result.loc[index_c,"chr_human"]
		all_type_result_cluster.loc[index_c,"start_human_chimeric"]=type3_result.loc[index_c,"start_human"]
		all_type_result_cluster.loc[index_c,"end_human_chimeric"]=type3_result.loc[index_c,"end_human"]
		all_type_result_cluster.loc[index_c,"chr_virus_chimeric"]=type3_result.loc[index_c,"chr_virus"]
		all_type_result_cluster.loc[index_c,"start_virus_chimeric"]=type3_result.loc[index_c,"start_virus"]
		all_type_result_cluster.loc[index_c,"end_virus_chimeric"]=type3_result.loc[index_c,"end_virus"]
		all_type_result_cluster.loc[index_c,"support_reads_chimeric"]=int(type3_result.loc[index_c,"supported"])
		all_type_result_cluster.loc[index_c,"total_reads"]=int(type3_result.loc[index_c,"supported"])
		all_type_result_cluster.loc[index_c,"repeat_region"]="No"
		all_type_result_cluster.loc[index_c,"repeat_mark"]="No"
	#all_type_result_cluster=all_type_result_cluster[(all_type_result_cluster["total_reads"]>=threshold)]

################################# if breakpoints  exist and cluster does not exist #################################
if (len(all_type_result.index)>0) and (len(type3_result.index)==0):
	all_type_result=all_type_result_support_read
	all_type_result["breakpoint_virus"]=all_type_result["breakpoint_virus"].astype("int")
	all_type_result["key"]=all_type_result[["chr_human","chr_virus","breaktype","virus_direction"]].apply(lambda x: '_'.join(x), axis=1)
	key_uniq=pd.unique(all_type_result["key"])
	
	
	######### this region is new ####################33
	
	all_type_result_cluster=Cluster_result(key_uniq,all_type_result,threshold)
	all_type_result_cluster=all_type_result_cluster.assign(chr_human_chimeric=['0']*len(all_type_result_cluster.index),start_human_chimeric=[0]*len(all_type_result_cluster.index),end_human_chimeric=[0]*len(all_type_result_cluster.index),chr_virus_chimeric=["0"]*len(all_type_result_cluster.index),start_virus_chimeric=[0]*len(all_type_result_cluster.index),end_virus_chimeric=[0]*len(all_type_result_cluster.index),support_reads_chimeric=[0]*len(all_type_result_cluster.index),total_reads=[0]*len(all_type_result_cluster.index))
	all_type_result_cluster.index=range(len(all_type_result_cluster.index))
	all_type_result_cluster["total_reads"]=all_type_result_cluster["support_reads_softclip"].astype(int)
	######### this region is new ####################33

	#all_type_result_cluster=all_type_result_cluster[(all_type_result_cluster["total_reads"]>=threshold)]
	all_type_result_cluster=all_type_result_cluster.sort_values(by=['chr_human','breakpoint_human'])
	all_type_result_cluster=all_type_result_cluster.drop(columns=["key"])
	all_type_result_cluster=all_type_result_cluster[["chr_human","breakpoint_human","chr_virus","breakpoint_virus","breaktype","virus_direction","chr_human_chimeric","start_human_chimeric","end_human_chimeric","chr_virus_chimeric","start_virus_chimeric","end_virus_chimeric","support_reads_softclip","support_reads_chimeric","total_reads","repeat_region","repeat_mark"]]

if (len(all_type_result.index)==0) and (len(type3_result.index)==0):
	all_type_result_cluster=pd.DataFrame(columns=["chr_human","breakpoint_human","chr_virus","breakpoint_virus","breaktype","virus_direction","chr_human_chimeric","start_human_chimeric","end_human_chimeric","chr_virus_chimeric","start_virus_chimeric","end_virus_chimeric","support_reads_softclip","support_reads_chimeric","total_reads","repeat_mark","region_human","gene_human","region_virus","gene_virus"])
	all_type_result_cluster.to_csv("%s/all_type_result.out"%(outdir),sep='\t',index=False)
	all_type_result_cluster_anno=pd.DataFrame()
	all_type_result_cluster_anno.to_csv(outdir+"/"+"human_annovar.avinput",index=False,sep='\t',header=False)

else:

	all_type_result_cluster=all_type_result_cluster[~((all_type_result_cluster["repeat_mark"]=="No")&(all_type_result_cluster["total_reads"]<threshold))]
	all_type_result_cluster.to_csv("%s/all_type_result.out"%(outdir),sep='\t',index=False)

	all_type_result_cluster2=all_type_result_cluster.copy()
	all_type_result_cluster2_region=all_type_result_cluster2[(all_type_result_cluster2["breakpoint_human"]==0)]
	all_type_result_cluster2=all_type_result_cluster2[(all_type_result_cluster2["breakpoint_human"]>0)]
	all_type_result_cluster2_region=all_type_result_cluster2_region[["chr_human_chimeric","start_human_chimeric","start_human_chimeric"]]
	all_type_result_cluster2_region["variant"]=[0]*len(all_type_result_cluster2_region.index)
	all_type_result_cluster2_region["ref"]=[0]*len(all_type_result_cluster2_region.index)
	all_type_result_cluster2=all_type_result_cluster2[["chr_human","breakpoint_human","breakpoint_human"]]
	all_type_result_cluster2['variant']=[0]*len(all_type_result_cluster2.index)
	all_type_result_cluster2['ref']=[0]*len(all_type_result_cluster2.index)
	all_type_result_cluster2=all_type_result_cluster2.drop_duplicates(subset=["chr_human","breakpoint_human"],keep='first')
	all_type_result_cluster2_region.columns=all_type_result_cluster2.columns
	all_type_result_cluster_anno=pd.concat([all_type_result_cluster2,all_type_result_cluster2_region])
	all_type_result_cluster_anno["variant"]=all_type_result_cluster_anno["variant"].astype("int")
	all_type_result_cluster_anno["ref"]=all_type_result_cluster_anno["ref"].astype("int")
	all_type_result_cluster_anno.to_csv(outdir+"/"+"human_annovar.avinput",index=False,sep='\t',header=False)














