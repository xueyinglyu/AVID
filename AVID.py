import os,sys
import argparse
from optparse import OptionParser
import pysam
import numpy as np
import statistics
import math
import subprocess
import importlib.util




################################################isolate virus related fastq file ########################################################
### remember to add code for config.txt in workdir directory

def main ():

	usage="python AVID.py -1 test_R1.fq.gz -2 test_R2.fq.gz -d /output/dir/ -s test -r virus_bwa_index.fa -g hg19_blastn_index -p virus_blastn_index -H hg19_bwa_index.fa -l 10 -q 10 -t 3 -@ 8 -v hg19 -a virus_annotation.bed -R read_len -I insert_size -F True -S True"
	parser=OptionParser(usage=usage)
	parser.add_option("-1","--fastq1",dest="fastq1",help="fastq file read1")
	parser.add_option("-2","--fastq2",dest="fastq2",help="fastq file read2")
	parser.add_option("-d","--directory",dest="directory",help="output folder")
	parser.add_option("-s","--sample",dest="sample",help="sample id")
	parser.add_option("-r","--reference",dest="reference",help="virus reference genome(s) fasta file for BWA")
	parser.add_option("-p","--pathogen",dest="pathogen",help="the virus reference genome for BlastN")
	parser.add_option("-H","--Human",dest="human_index",help="human genome bwa index")
	parser.add_option("-l","--len",dest="minlen",default="10",help="minimum length for softclip reads")
	parser.add_option("-q","--quality",dest="quality",help="minimum quality for BWA mapping")
	parser.add_option("-t","--threshold",dest="threshold",default=3,help="minimum reads (softclip and chimeric reads) supported the breakpoint")
	parser.add_option("-@","--threads",dest="threads",help="threads for BWA")
	parser.add_option("-v","--version",dest="human_version",default="hg38",help="human genome version. eg hg19 hg38")
	parser.add_option("-a","--annotation",dest="virus_annotation",help="virus annotation file, chr start end region strand gene")
	parser.add_option("-R","--Readlen",dest="read_len",help="read length")
	parser.add_option("-I","--insert_size",dest="insert",help="insert size")
	parser.add_option("-T","--data_type",dest="data_type",default="DNA",help="DNA or RNA")
	parser.add_option("-F","--chr-incompatibility",dest="chr_incompatibility",default="True",help="Activate the filtering of chromosomal incompatibility (False or True)")
	parser.add_option("-S","--filter-insertsize",dest="filter_insert",default="True",help="Activate the filtering of DNA fragment insert size (False or True)")
	(options,args)=parser.parse_args()
	fastq1=options.fastq1
	fastq2=options.fastq2
	sample=options.sample
	directory=options.directory
	virus_bwa=options.reference
	pathogen_blastn=options.pathogen
	host_bwa=options.human_index
	minlen=options.minlen
	quality=options.quality
	threshold=options.threshold
	threads=options.threads
	human_version=options.human_version
	virus_annotation=options.virus_annotation
	read_len=int(float(options.read_len))
	insert_size=int(options.insert)
	data_type=options.data_type
	chr_incompatibility=options.chr_incompatibility
	filter_insert=options.filter_insert


	Virusclip_PATH = os.path.dirname(os.path.realpath(__file__))
	print(Virusclip_PATH)

	

	def cal_virus_len(virus_fasta):
		virus_file=open(virus_bwa,'r')
		ref_seq=''
		for line in virus_file:
			if ">" in line:
				continue
			ref_seq=ref_seq+line.strip()
		virus_len=len(ref_seq)
		return virus_len

	def BWA_fasta(fasta,reference,outfile,sample,threads):
		mapping_cmd="%s/bwa mem -t %s -R \"@RG\\tID:%s\\tSM:%s\" %s %s -o %s"%(Virusclip_PATH,threads,sample,sample,reference,fasta,outfile)
		os.system(mapping_cmd)

	def STAR_fasta(fasta,outfile,sample,threads,outdir):
		mapping_cmd="%s/STAR --runThreadN %s --genomeDir /pathowh01/disk1/lyuxy/HBV/reference/Human/STAR --outFileNamePrefix %s/%s --readFilesIn %s --outSAMtype SAM --outSAMmode SAM"%(Virusclip_PATH,threads,outdir,sample,fasta)
		os.system(mapping_cmd)
		os.system("mv %s/%sAligned.out.sam %s"%(outdir,sample,outfile))

	def BLASTN(fasta,reference,outfile,threads):
		#print("blastn -db %s -query %s -out %s -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand\" -max_target_seqs 5 -word_size 11 -best_hit_overhang 0.1 -num_threads %s"%(reference,fasta,outfile,threads))
		os.system("%s/blastn -db %s -query %s -out %s -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand\" -max_target_seqs 5 -word_size 11 -best_hit_overhang 0.1 -num_threads %s"%(Virusclip_PATH,reference,fasta,outfile,threads))

	def InsertSize(samfile): # type3 paired end mapped bam file
		with pysam.AlignmentFile(samfile,"r") as fh:
				insize_list=[]
				for record in fh:
					if(record.template_length<0):
						continue
					insize_list.append(record.template_length)
				average_insize=round(np.mean(insize_list),2)
				sd_insize=statistics.stdev(insize_list)
				insert_size=average_insize+3*sd_insize
				insert_size=round(insert_size,2)
				return insert_size

	################################################################################## check tools ###########################################################################################################

	def is_package_installed(package_name):
		return importlib.util.find_spec(package_name) is not None

	#if not os.path.exists("%s/Rscript"%(Virusclip_PATH)):
	if not os.path.exists("%s/Rscript"%(Virusclip_PATH)):
		print("Error: Rscript has not intalled!")
		sys.exit()

	if not is_package_installed('numpy'):
		print("Error: python package numpy is not installed.")
		sys.exit()

	if not is_package_installed('pysam'):
		print("Error: python package pysam is not installed.")
		sys.exit()

	if not is_package_installed('statistics'):
		print("Error: python package statistics is not installed.")
		sys.exit()

	if not is_package_installed('math'):
		print("Error: python package math is not installed.")
		sys.exit()

	if not is_package_installed('subprocess'):
		print("Error: python package subprocess is not installed.")
		sys.exit()

	if not is_package_installed('importlib'):
		print("Error: python package importlib is not installed.")
		sys.exit()

	if not is_package_installed('pandas'):
		print("Error: python package pandas is not installed.")
		sys.exit()

	if not is_package_installed('optparse'):
		print("Error: python package optparse is not installed.")
		sys.exit()

	if not is_package_installed('argparse'):
		print("Error: python package argparse is not installed.")
		sys.exit()

	if not is_package_installed('re'):
		print("Error: python package re is not installed.")
		sys.exit()
		
	print("###############################################################\n")
	print("All tools are installed !\n")
	print("###############################################################\n")
	print("STEP1: Starting to extract virus related reads......\n")
	print("###############################################################")



	change_chmod="chmod 777 %s/*"%(Virusclip_PATH)
	os.system(change_chmod)

	fh=open("%s/config.txt"%(directory),'w')
	fh.write("threads %s\nmin_sc_size 20\nmax_sc_dist 10\nread_len %s"%(int(threads),read_len))
	fh.close()
	os.system("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s"%(Virusclip_PATH))
	isolate_virus_cmd="%s/isolate_relevant_pairs_fq %s %s %s %s %s %s" %(Virusclip_PATH,fastq1,fastq2,host_bwa,virus_bwa,directory,directory)
	os.system(isolate_virus_cmd)

	if os.path.exists("%s/retained-pairs_1.fq"%(directory)):
		print("\nSTEP1: Done!\n")
	else:
		print("STEP1: ERROR!")
		sys.exit()

	################################################# mapping reads to virus genome ###########################################################3
	print("###############################################################\n")
	print("STEP2: Starting to align to virus genome using BWA-0.7.17......\n")
	print("###############################################################")
	virus_bwa_cmd="%s/bwa mem -t %s -R \"@RG\\tID:%s\\tSM:%s\" %s %s/retained-pairs_1.fq %s/retained-pairs_2.fq -o %s/%s.virus.sam"%(Virusclip_PATH,threads,sample,sample,virus_bwa,directory,directory,directory,sample)
	os.system(virus_bwa_cmd)


	if os.path.exists("%s/%s.virus.sam"%(directory,sample)):
		print("\nSTEP2: Done!\n")
	else:
		print("STEP2: ERROR!")
		sys.exit()



	######################################################### subtype reads #############################################################################
	print("###############################################################\n")
	print("STEP3: Starting to subtype reads......\n")
	print("###############################################################")
	os.system("python %s/Summary_BWA.py -i %s/%s.virus.sam -o %s/%s.bwa.virus.out -q %s"%(Virusclip_PATH,directory,sample,directory,sample,quality))
	os.system("python %s/Subtype.py -i %s/%s.bwa.virus.out -d %s -s HBV -l %s"%(Virusclip_PATH,directory,sample,directory,minlen))
	os.system("python %s/Extract_mapping_seq.py -i %s/%s.bwa.virus.out -o %s/%s.virus.mapping.fasta"%(Virusclip_PATH,directory,sample,directory,sample))


	######################################################### mapping to human genome#######################################################################
	if(data_type=="DNA"):
		BWA_fasta("%s/HBV_Unmapped.fasta"%(directory),host_bwa,"%s/HBV_Unmapped_human.sam"%(directory),sample,threads)
		BWA_fasta("%s/HBV_Partial_BWA.fasta"%(directory),host_bwa,"%s/HBV_Partial_human.sam"%(directory),sample,threads)
		BWA_fasta("%s/%s.virus.mapping.fasta"%(directory,sample),host_bwa,"%s/Human_check.sam"%(directory),sample,threads)


	if(data_type=="RNA"):
		STAR_fasta("%s/HBV_Unmapped.fasta"%(directory),"%s/HBV_Unmapped_human.sam"%(directory),sample,threads,directory)
		STAR_fasta("%s/HBV_Partial_BWA.fasta"%(directory),"%s/HBV_Partial_human.sam"%(directory),sample,threads,directory)
		STAR_fasta("%s/%s.virus.mapping.fasta"%(directory,sample),"%s/Human_check.sam"%(directory),sample,threads,directory)


	####################################################### summary human alignment ################################################

	os.system("python %s/Summary_BWA.py -i %s/HBV_Unmapped_human.sam -o %s/HBV_Unmapped_human.out -q %s"%(Virusclip_PATH,directory,directory,quality))
	os.system("python %s/Summary_BWA.py -i %s/HBV_Partial_human.sam -o %s/HBV_Partial_human.out -q %s"%(Virusclip_PATH,directory,directory,quality))
	os.system("python %s/Subtype.py -i %s/HBV_Unmapped_human.out -d %s -s Unmap_HBV -l %s"%(Virusclip_PATH,directory,directory,minlen))
	os.system("python %s/Subtype.py -i %s/HBV_Partial_human.out -d %s -s Partial_HBV -l %s"%(Virusclip_PATH,directory,directory,minlen))

	# ######################################################## check virus and human alignment ############################################

	os.system("python %s/Summary_BWA.py -i %s/Human_check.sam -o %s/Human_check.out -q %s"%(Virusclip_PATH,directory,directory,quality))
	subprocess.call(['bash','-c',"grep -v -F %s/HBV_Partial.out -f <(grep -w 'F' %s/Human_check.out | awk '{print $1}') > %s/HBV_Partial.out2"%(directory,directory,directory)])
	subprocess.call(['bash','-c',"grep -v -F %s/HBV_Full.out -f <(grep -w 'F' %s/Human_check.out | awk '{print $1}') > %s/HBV_Full.out2"%(directory,directory,directory)])
	subprocess.call(['bash','-c',"cat %s/HBV_Unmapped_human.out <(grep -w 'F' %s/Human_check.out) > %s/HBV_Unmapped_human.out2"%(directory,directory,directory)])
	os.system("mv %s/HBV_Partial.out2 %s/HBV_Partial.out"%(directory,directory))
	os.system("mv %s/HBV_Full.out2 %s/HBV_Full.out"%(directory,directory))
	os.system("mv %s/HBV_Unmapped_human.out2 %s/HBV_Unmapped_human.out"%(directory,directory))

	if(os.path.exists("%s/Unmap_HBV_Partial.fasta"%(directory)) and (os.path.exists("%s/HBV_Full.out"%(directory))) and (os.path.exists("%s/HBV_Partial.out"%(directory))) and (os.path.exists("%s/HBV_Unmapped.out"%(directory))) and (os.path.exists("%s/Partial_HBV_Full.out"%(directory))) and (os.path.exists("%s/Partial_HBV_Partial.out"%(directory))) and (os.path.exists("%s/Partial_HBV_Unmapped.out"%(directory))) and (os.path.exists("%s/Unmap_HBV_Full.out"%(directory))) and (os.path.exists("%s/Unmap_HBV_Partial.out"%(directory))) and os.path.exists("%s/Unmap_HBV_Unmapped.out"%(directory))):
		print("\nSTEP3: Done!\n")
	else:
		print("STEP3: ERROR!")
		sys.exit()

	################################################# mapping reads to virus genome ###########################################################3
	print("###############################################################\n")
	print("STEP4: Starting to identify viral integration breakpoints......\n")
	print("###############################################################")
	virus_len=cal_virus_len(virus_bwa)
	BLASTN("%s/Unmap_HBV_Partial.fasta"%(directory),pathogen_blastn,"%s/Unmap_HBV_Partial_HBV_blastn.out"%(directory),threads)
	os.system("python %s/Find_breakpoint3.py -d %s -s %s -t %s -l %s -v %s -S %s -F %s"%(Virusclip_PATH,directory,insert_size,threshold,read_len,str(virus_len),filter_insert,chr_incompatibility))
	
	if(os.path.exists("%s/human_annovar.avinput"%(directory)) and (os.path.exists("%s/all_type_result.out"%(directory)))):
		print("\nSTEP4: Done!\n")
	else:
		print("STEP4: ERROR!")
		sys.exit()
	################################################## annotation #################################################################################
	print("###############################################################\n")
	print("STEP5: Starting to do annotation......\n")
	print("###############################################################")





	if os.stat("%s/human_annovar.avinput"%(directory)).st_size>0:
		print("perl %s/annotate_variation.pl -build %s -out %s/%s %s/human_annovar.avinput %s/annotation/annovar/%s"%(Virusclip_PATH,human_version,directory,sample,directory,Virusclip_PATH,human_version))
		os.system("perl %s/annotate_variation.pl -build %s -out %s/%s %s/human_annovar.avinput %s/annotation/annovar/%s"%(Virusclip_PATH,human_version,directory,sample,directory,Virusclip_PATH,human_version))
		os.system("python %s/merge_annotation.py -H %s/%s.variant_function -v %s -i %s/all_type_result.out -o %s/%s.AVID.txt"%(Virusclip_PATH,directory,sample,virus_annotation,directory,directory,sample))
		os.system("python %s/merge_repeat.py -i %s/%s.AVID.txt -o %s/%s.AVID.final.txt -t %s"%(Virusclip_PATH,directory,sample,directory,sample,threshold))
########################## install samtools	########################## 	
		os.system("ls %s/*.virus.sam | while read f;do(outfile=${f/sam/insert}; %s/samtools stat $f | grep ^IS | cut -f 2- > $outfile);done"%(directory,Virusclip_PATH))
		
		if(os.path.exists("%s/%s.AVID.final.txt"%(directory,sample))):
			print("\nSTEP5: Done!\n")
		else:
			print("STEP5: ERROR!")
			sys.exit()

		print("###############################################################\n")
		print("STEP6: Starting to visualize result......\n")
		print("###############################################################")


		os.system("Rscript %s/visualization.R %s %s %s"%(Virusclip_PATH,directory,virus_len,human_version))
	else:
		os.system("touch %s/%s.AVID.txt"%(directory,sample))
		os.system("touch %s/%s.AVID.final.txt"%(directory,sample))

	os.system("rm %s/*.out"%(directory))
	os.system("rm %s/*.sam"%(directory))
	os.system("rm %s/*.fasta"%(directory))
	os.system("rm %s/type*"%(directory))
	os.system("rm %s/*.AVID.txt"%(directory))
	os.system("rm %s/*.virus.insert"%(directory))

	if(os.path.exists("%s/plot/circos_plot.pdf"%(directory))):
		print("\nSTEP6: Done!\n")
	else:
		print("STEP6: ERROR!")
		sys.exit()

	print("###############################################################\n")
	print("\n\nAVID COMPLETED!!!\n\nPlease check the result in %s\n"%(directory))



if __name__ == '__main__':
	main()






























