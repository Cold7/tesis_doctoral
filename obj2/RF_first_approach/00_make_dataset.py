#importing libraries
# bed format https://genome.ucsc.edu/FAQ/FAQformat.html#format1

import argparse   # arguments parser
from glob import glob
import sys
import os
from Bio import SeqIO
import tempfile
import random

# importing my own libraries
sys.path.append(str(os.path.realpath(__file__))[:-8]+"/libs")
from getFeatures import getFeatures
from tfbs import ensembl


from save import saveDataset

if __name__=="__main__":

	#input
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genome", help="Path to genome file (in fasta format)", required = True)
	parser.add_argument("-f","--fragment_size", help="number of bases to fragment DNA sequence. Default: 500", type = float, default=500)
	parser.add_argument("-p","--percentage", help="For some experiments you may need to merge them, so to do it you need to incluye a percentage of files that contain the feature. Default: 100", type = float, default = 100)
	parser.add_argument("-o","--output", help="folder where results will be placed. Default: ./", default = "./")
	parser.add_argument("-t","--tmp", help="folder for temporal files. Default: ./temp/", default = "./temp/")
	parser.add_argument("-dttu","--dataTypeToUse", help="Definition to use: \"bedscore\" or \"presence\". Default: bedscore", default = "bedscore")
	parser.add_argument("-fm","--filling_mode", help="If you wish your data in  binary (is or not present) or in number (average characteristic in the current fragment or percentage). Options are \"binary\", \"percentage\" or \"average\". Default: binary", default="binary")
	parser.add_argument("-n","--nproc", help="Number of processors to use to compute TFBS active/inactive. Default: 1", type = int, default=1)

	####################################
	##
	## DNase-seq arguments
	##
	####################################
	parser.add_argument("-ds","--DNase",help="Option to generate features for DNAse.", action="store_true")
	parser.add_argument("-dsf","--DNase_file", help="DNase files.")

	####################################
	##
	## DNA methylation arguments
	##
	####################################
	parser.add_argument("-dm","--dna_methylation",help="Option to generate features for dna methylation.", action="store_true")
	parser.add_argument("-mf","--methylation_File", help="Folder where cpg DNAme files (in bed format) are located.")
	
	####################################
	##
	## ctcf arguments
	##
	####################################
	parser.add_argument("-c","--ctcf",help="Option to generate features for ctcf.", action="store_true")
	parser.add_argument("-cFile","--ctcf_File", help="Folder where ctcf results (in bed file) are located.")

	####################################
	##
	## yy1 arguments
	##
	####################################
	parser.add_argument("-y","--yy1",help="Option to generate features for yy1.", action="store_true")
	parser.add_argument("-yFile","--yy1_File", help="Folder where yy1 results (in bed file) are located.")

	####################################
	##
	## PolR2a arguments
	##
	####################################
	parser.add_argument("-pol","--polr2a",help="Option to generate features for polr2a.", action="store_true")
	parser.add_argument("-pFile","--polr2a_File", help="Folder where polr2a results (in bed file) are located.")

    ####################################
    ##
    ## EP300 arguments
    ##
    ####################################
	parser.add_argument("-ep","--ep300",help="Option to generate features for ep300.", action="store_true")
	parser.add_argument("-epFile","--ep300_File", help="Folder where polr2a results (in bed file) are located.")

	####################################
	##
	## Histone marks arguments
	##
	####################################
	parser.add_argument("-hm","--histonemarks",help="Option to generate features for histone marks.", action="store_true")
	parser.add_argument("-hmf","--histoneMarksFiles", nargs="*", help="Folder where histone marks ChIP-seq results (in bed file) are located.")

	####################################
	##
	## TFBSac arguments for ensembl
	##
	####################################
	parser.add_argument("-tf","--tf_folder", help="text file with a list of TF name and bed file path.", required = True)
	parser.add_argument("-ens","--ensembl_gtf", help="Ensembl gtf with tfbs sites.", required = True)
#	parser.add_argument("-bin","--binarize",help="Option to keep those classes 1 or -1.", action="store_true")

	#parsing arguments
	args = parser.parse_args()

	tempFilenames = [] # to save temporal data
	abc = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
	ID ="" # a random identifier to save data
	for i in range (10):
		ID += random.choice(abc)

	if args.fragment_size <= 0:
		print ("The number of kb to parse dna (--DNA_split option) must be higher than 0. Exiting")
		exit()
	#for each chr
	for record in SeqIO.parse(args.genome, "fasta"):
		currentCHR = record.id
		#currentCHR = "chr2" #delete or comment this line  for full dataset creation toghether with the exit line at the end of the script

		###################################
		##
		## getting vectors
		##
		###################################

		#for dnase
		DNase = None
		if args.DNase:
			print("DNAse on chr", currentCHR)
			DNase = getFeatures(args.DNase_file, args.filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR, args.dataTypeToUse)
			tempFilenames.append(args.tmp+"/DNase_"+ID)
			f = open(args.tmp+"/DNase_"+ID,"w")
			f.write("DNase\n")
			for i in DNase:
				f.write(str(i)+"\n")
			f.close()
			del(DNase)

		methylation = None
		#for dna methylation
		if args.dna_methylation:
			print("methylation on chr", currentCHR)
			#methylation = getFeatures(args.methylation_Folder, args.DNA_meth_filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR)
			methylation = getFeatures(args.methylation_File, "percentage", args.fragment_size, args.percentage, args.genome, currentCHR, "presence")
			tempFilenames.append(args.tmp+"/methylation_"+ID)
			f = open(args.tmp+"/methylation_"+ID,"w")
			f.write("methylation\n")
			for i in methylation:
				f.write(str(i)+"\n")
			f.close()
			del(methylation)

		#for ctcf
		if args.ctcf:
			print("ctcf for chr", currentCHR)
			ctcf = getFeatures(args.ctcf_File, args.filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR, args.dataTypeToUse)
			tempFilenames.append(args.tmp+"/ctcf_"+ID)
			f = open(args.tmp+"/ctcf_"+ID,"w")
			f.write("ctcf\n")
			for i in ctcf:
				f.write(str(i)+"\n")
			f.close()
			del(ctcf)

		#for yy1
		if args.yy1:
			print("yy1 for chr", currentCHR)
			yy1 = getFeatures(args.yy1_File, args.filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR, args.dataTypeToUse)
			tempFilenames.append(args.tmp+"/yy1_"+ID)
			f = open(args.tmp+"/yy1_"+ID,"w")
			f.write("yy1\n")
			for i in yy1:
				f.write(str(i)+"\n")
			f.close()
			del(yy1)

		#for PolR2a
		if args.polr2a:
			print("polr2a for chr", currentCHR)
			polr2a = getFeatures(args.polr2a_File, args.filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR, args.dataTypeToUse)
			tempFilenames.append(args.tmp+"/polr2a_"+ID)
			f = open(args.tmp+"/polr2a_"+ID,"w")
			f.write("polr2a\n")
			for i in polr2a:
				f.write(str(i)+"\n")
			f.close()
			del(polr2a)	

		#for EP300
		if args.ep300:
			print("ep300 for chr", currentCHR)
			ep300 = getFeatures(args.ep300_File, args.filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR, args.dataTypeToUse)
			tempFilenames.append(args.tmp+"/ep300_"+ID)
			f = open(args.tmp+"/ep300_"+ID,"w")
			f.write("ep300\n")
			for i in ep300:
				f.write(str(i)+"\n")
			f.close()
			del(ep300)

		#for histone marks
		if args.histonemarks:
			print("histones on chr", currentCHR)
			for HMFile in args.histoneMarksFiles:
				print("\tusing", HMFile)
				histone_mark = getFeatures(HMFile, args.filling_mode, args.fragment_size, args.percentage, args.genome, currentCHR, args.dataTypeToUse)
				HMName = HMFile.split("/")[-1].split("_")[0]
				tempFilenames.append(args.tmp+"/"+HMName+"_"+ID)
				f = open(args.tmp+"/"+HMName+"_"+ID,"w")
				f.write(HMName+"\n")
				for i in histone_mark:
					f.write(str(i)+"\n")
				f.close()
				del(histone_mark)

#		#for active TFBS using ensembl data
#		if args.tf_folder:
#			print("TFBSac using ensembl on chr", currentCHR)
#			TFBSac_ensembl = ensembl(args.genome, currentCHR, args.fragment_size, args.percentage, args.filling_mode, args.tf_folder, args.tmp, args.nproc)
#			tempFilenames.append(args.tmp+"/TFBSac_ensembl_"+ID)
#			f = open(args.tmp+"/TFBSac_ensembl_"+ID,"w")
#			f.write("TFBSac_ensembl\n")
#			for i in TFBSac_ensembl:
#				f.write(str(i)+"\n")
#			f.close()
#			del(TFBSac_ensembl)

		if args.tf_folder and args.ensembl_gtf != None:
			print("TFBSac using ensembl on chr", currentCHR)
			TFBSac_ensembl = ensembl(args.genome, currentCHR, args.fragment_size, args.percentage, args.filling_mode, args.tf_folder, args.tmp, args.nproc, args.ensembl_gtf)
			tempFilenames.append(args.tmp+"/TFBSac_ensembl_"+ID)
			f = open(args.tmp+"/TFBSac_ensembl_"+ID,"w")
			f.write("TFBSac_ensembl\n")
			for i in TFBSac_ensembl:
				f.write(str(i)+"\n")
			f.close()
			del(TFBSac_ensembl)


		print("saving chr", currentCHR)
		saveDataset(tempFilenames, args.output, currentCHR)
		print("deleting temporal files for chr", currentCHR)
		for file in tempFilenames:
			os.system("rm "+file)


#python3 main.py  -g ../k562_data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fa -dttu bedscore -fm average -ds -dsf ../k562_data/dnase/dnase_enrichment_hotspot_ENCFF828WSI.bed -dm -dname ../k562_data/wgbs/salida.bed -c -cFile ../k562_data/ctcf/CTCF_optimal_ENCFF651WYF.bed -y -yFile ../k562_data/yy1/ENCFF629EQZ_optimal.bed -pol -pFile ../k562_data/polr2a/POLR2A_ENCFF130WVT.bed -hm -hmf ../k562_data/HM/H3K27ac_ENCFF045OHM.bed ../k562_data/HM/H3K36me3_ENCFF753VAK.bed ../k562_data/HM/H3K4me3_ENCFF465RJJ.bed ../k562_data/HM/H3K27me3_ENCFF233ODK.bed ../k562_data/HM/H3K4me1_ENCFF359WWB.bed ../k562_data/HM/H3K9me3_ENCFF361WTS.bed -tf ../k562_data/ensembl_chipseq/ -ens ../k562_data/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff





