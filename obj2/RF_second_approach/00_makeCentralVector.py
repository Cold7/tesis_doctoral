import argparse   # arguments parser
from Bio import SeqIO
import pandas as pd
import numpy as np
import multiprocessing as mp
import os

def readFourCols(currentCHR, file, length, coordBased):
	vector = np.zeros(int(length/args.fragment)+1, dtype=np.int32)
	countVector = np.zeros(int(length/args.fragment)+1, dtype=np.int32) #to compute an average if two or more times data was append in the same position
	f = open(file,"r")
	for line in f:
		aux = line[:-1].split("\t")
		if currentCHR == aux[0]:
			init = int(aux[1])
			end = int(aux[2]) - 1
			score = float(aux[4])
			if coordBased == 1:
				init -= 1
			init = int(init/args.fragment)
			end = int(end/args.fragment) 
			while init <= end:
				if args.filling_mode == "score":
					vector[init] += score
					countVector[init] += 1
				else:
					vector[init] += 1
					countVector[init] += 1
				init += 1

	for i in range(len(vector)):
		if vector[i] != 0:
			vector[i] /= countVector[i]
	f.close()
	return vector


def doData(data):
	currentCHR, length = data
	dataset = pd.DataFrame()
	
	##########
	## dnase
	#########
	if args.DNase:
#		print("\tLoading DNAse for",currentCHR)
		dataset["DNAse-seq"] = readFourCols(currentCHR, args.DNase_file, length, args.DNase_base)

	####################
	## dna methylation
	###################
	if args.dna_methylation:
#		print("\tLoading methylation for",currentCHR)
		dataset["methylation"] = readFourCols(currentCHR, args.dna_methylation_File, length, args.dna_methylation_base)

	##########
	## %CG
	#########
	if args.dna_CG:
		#looking for citocines in the genome
		auxVector = np.zeros(int(length/args.fragment)+1, dtype=np.int32)
		for record in SeqIO.parse(genomeFile, "fasta"):
			if record.id == currentCHR:
				current_position = 0
				while current_position < length:
					substring = record.seq[current_position:current_position+args.fragment].lower()
					auxVector[int(current_position/args.fragment)] = substring.count("c")+substring.count("g")
					current_position += args.fragment
		dataset["CG"] = auxVector

	###########
	## ctcf
	##########
	if args.ctcf:
#		print("\tLoading CTCF for",currentCHR)
		dataset["CTCF"] = readFourCols(currentCHR, args.ctcf_File, length, args.ctcf_base)
	
	###############
	## for yy1
	##############
	if args.yy1:
#		print("\tLoading yy1 for",currentCHR)
		dataset["YY1"] = readFourCols(currentCHR, args.yy1_File, length, args.yy1_base)

	################
	## PolR2a
	###############
	if args.polr2a:
#		print("\tLoading polr2a for",currentCHR)
		dataset["POLR2A"] = readFourCols(currentCHR, args.polr2a_File, length, args.polr2a_base)

	#############
	## EP300
	############
	if args.ep300:
#		print("\tLoading ep300 for",currentCHR)
		dataset["ep300"] = readFourCols(currentCHR, args.ep300_File, length, args.ep300_base)

	###############
	## cbp
	##############
	if args.crebbp:
#		print("\tLoading cbp for",currentCHR)
		dataset["cbp"] = readFourCols(currentCHR, args.crebbp_File, length, args.crebbp_base)

	#################
	## for rnaseq
	################
	if args.rna:
#		print("\tLoading rna for",currentCHR)
		dataset["rnaseq"] = readFourCols(currentCHR, args.rna_File, length, "1")  #1 is to use 1-base

	################
	## TSS & Gene
	##############
	if args.gene:
		vectorGene = np.zeros(int(length/args.fragment)+1, dtype=np.int32)
		vectorTSS = np.zeros(int(length/args.fragment)+1, dtype=np.int32)
		gtf = open(args.gene_GTF,"r")
		for line in gtf:
			aux = line[:-1].split("\t")
			if aux[0] == currentCHR and aux[2] == "gene":
				init = int(aux[3])-1 #assuming 1-base
				end = int(aux[4])-1 #same thing
				init = int(init/args.fragment)
				end = int(end/args.fragment)
				vectorTSS[init] = 1
				while init <= end:
					vectorGene[init] += 1
					init += 1
		dataset["TSS"] = vectorTSS
		dataset["Gene"] = vectorGene

	#######################
	## histone marks
	######################
	HMs = {}
	if args.histonemarks:
		print("\tLoading histones for",currentCHR)
		for HMFile in args.histoneMarks_Files:
			name = HMFile.split("/")[-1]
#			print("\t\tHM: ", name)
			dataset[name] =  readFourCols(currentCHR, HMFile, length, args.histoneMarks_base)


	######################
	## current Class
	#####################

#	print("\tLoading "+args.Class+" active for",currentCHR)

	##doing a list of TF for I have data, then I will add motifs just for these chipseq
	f = open(args.class_File,"r")
	tfList = []
	for line in f:
		tfList.append(line.split("\t")[0].upper())
	f.close()


	##class definition##
	f = open(args.class_definition,"r")
	dictSite = {}
	count = 0
	for line in f:
		aux = line[:-1].split("\t")
		if aux[1] == currentCHR:
			motifs = aux[4].split(",")
			motifs = [x.upper() for x in motifs] # motifs in uppper
			dictSite[count] = {"init":int(aux[2])-1,"end":int(aux[3])-1, "status":False, "motifs":[]} #-1 assuming 1-base

			for motif in motifs:
				if motif in tfList:
					dictSite[count]["motifs"].append(motif)
			count += 1
	f.close()

	##chip overlapping class##
	auxVector = np.zeros(int(length/args.fragment)+1, dtype=np.int8)
	class_File = open(args.class_File,"r")
	used = [] #to know used motifs
	for line in class_File:
		data = line[:-1].split("\t")
		name= data[0].upper()
		path = data[1]
		f = open(path, "r")
		for line in f:
			aux = line.split("\t")
			if aux[0] == currentCHR:
				class_init = int(aux[1])
				class_end = int(aux[2]) - 1
				if args.class_base == "1":
					class_init -= 1
				for i in range(len(dictSite)):
					if name in dictSite[i]["motifs"]:
						if class_end > dictSite[i]["init"] and class_init < dictSite[i]["end"]:
							dictSite[i]["status"] = True
							if name not in used:
								used.append(name)
		f.close()
	class_File.close()

	##looking for no used tfs
	noUsed = []
	for tf in tfList:
		if tf not in used:
			noUsed.append(tf)
	##filling vector with active or inactive class##
	for site in dictSite:
		init = int(dictSite[site]["init"]/args.fragment)
		end = int(dictSite[site]["end"]/args.fragment)
		toFill = -1
		if dictSite[site]["status"] == True:
			toFill = 1

		while init <= end:
			auxVector[init] = toFill
			init += 1
	dataset[args.Class] = auxVector

	#########################
	## saving dataset
	########################
#	print("\tsaving dataset for",currentCHR)
	dataset.to_csv(args.output+"/centralVector_"+currentCHR+".tsv", sep="\t", index = False)

	print("\tDone for chr "+currentCHR+" to be saved on "+args.output+" with fragment length "+str(args.fragment)+". Used TFs are: "+str(used)+". and no used TFs: "+str(noUsed))

if __name__ == "__main__":
	global args, genomeFile

	#input
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genome", help="Path to genome file (in fasta format)", required = True)
	parser.add_argument("-f","--fragment", help="number of bases to fragment DNA sequence. Default: 50", type = int, default=50)
	parser.add_argument("-o","--output", help="folder where results will be placed. Default: ./", default = "./")
	parser.add_argument("-fm","--filling_mode", help="form as data is presented. Options are score to use the score in tsv file or percentage to use percentage of occupancy. Default: score", default="score",choices=["score","percentage"])
	parser.add_argument("-n","--nproc", help="Number of processors to use to compute TFBS active/inactive. Default: 1", type = int, default=1)

	####################################
	##
	## DNase-seq arguments
	##
	####################################
	parser.add_argument("-ds","--DNase",help="Option to generate features for DNAse.", action="store_true")
	parser.add_argument("-dsf","--DNase_file", help="DNase files.")
	parser.add_argument("-dsb","--DNase_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## DNA methylation arguments
	##
	####################################
	parser.add_argument("-dm","--dna_methylation",help="Option to generate features for dna methylation.", action="store_true")
	parser.add_argument("-dmf","--dna_methylation_File", help="Folder where cpg DNAme files (in bed format) are located.")
	parser.add_argument("-dmb","--dna_methylation_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")
	parser.add_argument("-dc","--dna_CG",help="Option to generate percentage of C-G or G-C in the fragment.", action="store_true")

	####################################
	##
	## ctcf arguments
	##
	####################################
	parser.add_argument("-c","--ctcf",help="Option to generate features for ctcf.", action="store_true")
	parser.add_argument("-cFile","--ctcf_File", help="Folder where ctcf results (in bed file) are located.")
	parser.add_argument("-cb","--ctcf_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## yy1 arguments
	##
	####################################
	parser.add_argument("-y","--yy1",help="Option to generate features for yy1.", action="store_true")
	parser.add_argument("-yFile","--yy1_File", help="Folder where yy1 results (in bed file) are located.")
	parser.add_argument("-yb","--yy1_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## PolR2a arguments
	##
	####################################
	parser.add_argument("-pol","--polr2a",help="Option to generate features for polr2a.", action="store_true")
	parser.add_argument("-pFile","--polr2a_File", help="Folder where polr2a results (in bed file) are located.")
	parser.add_argument("-pb","--polr2a_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## EP300 arguments
	##
	####################################
	parser.add_argument("-ep","--ep300",help="Option to generate features for ep300.", action="store_true")
	parser.add_argument("-epFile","--ep300_File", help="Folder where ep300 results (in bed file) are located.")
	parser.add_argument("-eb","--ep300_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## cbp arguments
	##
	####################################
	parser.add_argument("-cbp","--crebbp",help="Option to generate features for cbp.", action="store_true")
	parser.add_argument("-cbpFile","--crebbp_File", help="Folder where cbp results (in bed file) are located.")	
	parser.add_argument("-cbpb","--crebbp_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## Histone marks arguments
	##
	####################################
	parser.add_argument("-hm","--histonemarks",help="Option to generate features for histone marks.", action="store_true")
	parser.add_argument("-hmf","--histoneMarks_Files", nargs="*", help="path for each histone mark (with title as name of the current HM).")
	parser.add_argument("-hmb","--histoneMarks_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="0")

	####################################
	##
	## Gene & RNA data
	##
	####################################
	parser.add_argument("-r","--rna",help="Option to generate features for RNA count.", action="store_true")
	parser.add_argument("-rf","--rna_File",help="Option to generate features for RNA. Input required is a tsv or bed file with 5 columns (chr init end strand and count)")
	parser.add_argument("-gen","--gene",help="Option to generate features for genes and TSS.", action="store_true")
	parser.add_argument("-gg","--gene_GTF", help="path for gene GTF.")

	####################################
	##
	## Class arguments
	##
	####################################
	parser.add_argument("-C","--Class", help="class to use. It will be searched in file defined in --class_definition option.", required = True)
	parser.add_argument("-cf","--class_File", help="file with path to 5 columns chip-seq results for the current class.", required = True)
	parser.add_argument("-cd","--class_definition", help="tsv file with site definition (format is id, chr, init, end and  motifs).", required = True)
	parser.add_argument("-Cb","--class_base", help="to use 0-based or 1-based coord system. Default: 0", choices=["0","1"], default="1")

	args = parser.parse_args()
	if os.path.exists(args.output) == False:
		os.makedirs(args.output)
	chrs = []
	genomeFile = args.genome
	for record in SeqIO.parse(args.genome, "fasta"):
		dataset = pd.DataFrame()
		if "_" not in record.id and "chrM" not in record.id and "chrE" not in record.id:
			#doData([record.id, len(record.seq)])
			chrs.append([record.id, len(record.seq)])
	pool = mp.Pool(processes=args.nproc)
	pool.map(doData, chrs, chunksize=1)
