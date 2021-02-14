import argparse   # arguments parser
from random import randint
from Bio import SeqIO
import numpy as np
import multiprocessing as mp
import os
import tempfile
from glob import glob


def parseBed(filename,base):
	featureDict = {}
	f = open(filename,"r")
	for line in f:
		line = line[:-1].split("\t")
		if line[0] not in featureDict:
			featureDict[line[0]] = []
		if  base == 1:
			featureDict[line[0]].append((int(line[1])-1,int(line[2])-1, float(line[4]))) #coords and score
		else:
			featureDict[line[0]].append((int(line[1]),int(line[2])-1, float(line[4])))
	f.close()

	return featureDict


def determineOrientation(middleCoord, chr):
	closerDist = float("inf") # the closer TSS, the infinite is to search the closer one so it will change when founding a closer one
	tssValue = -1
	for TSS in TSSpos[chr]:
		if abs(TSS-middleCoord)< closerDist:
			closerDist = abs(TSS-middleCoord)
			tssValue = TSS
	if tssValue < middleCoord:
		#this mean that the gene is upstream of the promoter and we need to change orientation
		return True
	else:
		#this mean that the tss is downstream of the promoter and we do not need to change the orientation
		return False

def describeBin(binCoords, chr):
	dictReturn = {}
	if args.DNase:
		currentValue = 0
		for data in dnaseqDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["DNase-seq"] = currentValue
	if args.dna_methylation:
		currentValue = 0
		for data in methDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["methylation"] = currentValue
	if args.ctcf:
		currentValue = 0
		for data in ctcfDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["CTCF"] = currentValue
	if args.yy1:
		currentValue = 0
		for data in yy1Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["YY1"] = currentValue
	if args.polr2a:
		currentValue = 0
		for data in polr2aDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["PolR2a"] = currentValue
	if args.ep300:
		currentValue = 0
		for data in ep300Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["EP300"] = currentValue
	if args.crebbp:
		currentValue = 0
		for data in cbpDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["CBP"] = currentValue
	if args.h3k27ac:
		currentValue = 0
		for data in h3k27acDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K27ac"] = currentValue
	if args.h3k9me3:
		currentValue = 0
		for data in h3k9me3Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K9me3"] = currentValue
	if args.h3k27me3:
		currentValue = 0
		for data in h3k27me3Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K27me3"] = currentValue
	if args.h3k4me3:
		currentValue = 0
		for data in h3k4me3Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K4me3"] = currentValue
	if args.h3k4me1:
		currentValue = 0
		for data in h3k4me1Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K4me"] = currentValue
	if args.h3k4me2:
		currentValue = 0
		for data in h3k4me2Dict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K4me2"] = currentValue
	if args.h3k9ac:
		currentValue = 0
		for data in h3k9acDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["H3K9ac"] = currentValue
	if  args.rna:
		currentValue = 0
		for data in rnaDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["RNAseq"] = currentValue
	if args.gene:
		currentValue = 0
		for data in geneDict[chr]:
			if binCoords[1] > data[0] and binCoords[0] < data[1]:
				currentValue = data[2]
		dictReturn["gene"] = currentValue
	return dictReturn

def doDataset(siteID):
	#the first step is looking if the current site have motif for the collection of chipseq
	motifs = sites[siteID][-1].split(",")
	motifCount = False
	for motif in motifs:
		if motif in tfPaths:
			motifCount = True
	if motifCount == False:
		return

	#so as we have motif and TF we evaluate the status of the site
	active = False
	for tf in tfPaths:
		if tf in sites[siteID][-1].split(","):
			file = open(tfPaths[tf],"r")
			for line in file:
				line = line[:-1].split("\t")
				if line[0] == sites[siteID][0]:
				#if end of the tf >= init of site  and init of the motif is <= that the end of the site
					if int(line[2]) >= int(sites[siteID][1]) and int(line[1]) <= int(sites[siteID][2]):
						active = True
			file.close()
	#creating a dict to save each site description, so first at all is to know orientation for the vector to the TSS
	#getting pseudo coords
	binSize = args.bins
	middleCoord = int(sites[siteID][1])
	middleCoord += int((int(sites[siteID][2])-int(sites[siteID][1]))/2)
	pseudoInit = middleCoord - int(binSize/2)
	pseudoEnd = middleCoord + int(binSize/2)


	#Now I will describe each fragment depending on the needChangeOrientation variable
	#the first step is to create a list up and down orientateless to the TSS
	#the orientation will be given after this step
	listOfBinCoords = []
	#upstream of the bin
	auxCoord = pseudoInit
	for i in range(args.flanking):
		auxInit = auxCoord - args.bins
		auxEnd = auxCoord
		listOfBinCoords.append((auxInit, auxEnd))
		auxCoord = auxInit
	listOfBinCoords.reverse()
	#middle bin
	listOfBinCoords.append((pseudoInit, pseudoEnd))
	#downstream of the central bin
	auxCoord = pseudoEnd
	for i in range(args.flanking):
		auxInit = auxCoord
		auxEnd = auxCoord + args.bins
		listOfBinCoords.append((auxInit, auxEnd))
		auxCoord = auxEnd
	#looking if the last coord is not bigger than the length of the current chr
	chr = sites[siteID][0]
	if listOfBinCoords[-1][1] < chrs[chr]: # if it is, we proceed, else we do not save nothing
		f = None
		#if site is active it will be displayed in the file title
		if active == False:
			f = open(args.temp+"/"+id+"_"+siteID+"_inactive","w")
		else:
			f = open(args.temp+"/"+id+"_"+siteID+"_active","w")
		#comparing the middle coord with TSS
		needChangeOrientation = determineOrientation(middleCoord, sites[siteID][0])
		#now we look for the orientation
		if needChangeOrientation == True:
			listOfBinCoords.reverse()

		#finally,  for each bin we will describe it and saving in a string
		fullVectorDict = {}
		fullVector = ""
		for i in range(len(listOfBinCoords)):
			auxDict = describeBin(listOfBinCoords[i],chr) #for each position it will be described as a dict of features
			for experiment in auxDict:
				fullVectorDict[experiment+"_"+str(i)] = auxDict[experiment]
		for i in range(len(listOfBinCoords)):
			if args.DNase:
				fullVector += str(fullVectorDict["DNase-seq_"+str(i)])+"\t"
			if args.dna_methylation:
				fullVector += str(fullVectorDict["methylation_"+str(i)])+"\t"
			if args.ctcf:
				fullVector += str(fullVectorDict["CTCF_"+str(i)])+"\t"
			if args.yy1:
				fullVector += str(fullVectorDict["YY1_"+str(i)])+"\t"
			if args.polr2a:
				fullVector += str(fullVectorDict["PolR2a_"+str(i)])+"\t"
			if args.ep300:
				fullVector += str(fullVectorDict["EP300_"+str(i)])+"\t"
			if args.crebbp:
				fullVector += str(fullVectorDict["CBP_"+str(i)])+"\t"
			if args.h3k27ac:
				fullVector += str(fullVectorDict["H3K27ac_"+str(i)])+"\t"
			if args.h3k9me3:
				fullVector += str(fullVectorDict["H3K9me3_"+str(i)])+"\t"
			if args.h3k27me3:
				fullVector += str(fullVectorDict["H3K27me3_"+str(i)])+"\t"
			if args.h3k4me3:
				fullVector += str(fullVectorDict["H3K4me3_"+str(i)])+"\t"
			if args.h3k4me1:
				fullVector += str(fullVectorDict["H3K4me_"+str(i)])+"\t"
			if args.h3k4me2:
				fullVector += str(fullVectorDict["H3K4me2_"+str(i)])+"\t"
			if args.h3k9ac:
				fullVector += str(fullVectorDict["H3K9ac_"+str(i)])+"\t"
			if  args.rna:
				fullVector += str(fullVectorDict["RNAseq_"+str(i)])+"\t"
			if args.gene:
				fullVector += str(fullVectorDict["gene_"+str(i)])+"\t"
		#and adding the current class
		if active == False:
			fullVector += "0"
		else:
			fullVector += "1"

		f.write(fullVector)
		f.close()
	return

if __name__ == "__main__":
	global id, args, chrs, TSSpos, sites, tfPaths, dnaseqDict, methDict, ctcfDict, yy1Dict, polr2aDict, ep300Dict, cbpDict, h3k27acDict, h3k9me3Dict, h3k27me3Dict, h3k4me3Dict, h3k4me1Dict, h3k4me2Dict, h3k9acDict, geneDict, rnaDict
	id = ""
	#generating a random ID
	chars = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","1","2","3","4","5","6","7","8","9","0"]
	for i in range(10):
		value = randint(0, len(chars)-1)
		id += chars[value]

	#input
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--genome", help="Path to genome file (in fasta format)", required = True)
	parser.add_argument("-gg","--gene_GTF", help="path for gene GTF (to get TSS).", required = True)
	parser.add_argument("-b","--bins", help="number of bases covering each bin. Default: 100", type = int, default=100)
	parser.add_argument("-fl","--flanking", help="number of vector up and downstream of the central bins. Default: 50", type = int, default=50)
	parser.add_argument("-o","--output", help="folder where results will be placed. Default: ./", default = "./")
	parser.add_argument("-na","--output_name", help=". Default: dataset", default = "dataset")
	parser.add_argument("-n","--nproc", help="Number of processors to use to compute TFBS active/inactive. Default: 1", type = int, default=1)
	parser.add_argument("-t","--temp", help="path for temporal folder", default=tempfile.gettempdir())

	####################################
	##
	## DNase-seq arguments
	##
	####################################
	parser.add_argument("-ds","--DNase",help="Option to generate features for DNAse.", action="store_true")
	parser.add_argument("-dsf","--DNase_file", help="DNase files.")
	parser.add_argument("-dsb","--DNase_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)

	####################################
	##
	## DNA methylation arguments
	##
	####################################
	parser.add_argument("-dm","--dna_methylation",help="Option to generate features for dna methylation.", action="store_true")
	parser.add_argument("-dmf","--dna_methylation_File", help="Folder where cpg DNAme files (in bed format) are located.")
	parser.add_argument("-dmb","--dna_methylation_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)
	parser.add_argument("-dc","--dna_CG",help="Option to generate percentage of C-G or G-C in the fragment.", action="store_true")

	####################################
	##
	## ctcf arguments
	##
	####################################
	parser.add_argument("-c","--ctcf",help="Option to generate features for ctcf.", action="store_true")
	parser.add_argument("-cFile","--ctcf_File", help="Folder where ctcf results (in bed file) are located.")
	parser.add_argument("-cb","--ctcf_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)

	####################################
	##
	## yy1 arguments
	##
	####################################
	parser.add_argument("-y","--yy1",help="Option to generate features for yy1.", action="store_true")
	parser.add_argument("-yFile","--yy1_File", help="Folder where yy1 results (in bed file) are located.")
	parser.add_argument("-yb","--yy1_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)

	####################################
	##
	## PolR2a arguments
	##
	####################################
	parser.add_argument("-pol","--polr2a",help="Option to generate features for polr2a.", action="store_true")
	parser.add_argument("-pFile","--polr2a_File", help="Folder where polr2a results (in bed file) are located.")
	parser.add_argument("-pb","--polr2a_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)

	####################################
	##
	## EP300 arguments
	##
	####################################
	parser.add_argument("-ep","--ep300",help="Option to generate features for ep300.", action="store_true")
	parser.add_argument("-epFile","--ep300_File", help="Folder where ep300 results (in bed file) are located.")
	parser.add_argument("-eb","--ep300_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)

	####################################
	##
	## cbp arguments
	##
	####################################
	parser.add_argument("-cbp","--crebbp",help="Option to generate features for cbp.", action="store_true")
	parser.add_argument("-cbpFile","--crebbp_File", help="Folder where cbp results (in bed file) are located.")
	parser.add_argument("-cbpb","--crebbp_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	## h3k27ac arguments
	##
	####################################
	parser.add_argument("-h1","--h3k27ac",help="Option to generate features for h3k27ac.", action="store_true")
	parser.add_argument("-h1File","--h3k27ac_File", help="Folder where h3k27ac	 results (in bed file) are located.")
	parser.add_argument("-h1b","--h3k27ac_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	## h3k9me3 arguments
	##
	####################################
	parser.add_argument("-h2","--h3k9me3",help="Option to generate features for h3k9me3.", action="store_true")
	parser.add_argument("-h2File","--h3k9me3_File", help="Folder where h3k9me3 results (in bed file) are located.")
	parser.add_argument("-h2b","--h3k9me3_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	##  h3k27me3 arguments
	##
	####################################
	parser.add_argument("-h3","--h3k27me3",help="Option to generate features for h3k27me3.", action="store_true")
	parser.add_argument("-h3File","--h3k27me3_File", help="Folder where h3k27me3 results (in bed file) are located.")
	parser.add_argument("-h3b","--h3k27me3_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	##  h3k4me3 arguments
	##
	####################################
	parser.add_argument("-h4","--h3k4me3",help="Option to generate features for h3k4me3.", action="store_true")
	parser.add_argument("-h4File","--h3k4me3_File", help="Folder where h3k4me3 results (in bed file) are located.")
	parser.add_argument("-h4b","--h3k4me3_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	##  h3k4me1 arguments
	##
	####################################
	parser.add_argument("-h5","--h3k4me1",help="Option to generate features for h3k4me1.", action="store_true")
	parser.add_argument("-h5File","--h3k4me1_File", help="Folder where h3k4me1 results (in bed file) are located.")
	parser.add_argument("-h5b","--h3k4me1_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	##  h3k4me2 arguments
	##
	####################################
	parser.add_argument("-h6","--h3k4me2",help="Option to generate features for h3k4me2.", action="store_true")
	parser.add_argument("-h6File","--h3k4me2_File", help="Folder where h3k4me2 results (in bed file) are located.")
	parser.add_argument("-h6b","--h3k4me2_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	##  h3k9ac arguments
	##
	####################################
	parser.add_argument("-h7","--h3k9ac",help="Option to generate features for h3k9ac.", action="store_true")
	parser.add_argument("-h7File","--h3k9ac_File", help="Folder where h3k9ac results (in bed file) are located.")
	parser.add_argument("-h7b","--h3k9ac_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)

	####################################
	##
	## Gene & RNA data
	##
	####################################
	parser.add_argument("-r","--rna",help="Option to generate features for RNA count.", action="store_true")
	parser.add_argument("-rf","--rna_File",help="Option to generate features for RNA. Input required is a tsv or bed file with 5 columns (chr init end strand and count)")
	parser.add_argument("-rb","--rna_base", help="to use 0-based or 1-based coord system. Default: 0",type=int,  choices=[0,1], default=0)
	parser.add_argument("-gen","--gene",help="Option to generate features for genes.", action="store_true")

	####################################
	##
	## Class arguments
	##
	####################################
	parser.add_argument("-C","--Class", help="name for the class.", required = True)
	parser.add_argument("-cf","--class_file", help="file with path to 5 columns chip-seq results for the current class.", required = True)
	parser.add_argument("-cd","--class_definition", help="tsv file with site definition (format is id, chr, init, end and motifs).", required = True)
	parser.add_argument("-Cb","--class_base", help="to use 0-based or 1-based coord system. Default: 0", type=int,  choices=[0,1], default=0)

	args = parser.parse_args()
	if os.path.exists(args.output) == False:
		os.makedirs(args.output)
	chrs = {}
	#generating TSS index and geneDict
	geneDict = {}
	f = open(args.gene_GTF,"r")
	TSSpos = {} #for each chr I will save TSSs
	for line in f:
		line = line[:-1].split("\t")
		if line[2] == "gene" and "gene_type \"protein_coding\"" in line[8]:
			if line[0] not in TSSpos:
				TSSpos[line[0]] = []
			if line[6] == "+":
				TSSpos[line[0]].append(int(line[3])-1)
			else:
				TSSpos[line[0]].append(int(line[4])-1)
			if  args.gene:
				if line[0] not in geneDict:
					geneDict[line[0]] = []
				geneDict[line[0]].append((int(line[3])-1,int(line[4])-1,1)) #the last 1 is an "artifitial score"
	f.close()

	#getting chr lengths
	for record in SeqIO.parse(args.genome, "fasta"):
		if "_" not in record.id and "chrM" not in record.id and "chrE" not in record.id:
			chrs[record.id] = len(record.seq)

	#getting site data
	sites = {}
	f = open(args.class_definition,"r")
	for line in f:
		line = line[:-1].split("\t")
		sites[line[0]] = line[1:]
	f.close()
	siteIDs = list(sites.keys())

	#getting promoter to use
	f = open(args.class_file,"r")
	tfPaths = {}
	for line in f:
		line = line[:-1].split("\t")
		tfPaths[line[0].upper()] = line[1]
	f.close()

	####################################################################
	##loading into memory the necessary files for a faster search
	####################################################################
	if args.DNase:
		dnaseqDict = parseBed(args.DNase_file,args.DNase_base)
	if args.dna_methylation:
		methDict = parseBed(args.dna_methylation_File, args.dna_methylation_base)
	if args.ctcf:
		ctcfDict = parseBed(args.ctcf_File,args.ctcf_base)
	if args.yy1:
		yy1Dict = parseBed(args.yy1_File,args.yy1_base)
	if args.polr2a:
		polr2aDict = parseBed(args.polr2a_File,polr2a_base)
	if args.ep300:
		ep300Dict = parseBed(args.ep300_File, args.ep300_base)
	if args.crebbp:
		cbpDict = parseBed(args.crebbp_File,args.crebbp_base)
	if args.h3k27ac:
		h3k27acDict = parseBed(args.h3k27ac_File,args.h3k27ac_base)
	if args.h3k9me3:
		h3k9me3Dict = parseBed(args.h3k9me3_File, args.h3k9me3_base)
	if args.h3k27me3:
		h3k27me3Dict = parseBed(args.h3k27me3_File, args.h3k27me3_base)
	if args.h3k4me3:
		h3k4me3Dict = parseBed(args.h3k4me3_File, args.h3k4me3_base)
	if args.h3k4me1:
		h3k4me1Dict = parseBed(args.h3k4me1_File, args.h3k4me1_base)
	if args.h3k4me2:
		h3k4me2Dict = parseBed(args.h3k4me2_File, args.h3k4me2_base)
	if args.h3k9ac:
		h3k9acDict = parseBed(args.h3k9ac_File, args.h3k9ac_base)
	if  args.rna:
		rnaDict = parseBed(args.rna_File, rna_base)

	#and making the vector for each site in a parallel fashion
	pool = mp.Pool(processes=args.nproc)
	pool.map(doDataset, siteIDs, chunksize=1)

	#once all sites are described we will join all temp, make simple statistics and delete temp files
	files = glob(args.temp+"/"+id+"*_active")
	print("there are "+str(len(files))+" active sites")
	files = glob(args.temp+"/"+id+"*_inactive")
	print("there are "+str(len(files))+" inactive sites")
	files = glob(args.temp+"/"+id+"*")
	outFile = open(args.output+"/"+args.output_name,"w")
	titleLine = ""

	for i in range((args.flanking*2)+1): #both side of the central vector plus the central vector
		if args.DNase:
			titleLine += "DNase-seq_"+str(i)+"\t"
		if args.dna_methylation:
			titleLine += "methylation_"+str(i)+"\t"
		if args.ctcf:
			titleLine += "CTCF_"+str(i)+"\t"
		if args.yy1:
			titleLine += "YY1_"+str(i)+"\t"
		if args.polr2a:
			titleLine += "PolR2a_"+str(i)+"\t"
		if args.ep300:
			titleLine += "EP300_"+str(i)+"\t"
		if args.crebbp:
			titleLine += "CBP_"+str(i)+"\t"
		if args.h3k27ac:
			titleLine += "H3K27ac_"+str(i)+"\t"
		if args.h3k9me3:
			titleLine += "H3K9me3_"+str(i)+"\t"
		if args.h3k27me3:
			titleLine += "H3K27me3_"+str(i)+"\t"
		if args.h3k4me3:
			titleLine += "H3K4me3_"+str(i)+"\t"
		if args.h3k4me1:
			titleLine += "H3K4me_"+str(i)+"\t"
		if args.h3k4me2:
			titleLine += "H3K4me2_"+str(i)+"\t"
		if args.h3k9ac:
			titleLine += "H3K9ac_"+str(i)+"\t"
		if  args.rna:
			titleLine += "RNAseq_"+str(i)+"\t"
		if args.gene:
			titleLine += "gene_"+str(i)+"\t"
	titleLine += args.Class+"\n"
	outFile.write(titleLine)
	for file in files:
		f = open(file,"r")
		l = f.readline()
		f.close()
		outFile.write(l+"\n")
		os.remove(file)
	print("Done!!")

