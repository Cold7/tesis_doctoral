import pandas as pd
from os import system
from glob import glob
import tempfile

import argparse   # arguments parser

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	#the p option is o know where are the tsv files of all chromosomes
	parser.add_argument("-p", "--path", help="Path to tsv files", required = True)
	parser.add_argument("-C","--Class", help="class name, default: TFBSac_ensembl", default="TFBSac_ensembl")
	parser.add_argument("-s","--size", help="original vector size. Default: 50", type = int, default = 50)
	parser.add_argument("-ud","--up_down", help="up and down bp size. Default: 2000", type = int, default = 2000)
	parser.add_argument("-o","--output", help="folder where results will be placed. Default: ./", default = "./")
	parser.add_argument("-t","--tmp", help="folder for temporal files. Default: "+tempfile.gettempdir(), default = tempfile.gettempdir())
	parser.add_argument("-n","--nproc", help="number of processors to use. Default: 1", default = 1, type = int)

	parser.add_argument("-ds","--DNase",help="use DNAse.", action="store_true")
	parser.add_argument("-dm","--dna_methylation",help="use dna methylation.", action="store_true")
	parser.add_argument("-c","--ctcf",help="use ctcf.", action="store_true")
	parser.add_argument("-y","--yy1",help="use yy1.", action="store_true")
	parser.add_argument("-e","--ep300",help="use ep300.", action="store_true")
	parser.add_argument("-h1","--H3K27ac",help="use H3K27ac.", action="store_true")
	parser.add_argument("-h2","--H3K27me3",help="use H3K27me3.", action="store_true")
	parser.add_argument("-h3","--H3K36me3",help="use H3K36me3.", action="store_true")
	parser.add_argument("-h4","--H3K4me1",help="use H3K4me1.", action="store_true")
	parser.add_argument("-h5","--H3K4me3",help="use H3K4me3.", action="store_true")
	parser.add_argument("-h6","--H3K9me3",help="use H3K9me3.", action="store_true")
	parser.add_argument("-G1","--gene1",help="Use gene as a single strand", action="store_true")
	parser.add_argument("-G2","--gene2",help="Use gene for plus and minus strand", action="store_true")
	parser.add_argument("-E1","--expression1",help="use expression as a single strand", action="store_true")
	parser.add_argument("-E2","--expression2",help="Use gene for plus and minus strand", action="store_true")

	#parsing arguments
	args = parser.parse_args()

	size = args.size
	inputFile = "output.tsv"# input("Input file to use: ")
	outFolder = args.output
	nproc = args.nproc

	className = args.Class


	colNames = []
	if args.DNase:
		colNames.append("DNase")
	if args.dna_methylation:
		colNames.append("methylation")
	if args.ctcf:
		colNames.append("ctcf")
	if args.yy1:
		colNames.append("yy1")
	if args.ep300:
		colNames.append("ep300")
	if args.H3K27ac:
		colNames.append("H3K27ac")
	if args.H3K27me3:
		colNames.append("H3K27me3")
	if args.H3K36me3:
		colNames.append("H3K36me3")
	if args.H3K4me1:
		colNames.append("H3K4me1")
	if args.H3K4me3:
		colNames.append("H3K4me3")
	if args.H3K9me3:
		colNames.append("H3K9me3")
	if args.gene1:
		colNames.append("gene")
	else:
		if args.gene2:
			colNames.append("gene_min")
			colNames.append("gene_plus")
	if args.expression1:
		colNames.append("expression")
	else:
		if args.expression2:
			colNames.append("expression_strand_min")
			colNames.append("expression_strand_plus")

	#preparing file to save the combinations
	Name = "up_down_"
	for col in colNames:
		Name += col+"-"
	Name += args.Class+".out"

	up_down_vector_num = args.up_down/args.size

	outFile = open(Name,"w")
	#writting col names
	i = 0
	init = i - up_down_vector_num
	end = i + up_down_vector_num
	name = ""
	while init <= end:
		for col in colNames:
			name += col+"_"+str(int(init))+"\t"
		init += 1
	name += "class\n"
	outFile.write(name)
	tsvFiles = glob("*.tsv")

	for file in tsvFiles:
		tsv = pd.read_csv(file, sep="\t")
		for i in range(len(tsv[className])):
			if tsv[className][i] == 1 or tsv[className][i] == -1:
				init = i - up_down_vector_num
				end = i + up_down_vector_num
				line = ""
				while init <= end:
					for col in colNames:
						line += str(tsv[col][init])+"\t"
					init += 1
				line += str(int(tsv[className][i]))+"\n"
				outFile.write(line)

	outFile.close()
	system("sort -u -r --temporary-directory="+args.tmp+" --parallel="+str(args.nproc)+" "+Name+" -o "+Name+".sorted")
	
	f = open(Name+".sorted","r")

	
	data = []
	
	
	for line in f:
		if "_" not in line:
			aux = line[:-1].split("\t")
			data2 = []
			init = 0
			end = len(colNames)
			while end < len(aux): #is < and not <= cause the last one is the class
				data2.append([aux[-1]] + aux[init:end])
				init = end
				end += len(colNames)
			data.append(data2)
	
	f1 = open(Name+".brnn","w")
	f1.write(str(len(data))+" "+str(len(colNames))+" 2\n")
	for item in data:
		f1.write(str(len(item))+"\n")
		for item2 in item:
			line = ""
			for item3 in item2:
				line +=item3+" "
			f1.write(line[:-1]+"\n")
	f.close()
	f1.close()
