import sys
from os import system as s

#a simple script that take as input an ENSEMBL promoter (chr, init, end, ID and feature type like Promoter
#then the script will open the fimo output for the current chr and will look for overlapping motifs
#once all steps are done, the script will save a temporary file, so you need to concatenate all temporary files to generate the anotation for the whole genome
if __name__ == "__main__":

	chr =  sys.argv[1]
	init = int(sys.argv[2])
	end = int(sys.argv[3])
	ID = sys.argv[4]
	featureType = sys.argv[5]
	fimo = open(chr+"/fimo.tsv","r")
	cont = 0

	for line in fimo:
		if cont != 0:
			aux = line[:-1].split("\t")
			if len(aux) > 1:
				if aux[2] == chr:
					motif = aux[0].split("_")[0]
					initFimo  = int(aux[3])
					endFimo =  int(aux[4])
					if init < endFimo and end > initFimo:
						s("echo "+motif+" >> "+ID+".temp")
		else:
			cont += 1

	s("sort --temporary-directory=./temp -ur "+ID+".temp >sort_"+ID)
	s("rm -rf "+ID+".temp")
	f2 = open("sort_"+ID,"r")
	finalMotifs = ""
	for line in f2:
		finalMotifs +=line[:-1]+","
	finalMotifs = finalMotifs[:-1]
	f2.close()
	s("rm -rf sort_"+ID)

	f = open("final/."+featureType+"_"+chr+"_"+ID,"w")
	f.write(ID+"\t"+chr+"\t"+str(init)+"\t"+str(end)+"\t"+finalMotifs)
	f.close()
