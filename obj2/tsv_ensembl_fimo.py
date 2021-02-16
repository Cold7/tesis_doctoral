import multiprocessing as mp
from glob import glob
from os import system 
import sys

def doAll(id):
	initT = tfbsSites[id]["init"]
	endT = tfbsSites[id]["init"]
	chrT = tfbsSites[id]["chr"]
	motifsInTFBS = []
	
	for data in motifs:
		if initT < data[3] and endT > data[2] and chrT == data[1]:
			if data[0] not in motifsInTFBS:
				motifsInTFBS.append(data[0])
	#id chr init end tfbm type
	if motifsInTFBS != []:
		file = open(".chr_"+chr+"_"+id+".out","w")
		file.write(id+"\t"+tfbsSites[id]["chr"]+"\t"+str(tfbsSites[id]["init"])+"\t"+str(tfbsSites[id]["end"])+"\t"+str(motifsInTFBS)[1:-1]+"\t"+tfbsSites[id]["type"]+"\n")
		file.close()
if __name__ == "__main__":
	global tfbsSites, chr, motifs
	motifs = [] #will have sub vectors of [motif, chr, init, end]
	chr = sys.argv[1]
	nproc = 8
	#reading ensembl gff and making dictionary
	f = open("homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff","r")
	tfbsSites = {}

	for line in f:
		aux = line.split("\t")
		
		if (aux[2] == "TF_binding_site" or aux[2] == "enhancer") and ("G" not in aux[0] and "KI" not in aux[0]) and aux[0]==chr:
			id = aux[8].split("ID=")[1].split(":")[1].split(";")[0]
			tfbsSites[id] = {"chr":"chr"+aux[0], "init":int(aux[3]), "end":int(aux[4]), "TFBM":[], "type":aux[2]}
	f.close()

	#saving fimo for current chr
	fimo = open("fimo.tsv","r")
	
	for line in fimo:
		aux = line.split("\t")
		if len(aux) > 1:
			if "motif_id" not in line and aux[2] == "chr"+chr:
				data = [aux[1],aux[2], int(aux[3]), int(aux[4])]
				motifs.append(data)


	fimo.close()
	
	#merging ensembl and fimo data
	pool=mp.Pool(processes=nproc) #for multiprocessing
	pool.map(doAll, (list(tfbsSites.keys())))	
	
	files = glob(".chr_"+chr+"_*.out")
	toSave = open("chr"+chr+"_ensembl_fimo.tsv","w")
	for file in files:
		f = open(file,"r")
		for line in f:
			toSave.write(line)
		f.close()
		system("rm "+file)
		
