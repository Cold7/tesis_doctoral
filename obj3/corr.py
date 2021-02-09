import numpy as np
import multiprocessing as mp
def doCorr(i):
	gene1 = genes[i]
	v1 = exprDict[gene1]
	f = open("./temp/"+gene1+".corr","w")
	for j in range(len(genes)):
		if i != j:
			v2 = exprDict[genes[j]]
			r = round(np.corrcoef(v1,v2)[0][1],3)
			f.write(genes[j]+"\t"+str(r)+"\n")
	f.close()
	return

if __name__ == "__main__":
	global exprDict, genes
	f = open("CD4.tsv","r")
	f.readline()
	exprDict = {}
	for line in f:
		line = line[:-1].split("\t")
		vector = []
		for n in line[1:]:
			vector.append(int(n))
		exprDict[line[0]] = vector
	f.close()
	genes = list(exprDict)
	keynumbers = []
	for i in range(len(exprDict)):
		keynumbers.append(i)
	pool = mp.Pool(processes = 104)
	pool.map(doCorr, keynumbers, chunksize = 1)
