from glob import glob
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
	file = "gm12878_h1_hela_hepg2_k562.bin"
	f = open(file,"r")
	f1  = open("active","w")
	f2  = open("inactive","w")
	columns = ...
	for line in f:
		l = line
		line = line[:-1].split("\t")
		if line[-1] == "1.0":
			f1.write(l)
		elif line[-1] == "0.0":
			f2.write(l)
		else:
			f1.write(l)
			f2.write(l)
			columns = line[:-1]
	f.close()

	x = []
	pos = -5000
	while pos <= 5000:
		x.append(pos)
		pos+=100

	############################################
	##
	## uncomment the part that you wanna to plot
	##
	############################################
	############################ active part ###########################
#	df = pd.read_csv("active",sep="\t")
#	finalDict = {}
#	for column in columns:
#		c = column.split("_")[0]
#		finalDict[c] = []
#	for column in columns:
#		name = column.split("_")[0]
#		finalDict[name].append(sum(list(df[column]))/len(list(df[column])))
#		
#	fig, ax = plt.subplots()
#	for item in finalDict:
#		name = item
#		if item == "methylation":
#			name = "metilación"
#		ax.plot(x, finalDict[item], linestyle='-', lw=1,  label=name, alpha=.8)
#	plt.xlabel('Posición', fontsize=12)
#	plt.ylabel('Frecuencia', fontsize=12)
#	ax.legend()
#	plt.show()	
	############################ inactive part ###########################
	df = pd.read_csv("inactive",sep="\t")
	finalDict = {}
	for column in columns:
		c = column.split("_")[0]
		finalDict[c] = []
	for column in columns:
		name = column.split("_")[0]
		finalDict[name].append(sum(list(df[column]))/len(list(df[column])))
		
	fig, ax = plt.subplots()
	for item in finalDict:
		name = item
		if item == "methylation":
			name = "metilación"
		ax.plot(x, finalDict[item], linestyle='-', lw=1,  label=name, alpha=.8)
	plt.xlabel('Posición', fontsize=12)
	plt.ylabel('Frecuencia', fontsize=12)
	ax.legend()
	plt.show()	
