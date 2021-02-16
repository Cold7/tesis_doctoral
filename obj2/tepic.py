from glob import glob
import multiprocessing as mp

def compare(siteID):
	tepicFiles = glob("tepic/*")
	tfs = []
	for file in tepicFiles:
		tf = file.split("/")[1].split("_")[0].split("-")[1]
		if tf.upper() in dict[siteID]["motifs"]: #then I will use this site
			tfs.append(tf.upper())

	if len(tfs) == 0:
		return
	#as site have a motif with a prediction, then I will  compare the status of the site with the prediction
	realStatus = dict[siteID]["status"]
	predictedStatus = "0"  #in the next code section, this will be changed to "1" if it is ocuppied
	for tf in tfs:
		name = glob("tepic/K562-"+tf+"*")[0]
		f = open(name,"r")
		f.readline() #to avoid 1rst line
		for line in f:
			lineold = line
			line = line[2:-1].split("\t")[0].split("-")
			init = int(line[0].replace(":",""))
			end = int(line[1])
			chr = "chr"+lineold.split(":")[0]
			if init < dict[siteID]["end"] and end >dict[siteID]["init"] and chr == dict[siteID]["chr"]:
				predictedStatus = "1"
		f.close()
	f= open("temp/"+siteID,"w")
	f.write(realStatus+"\t"+predictedStatus)
	f.close()

#K562-TAF1_TEPIC_12_15_20_18_55_09_431528336_Affinity.txt
#   MAX
#1:1000080-1000260  1.891487929e-05
#1:100028931-100029038  1.533952867e-07


if __name__ == "__main__":
	global dict
	dict = {}
	f = open("promoter_hocomoco_status_for_k562","r")
	for line in f:
		line = line[:-1].split("\t")
#		if line[2] == "chr1":
		dict[line[0]] = {"status":line[1], "chr":line[2], "init":int(line[3]),"end":int(line[4]), "motifs":line[5].split(",")}
	pool = mp.Pool(processes = 104)
	pool.map(compare, list(dict),  chunksize=1)
	files = glob("temp/*")
	tp = 0
	tn = 0
	fp = 0
	fn = 0
	for file in files:
		f = open(file,"r")
		real,predict = f.readline().split("\t")
		f.close()
		if real == "1":
			if predict == "1":
				tp += 1
			else:
				fn += 1
		if real == "0":
			if predict == "0":
				tn += 1
			else:
				fp += 1

	print("tp","fp","fn","tn")
	print(tp,fp,fn,tn)
	print("tpr", tp/(tp+fn))
	print("fpr",fp/(fp+tn))
	print("P",tp/(tp+fp))
