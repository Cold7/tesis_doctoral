from glob import glob
import argparse
import pandas as pd
import multiprocessing as mp
import os

def returnStrand(currentPosition, vector): 
	#a function to knwo orientation of the full vector for the args.orientation option
	diff = 10000000 # a big number, cause it will have the abs value of current position less position in vector of tss
	pos = 0 # the position with the closer value to 0 in diff
	for i in range(len(vector)):
		if abs(currentPosition - vector[i]) < diff:
			diff = abs(currentPosition - vector[i])
			pos = vector[i]


	if pos >= currentPosition:
		return 1
	else:
		return -1


def doFullDataset(tsv):
	df1 = pd.read_csv(tsv,sep="\t")

	TSSList = [] # to make a list with TSS positions

	if args.orientation:
		for i in range(len(df1["TSS"])):
			if df1["TSS"][i] == 1:
				TSSList.append(i)


	dict = {}
	sorted_name = [] #to write with ordered index
	#verctor names for left and the list wich will be filled with data of the correspondig part
	i = args.left_vector * -1
	while i < 0:
		for name in args.feature_list:
			dict[name+"_"+str(i)] = []
			sorted_name.append(name+"_"+str(i))
		i += 1
	#central vector
	for name in args.feature_list:
		dict[name+"_0"] = []
		sorted_name.append(name+"_0")

	#right vectors
	i = 1
	while i<= args.right_vector:
		for name in args.feature_list:
			dict[name+"_"+str(i)] = []
			sorted_name.append(name+"_"+str(i))
		i += 1
	sorted_name.append(args.Class)
	
	#class
	dict[args.Class] = []


	#looping over dataset, and selecting vector if class is present or not (1 or -1)
	for i in range(len(df1[args.Class])):
		if i >args.left_vector and i <(len(df1[args.Class])-args.right_vector):
			if df1[args.Class][i] == 1 or df1[args.Class][i] == -1:
				orinetation = 1
				if args.orientation:
					orientation = returnStrand(i, TSSList)
				#looping column names
				for name in args.feature_list:
					data = df1[name][i-args.left_vector:i+args.right_vector+1].values.tolist()
					j = -1 * args.left_vector
					k = args.right_vector +1
					while j<0:
						dict[name+"_"+str(j*orientation)].append(data[j+args.left_vector])
						j += 1
					dict[name+"_0"].append(data[args.right_vector])
					j = 1
					while j<k:
						dict[name+"_"+str(j*orientation)].append(data[j+args.right_vector])
						j += 1
				dict[args.Class].append(df1[args.Class][i])

	
	df2 = pd.DataFrame(dict, columns= sorted_name)
	fileName = tsv.split("/")[-1]
	df2.to_csv(args.output+"/fullDataset_"+fileName,sep="\t", index= None)
	
	return
	



if __name__ == "__main__":
	global args

	parser = argparse.ArgumentParser() 
	parser.add_argument("-f", "--folder", help="Path to folder with tsv files with central vectors", required = True)
	parser.add_argument("-o","--output", help="folder where results will be placed. Default: ./", default = "./")
	parser.add_argument("-O","--orientation", help="to consider or not if full vector must be oriented to TSS.", default = "no", choices=["yes","no"])
	parser.add_argument("-n","--nproc", help="Number of processors to use. Default: 1", type = int, default=1)
	parser.add_argument("-lv","--left_vector", help="Number of vector at left of the central vetor. Default: 100", type = int, default=100)
	parser.add_argument("-rv","--right_vector", help="Number of vector at right of the central vetor. Default: 100", type = int, default=100)
	parser.add_argument("-fl","--feature_list", nargs="*", help="a list with features to consider.", required = True)
	parser.add_argument("-C", "--Class", help="feature to define as class", required = True)
	args = parser.parse_args()



	if os.path.exists(args.output) == False:
		os.makedirs(args.output)
	files = glob(args.folder+"/*.tsv")

	pool = mp.Pool(processes = args.nproc)
	pool.map(doFullDataset, files, chunksize = 1)
