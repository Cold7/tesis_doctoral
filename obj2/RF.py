from glob import glob
import pandas as pd
from sklearn.ensemble import RandomForestClassifier as rfc
from joblib import dump

import argparse   # arguments parser
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-n","--ntree",help="number of trees. default 500", default=500, type=int)
	parser.add_argument("-f","--file",help="file to use", required=True)
	parser.add_argument("-o","--outname",help="file to use", required=True)
	args = parser.parse_args()
	
	f = open(args.file,"r")
	header = f.readline()[:-1].split("\t")[:-1]
	f.close()
	X = pd.read_csv(args.file,sep="\t",usecols=header)
	y = pd.read_csv(args.file,sep="\t",usecols=["promoter"])

	criterias = ["gini","entropy"]
	for criteria in criterias:
		clf = rfc(n_estimators=args.ntree, oob_score = True, n_jobs = 28, random_state=500, class_weight='balanced', criterion=criteria)
		model = clf.fit(X,y)
		print("using file "+args.file+" with "+str(args.ntree)+" trees got an oob score of "+str(model.oob_score_))
		dump(model, args.outname+"_"+criteria+".pkl", compress=3)

