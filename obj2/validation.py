from glob import glob
import pandas as pd
from joblib import load
from sklearn import metrics

import argparse   # arguments parser
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-m","--model",help="model to use in pkl format", required = True)
	parser.add_argument("-f","--file",help="file to use", required=True)
	parser.add_argument("-o","--outname",help="file to use", required=True)
	args = parser.parse_args()
	
	f = open(args.file,"r")
	header = f.readline()[:-1].split("\t")[:-1]
	f.close()
	X = pd.read_csv(args.file,sep="\t",usecols=header)
	y = pd.read_csv(args.file,sep="\t",usecols=["promoter"])

	model = load(args.model)
	pred = model.predict(X)
	fpr,tpr,thresholds = metrics.roc_curve(y,pred)
	auc = metrics.auc(fpr,tpr)
	acc = metrics.accuracy_score(y,pred)
	acc_balanced = metrics.balanced_accuracy_score(y,pred)
	f1 = metrics.f1_score(y,pred)
	p = metrics.precision_score(y,pred)
	r = metrics.recall_score(y,pred)
	tn, fp, fn, tp = metrics.confusion_matrix(y, pred).ravel()
	f = open(args.outname,"w")
	f.write(str(args.model)+"\t"+args.file+"\t"+str(tp)+"\t"+str(fp)+"\t"+str(fn)+"\t"+str(fp)+"\t"+str(tp/(tp+fn))+"\t"+str(fp/(fp+tn))+"\t"+str(p)+"\t"+str(f1)+"\t"+str(auc)+"\t"+str(acc)+"\t"+str(acc_balanced))
	f.close()

