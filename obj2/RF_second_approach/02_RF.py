import argparse   # arguments parser
import pandas as pd
import numpy as np
import joblib
from sklearn import ensemble
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import confusion_matrix as cm


if __name__ == "__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-t","--ntree", help="number of trees. default 2000", default=2000, type=int)
	parser.add_argument("-p","--nproc", help="number of processors. default 16", default=16, type=int)
	parser.add_argument("-C","--Class", help="current class to predict", required=True)
	parser.add_argument("-f","--features",help="features to consider to create the dataframe", nargs="+", required=True)
	parser.add_argument("-ft","--file_train", help="file to train algorithm", required=True)
	parser.add_argument("-st","--separator_train", help="separator character for file to train algorithm. default=space", default = "space")
	parser.add_argument("-fT","--file_test", help="file to test algorithm", required=True)
	parser.add_argument("-sT","--separator_test", help="separator character for file to test algorithm. default=space", default ="space")
	parser.add_argument("-o","--output", help="folder where place final models. default ./", default="./")

	args = parser.parse_args()

	#opening train and test files

	separator = None;
	if args.separator_train == "space":
		separator = " "
	if args.separator_train == "tab":
		separator = "\t"
	X1 = pd.read_csv(args.file_train, sep = separator, dtype=np.float32)[args.features] #float32 cause sklearn require it
	y1 = pd.read_csv(args.file_train, sep = separator, dtype=np.float32)[[args.Class]] #double [ cause it requires a list

	if args.separator_test == "space":
		separator = " "
	if args.separator_test == "tab":
		separator = "\t"
	
	X2 = pd.read_csv(args.file_test, sep = separator, dtype=np.float32)[args.features]
	y2 = pd.read_csv(args.file_test, sep = separator, dtype=np.float32)[[args.Class]]


	#getting basic statistics for X1 and X2
	
	print("shape of X1:", X1.shape)
	print("number of value for each classes in y1:")
	print(y1[args.Class].value_counts())
	print("shape of X2:", X2.shape)
	print("number of value for each classes in y2:")
	print(y2[args.Class].value_counts())


	#criterions to train RF
	criterions = ["gini"]

	for algorithm in criterions:
		print("===== working with "+algorithm+" =====")
		model = ensemble.RandomForestClassifier(
				criterion=algorithm,
				n_estimators  = args.ntree,
				n_jobs = args.nproc)
		model.fit(X1,y1.values.ravel())
		name = args.output+"/"+args.file_train.split("/")[-1].replace(".tsv","_"+algorithm)+".pkl"
		joblib.dump(model,name, compress=3)

		mse = mean_absolute_error(y1.values.ravel(), model.predict(X1))

		print("Training set mean absolute error with "+algorithm+": "+str(mse)+"\n")
		#validation
		validation = model.predict(X2)
		mse = mean_absolute_error(y2.values.ravel(), validation)
		print("Validation set mean absolute error: "+str(mse)+"\n")
		tn, fp, fn, tp = cm(y2, validation).ravel()
		print("TN "+str(tn)+"\n")
		print("FP "+str(fp)+"\n")
		print("FN "+str(fn)+"\n")
		print("TP "+str(tp)+"\n")
		precision = 0
		recall = 0
		try:
			recall = tp/(tp+fn)
			print("TPR "+str(tp/(tp+fn))+"\n")
		except:
			print("TPR error\n")
		try:
			print("FPR "+str(fp/(fp+tn))+"\n")
		except:
			print("FPR error\n")
		try:
			precision = tp/(tp+fp)
			print("precision "+str(tp/(tp+fp))+"\n\n")
		except:
			print("precision error\n\n")
		try:
			print("F1 "+str(2*(precision*recall)/(precision + recall)))
		except:
			print("F1 error")
		print("\nimportance:\n")
		pos = 0
		for importance in model.feature_importances_:
			print(X1.columns[pos]+"\t"+str(importance)+"\n")
			pos += 1
