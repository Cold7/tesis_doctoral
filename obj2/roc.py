from glob import glob
import joblib
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix as cm
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import auc

import warnings
warnings.filterwarnings("ignore")

def makeFeatures(file):
	f = open(file,"r")
	data = f.readline().split("\t")[:-1]
	f.close()
	return data


if __name__ == "__main__":
	modelList = None
	model = None
	testSet = None
	features = None
	X1 = None
	y1 = None
	test = None
	separator = "\t"

	features = makeFeatures("K562_100.tsv")
	X1 = pd.read_csv("K562_100.tsv", sep = separator, dtype=np.float32, usecols=features)
	y1 = pd.read_csv("K562_100.tsv", sep = separator, dtype=np.float32, usecols=["promoter"])

	fig, ax = plt.subplots()
	ax.set_aspect('equal', adjustable='box')
	fig.set_size_inches(11,8)
	print("\tloading K562 test set")
	features = makeFeatures("K562_100.tsv")
	print("\tLooking for models")
	modelList = glob("*.pkl")

	#ROC
	roc = None
	mean_fpr = np.linspace(0, 1, 100)
	tprs = []
	aucs = []
	for m in modelList:
		print("\t\tloading model: "+m)
		without = ""
		cells = {"GM12878":0,"HeLa-S3":0,"H1 hESC":0,"HepG2":0}
		dataSplited = m.split("_")
		for data in dataSplited:
			print(data)
			if data == "GM12878":
				cells["GM12878"] = 1
			elif data == "H1":
				cells["H1 hESC"] = 1
			elif data == "HepG2":
				cells["HepG2"] = 1
			elif data == "HeLa":
				cells["HeLa-S3"] = 1

		tree = "1000"
		name = "without "
		print(cells)
		for cell in cells:
			if cells[cell] == 0:
				name += cell+" "
		model = joblib.load(m)
		roc = plot_roc_curve(model, X1, y1, ax = ax, name = name, alpha = .3)

		#ROC
		interp_tpr = np.interp(mean_fpr, roc.fpr, roc.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)

	ax.set_title("ROC curves", fontsize = 14) 
	ax.plot([0, 1], [0, 1], linestyle='--', lw=1, color='r', label='Chance', alpha=.8)
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	ax.plot(mean_fpr, mean_tpr, color='b', label=r'Mean ROC (AUC = %0.2f)' % (mean_auc), lw=1, alpha=.8)
	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)

	ax.legend()

	print("saving figure...")
	plt.savefig("roc.jpg")
		
