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
from sklearn.metrics import confusion_matrix

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

	#PR
	pr = None
	mean_recall = np.linspace(0,1, 100)
	recalls = []
	aucs_pr = []
	noSkill = []
	for m in modelList:
		print("\t\tloading model: "+m)
		tree = "1000"
		cells = {"GM12878":0,"HeLa-S3":0,"H1 hESC":0,"HepG2":0}
		dataSplited = m.split("_")
		for data in dataSplited:
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
		for cell in cells:
			if cells[cell] == 0:
				name += cell+" "		
		model = joblib.load(m)
		pr = plot_precision_recall_curve(model, X1, y1, ax=ax, name = name, alpha = .3)
		tn, fp, fn, tp = confusion_matrix(model.predict(X1),y1).ravel()
		noSkill.append(tp/(tp+tn))
		
		#PR
		interp_pr = np.interp(mean_recall, list(pr.recall)[::-1], list(pr.precision)[::-1])
		interp_pr[0] = 1.0
		recalls.append(interp_pr)

	ax.set_title("PR curves", fontsize = 14)
	mean_pr = np.mean(recalls, axis=0)
	mean_pr[0] = 1.0
	mean_pr[-1] = 0
	mean_auc_pr = auc(mean_recall, mean_pr)
	std_pr = np.std(recalls)
	prs_upper = np.minimum(mean_recall + std_pr, 1)
	prs_lower = np.maximum(mean_recall - std_pr, 0)
	ax.plot(mean_recall, mean_pr, color='b', label=r'Mean PR (AP = %0.2f)' % (mean_auc_pr),lw=1, alpha=.8)
	noSkillAvg = sum(noSkill)/len(noSkill)
	print("noSkillAvg",noSkillAvg)
	ax.plot([0, 1], [noSkillAvg, noSkillAvg], linestyle='--', lw=1, color='r', label='No skill', alpha=.8)

	ax.legend()

	print("saving figure...")
	plt.savefig("pr.jpg")
		
