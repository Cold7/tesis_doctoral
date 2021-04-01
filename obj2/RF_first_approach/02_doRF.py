from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import argparse   # arguments parser
import multiprocessing as mp
from glob import glob

def textInLine(line):
	l = line.lower()
	characters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"]
	for char in l:
		if char in characters:
			return True
	return False
    
def confusion_matrix(y_true, y_pred, true_class = "1.0"):
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	
	for i in range(len(y_pred)):
		if y_pred[i] == y_true[i]:
			if y_pred[i] == true_class:
				tp += 1
			else:
				tn += 1
		else:
			if y_pred[i] == true_class:
				fp += 1
			else:
				fn += 1
	print(tp, fp, fn, tn)
	return [tp, fp, fn, tn]
    
def doRF(file,tree, depth, param):
	print ("working with file: "+file+" using tree "+str(tree)+" depth "+str(depth)+" and params: "+str(param)+"\n")
	f = open(str(tree)+"_"+str(depth)+"_"+str(param)+"_"+file.replace("./","")+".log","w")
	f.write("with file: "+file+" using tree "+str(tree)+" depth "+str(depth)+" and params: "+str(param)+"\n")
	#print ("n of tree: "+str(tree)+" depth: "+str(depth))
	clf = RandomForestClassifier(n_estimators=tree, max_depth=depth, max_features=param, n_jobs=60)
	clf.fit(X_train, y_train)
	f.write(str(clf.estimators_)+"\n")
	print(str(clf.estimators_)+"\n")
	toPrint = ""
	for num in clf.feature_importances_:
		toPrint+=str(num)+"\t"
	toPrint = toPrint[:-1]+"\n"
	f.write(str(toPrint))
	print(toPrint)
	f.write(str(clf.get_params())+"\n")
	print(str(clf.get_params())+"\n")
	predictions = []
	for item in X_test:
		pred = clf.predict([item])
		predictions.append(pred)
	
	tp, fp, fn, tn = confusion_matrix(y_test, predictions)
	f.write("TP: "+str(tp)+"\n")
	f.write("TN: "+str(tn)+"\n")
	f.write("FP: "+str(fp)+"\n")
	f.write("FN: "+str(fn)+"\n\n")

	print("TP: "+str(tp)+"\n")
	print("TN: "+str(tn)+"\n")
	print("FP: "+str(fp)+"\n")
	print("FN: "+str(fn)+"\n\n")
	
	r = tp/(tp+fn)
	p = tp/(tp+fp)
	f1 = (2*p*r)/(p+r)
	f.write("TPR: "+str(r)+"\n")
	f.write("FPR: "+str(fp/(fp+tn))+"\n")
	f.write("P  : "+str(p)+"\n")
	f.write("F1 : "+str(f1)+"\n")
	print("TPR: "+str(r)+"\n")
	print("FPR: "+str(fp/(fp+tn))+"\n")
	print("P  : "+str(p)+"\n")
	print("F1 : "+str(f1)+"\n")
	print("##########################################\n############################")
	f.write("################################################################")
	f.close()
    
if __name__ == "__main__":
	global X_train, X_test, y_train, y_test
	files = glob("./*sorted")
	#files = glob("./*")
	
	aux = []
	for file in files:
		aux.append(0-len(file))
	
	for i in range(len(files)):
		j = i+1
		while j < (len(files)):
			if aux[j]< aux[i]:
				temp = aux[i]
				temp2 =files[i]
				aux[i] = aux[j]
				aux[j] = temp
				files[i] = files[j]
				files[j] = temp2
			j += 1
		
	for f in files:
		#parsing arguments
		X = []
		y = []
		file = open(f,"r")
		for line in file:
			if textInLine(line) == False:
				descriptors = line[:-1].split("\t")[:-1]
				Class = line[:-1].split("\t")[-1]
				X.append(descriptors)
				y.append(Class)
		X_train, X_test, y_train, y_test = train_test_split(X,y, test_size = 0.3)
		
		
		numberOfTrees=[1, 50, 100, 200, 300, 500]
		depths = [ 5, 10, 15]
		nParams = [5, 10, 15, 20]
	
		for tree in numberOfTrees:
			for depth in depths:
				for param in nParams:
					try:
						print("using file ", f, "with tree, depth, and param (max features)",tree, depth, param)
						doRF(f, tree, depth, param)
					except:
						print("error using file ", f, "with tree, depth, and param (max features)",tree, depth, param)

