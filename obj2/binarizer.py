from sklearn.preprocessing import Binarizer
from glob import glob
import pandas as pd


if __name__ == "__main__":
	binSizes = ["500","200","100"]
	for bin in binSizes:
		files = glob("*_"+bin+"/*")
		for file in files:
			dataset = pd.read_csv(file,sep="\t")
			columns = list(dataset.columns)
			dict = {}
			for column in columns:
				dict[column] = None

			transformer = Binarizer().fit(dataset)
			dataset = transformer.transform(dataset)
			for i in range(dataset.shape[1]):
				dict[columns[i]] = dataset[:,i]
			df = pd.DataFrame.from_dict(dict, orient="columns")
			df.to_csv(file+".bin", sep="\t", index=False)
		
