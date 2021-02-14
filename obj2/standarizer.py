from sklearn.preprocessing import StandardScaler
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

			transformer = StandardScaler().fit(dataset)
			dataset = transformer.transform(dataset)
			for i in range(dataset.shape[1]):
				dict[columns[i]] = dataset[:,i]
			for i in range(len(dict["promoter"])):
				if dict["promoter"][i] > 0:
					dict["promoter"][i] = 1
				else:
					dict["promoter"][i] = 0
			df = pd.DataFrame.from_dict(dict, orient="columns")
			df.to_csv(file+".stdv", sep="\t", index=False)
		
