from glob import glob
import pandas as pd

if __name__ == "__main__":

	binSizes = ["500","200","100"]
	for bin in binSizes:
									      ####
		files = glob("*_"+bin+"/*.bin") ################ change the extension
		for file in files:			     ####
			print("working with "+file)
			df = pd.read_csv(file,sep="\t")
			print("\tshape "+str(df.shape))
			counts = df.nunique()
			to_del = [i for i,v in enumerate(counts) if v == 1]
			print("\tcolumns to del: "+str(to_del))
			#drop useless columns
			df.drop(to_del,axis=1, inplace = True)
			print("\tshape "+str(df.shape))
			print("\tdeleting duplicate rows")
			df.drop_duplicates(inplace=True)
			print("\tshape "+str(df.shape))
			df.to_csv(file+".clean", sep="\t", index=False)
