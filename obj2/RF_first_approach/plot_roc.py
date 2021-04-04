import matplotlib.pyplot as plt

if __name__ == "__main__":
	fig, ax = plt.subplots()
	ax.set_aspect('equal', adjustable='box')
	
	ax.set_title("ROC space", fontsize = 14)
	plt.xlabel("False Positive Rate")
	plt.ylabel("True Positive Rate")
#	plt.xlabel("Recall")
#	plt.ylabel("Precision")	
	ax.plot([0, 1], [0, 1], linestyle='--', lw=1, color='r', label='Chance', alpha=.8)

	dict = {}
	dict["300_5"] =	[1, 1, 0.652]
	dict["200_5"] =	[1, 1, 0.652]
	dict["500_15"] = [1, 0.991, 0.654]
	dict["500_5"] =	[1, 1, 0.652]
	dict["300_10"] = [1, 1, 0.652]
	dict["500_10"] = [1, 0.999, 0.652]
	dict["50_10"] =	[1, 0.995, 0.653]
	dict["300_15"] = [1, 0.984, 0.655]
	dict["100_5"] =	[1, 1, 0.652]
	dict["200_10"] = [1, 0.999,0.652]
	dict["100_10"] = [1, 1, 0.652]
	dict["50_15"] =	[1, 0.983, 0.656]
	dict["100_15"] = [1, 0.985, 0.655]
	dict["50_5"] = [1, 0.994, 0.653]
	dict["500_20"] = [1, 0.984, 0.655]
	dict["200_15"] = [1, 0.991, 0.654]
	dict["20000_max"] = [0.811, 0.103, 0.963]	
	
	trees = ["50", "100","200","300","500","20000"]
	depths = ["5","10","15","max"]
	for tree in trees:
		for depth in depths:
			if tree+"_"+depth in dict:
				label = "number of trees: "+tree+", depth: "+depth
				ax.plot(dict[tree+"_"+depth ][1], dict[tree+"_"+depth ][0], marker="o", label=label)# lw=1, color='r', label='Chance', alpha=.8)
#	ax.plot(0, 0, marker="o", color='w')# lw=1, color='r', label='Chance', alpha=.8)
#	ax.plot(1, 1, marker="o", color='w')# lw=1, color='r', label='Chance', alpha=.8)


	ax.legend()
	plt.show()
#	print("saving figure...")
#	plt.savefig("roc.jpg")
