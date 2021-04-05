import networkx as nx
from pyloto import pyloto
from pyloto import REC
from pyloto import WREC
from pyloto import RGD

if __name__ == "__main__":
	n = 4 #number of processors to use
	dict = {}
	dict["consensus_intestino_control.gml"] = ["consensus_intestino_colitis.gml"]
	nodes = ["CAR1","HEXIM1","MYC","STAT3","BCL6","NFKB1","RELA","IRF4","IRF8","CD8A"]

	for item in dict:
		G = nx.read_gml(item)
		refLoTo =  pyloto(G, weight = "weight", nproc = n)
		for item2 in dict[item]:
			H =nx.read_gml(item2)

			infLoTo =  pyloto(H, weight = "weight", nproc = n)
			#getting Weighted REC (WREC)
			print("computing WREC for "+item+" and "+item2)
			myWREC = WREC(refLoTo, infLoTo,nproc = n) #weight 1 and 2 are the name for the weight attribute for network 1 and 2
			myWRGD = RGD(myWREC,list(refLoTo.network.nodes),nproc = n)#list to parse a list of string instead of nodeviews
			for node in nodes:
				if node in myWRGD:
					print(node+"\t"+str(myWRGD[node])+"\n")
