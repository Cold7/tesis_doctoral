import numpy as np
import networkx as nx
import multiprocessing as mp
import pandas as pd

def _parallel_RGD(node):
	RGD_ = 0
	total = 0
	for rec in REC_dict_:
		data = rec.split(",")
		if node in data: #if node in the graphlet
			RGD_ += REC_dict_[rec]
			total += 1
	if total != 0: #faster to ask if node is acting in any graphlet
		RGD_ /= total
	return [node, RGD_]
	
#function that compute the REC graphlet degree based on a dict of graphlet:REC and a single node or a list of node.
#return the RGD value
def RGD(REC_dict, node, nproc = 1):
	finalRGD = {}
	if type(node) == list:
		global REC_dict_
		REC_dict_ = REC_dict
		pool = mp.Pool(processes = nproc)
		temp_rgd = pool.map(_parallel_RGD, node, chunksize = 1)
		for data in temp_rgd:
			finalRGD[data[0]] = data[1]
	else:
		RGD_ = 0
		total = 0
		for rec in REC_dict:
			data = rec.split(",")
			if node in data: #if node in the graphlet
				RGD_ += REC_dict[rec]
				total += 1
		if total != 0:
			RGD_ /= total
		finalRGD = {node: RGD_}
	return finalRGD

def _parallel_WREC(currentGraphlet):
	net1 = np.zeros((3,3))
	net2 = np.zeros((3,3))
	#filling net 1 and 2 matrices
	for i in range(3):
		j = i +1 
		while j < 3:
			#the higher weight will be 1 and the lower will be (Wlower/Whigher)
			#net1
			if netRef1.has_edge(currentGraphlet[i],currentGraphlet[j]):
				net1[i][j] = netRef1.nodes[currentGraphlet[i]]["weight"]
			if netRef1.has_edge(currentGraphlet[j],currentGraphlet[i]):
				net1[j][i] = netRef1.nodes[currentGraphlet[j]]["weight"]
			#net2
			if netRef2.has_node(currentGraphlet[i]) and netRef2.has_node(currentGraphlet[j]):
				if netRef2.has_edge(currentGraphlet[i],currentGraphlet[j]):
					net2[i][j] = netRef2.nodes[currentGraphlet[i]]["weight"]
				if netRef2.has_edge(currentGraphlet[j],currentGraphlet[i]):
					net2[j][i] = netRef2.nodes[currentGraphlet[j]]["weight"]
			j += 1
	
	#comparing both matrices
	wrec = 0
	for i in range(3):
		for j in range(3):
			if i != j:
				minWeight = 1
				if net1[i][j]< net2[i][j]:
					minWeight = net1[i][j]/net2[i][j]
				if net1[i][j]> net2[i][j]:
					minWeight = net2[i][j]/net1[i][j]
				wrec += 1-minWeight
	wrec /= 6
	wrec = 1 - wrec
	return(currentGraphlet+[wrec])
						
#function that compute  the weighted graphlet reconstruction rate
#this function get two pyloto object and return a dict for graphlet and REC
#to compare
def WREC(referenceLoTo, comparingLoTo, nproc = 1):
	global netRef1, netRef2, max_min, min_, weight_1_, weight_2_

	netRef1 = referenceLoTo.network
	netRef2 = comparingLoTo.network

	pool = mp.Pool(processes = nproc)
	resultingWREC = pool.map(_parallel_WREC,referenceLoTo.graphlets[1]+referenceLoTo.graphlets[2]+referenceLoTo.graphlets[3]+referenceLoTo.graphlets[4]+referenceLoTo.graphlets[5]+referenceLoTo.graphlets[6]+referenceLoTo.graphlets[7]+referenceLoTo.graphlets[8]+referenceLoTo.graphlets[9]+referenceLoTo.graphlets[10]+referenceLoTo.graphlets[11]+referenceLoTo.graphlets[12]+referenceLoTo.graphlets[13], chunksize = 1)
	WREC_ = {}
	for data in resultingWREC:
		WREC_[data[0]+","+data[1]+","+data[2]] = data[3]
	return WREC_

def _parallel_REC(currentGraphlet):
	net1 = np.zeros((3,3))
	net2 = np.zeros((3,3))
	#filling net 1 and 2 matrices
	for i in range(3):
		j = i +1 
		while j < 3:
			#net1
			if netRef1.has_edge(currentGraphlet[i],currentGraphlet[j]):
				net1[i][j] = 1
			if netRef1.has_edge(currentGraphlet[j],currentGraphlet[i]):
				net1[j][i] = 1
			#net2
			if netRef2.has_node(currentGraphlet[i]) and netRef2.has_node(currentGraphlet[j]):
				if netRef2.has_edge(currentGraphlet[i],currentGraphlet[j]):
					net2[i][j] = 1
				if netRef2.has_edge(currentGraphlet[j],currentGraphlet[i]):
					net2[j][i] = 1
			j += 1
	
	#comparing both matrices
	rec = 0
	for i in range(3):
		for j in range(3):
			if i != j:
				rec += abs(net1[i][j] - net2[i][j])
	rec /= 6
	rec = 1 - rec
	return(currentGraphlet+[rec])
						

#function that compute  the graphlet reconstruction rate
#this function get two pyloto object and return a dict for graphlet and REC
#to compare
def REC(referenceLoTo, comparingLoTo, nproc = 1):
	global netRef1, netRef2, max_min, min_, weight_1_, weight_2_

	netRef1 = referenceLoTo.network
	netRef2 = comparingLoTo.network

	pool = mp.Pool(processes = nproc)
	resultingREC = pool.map(_parallel_REC,referenceLoTo.graphlets[1]+referenceLoTo.graphlets[2]+referenceLoTo.graphlets[3]+referenceLoTo.graphlets[4]+referenceLoTo.graphlets[5]+referenceLoTo.graphlets[6]+referenceLoTo.graphlets[7]+referenceLoTo.graphlets[8]+referenceLoTo.graphlets[9]+referenceLoTo.graphlets[10]+referenceLoTo.graphlets[11]+referenceLoTo.graphlets[12]+referenceLoTo.graphlets[13], chunksize = 1)
	REC_ = {}
	for data in resultingREC:
		REC_[data[0]+","+data[1]+","+data[2]] = data[3]
	return REC_

#private function to compute basic statics
def _basicStats(graph):
	nodes = list(graph)
	nGenes = len(nodes)
	TFs = []
	for node in graph.nodes():
		if graph.out_degree(node) > 0:
			TFs.append(node)
	nTFs =len(TFs)
	nEdges = len(graph.edges())
	selfLoops = len(list(nx.selfloop_edges(graph)))
	negativeLoops = (nGenes*(nGenes-1))	- nEdges + nGenes
	return [nTFs, nGenes, nEdges, negativeLoops]
	
#private function to compute graphlets in a parallel fashion
def _parallel_graphlets(data):
	node1, node2 = data
	graphlets = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],13:[]}
	for node3 in G.nodes():
		if node3 != node1 and node3 != node2:
			#G1
			if G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == False and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == False and nodes.index(node2)<nodes.index(node3):
				graphlets[1].append(node3)
			#G2
			elif G.has_edge(node1, node3) == False and G.has_edge(node3, node1) == True and G.has_edge(node2, node1) == False and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == False:
				graphlets[2].append(node3)
			#G3
			elif G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == False:
				graphlets[3].append(node3)
			#G4
			elif G.has_edge(node1, node3) == False and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == False and G.has_edge(node3, node2) == True and G.has_edge(node2, node3) == False and nodes.index(node1)<nodes.index(node3):
				graphlets[4].append(node3)
			#G5
			elif G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == False and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == True:
				graphlets[5].append(node3)
			#G6		
			elif G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == True  and nodes.index(node1)<nodes.index(node2): 
				graphlets[6].append(node3)
			#G7
			elif G.has_edge(node1, node3) == False and G.has_edge(node3, node1) == True and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == False:
				graphlets[7].append(node3)
			#G8
			if G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == True and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == False and nodes.index(node2)<nodes.index(node3):
				graphlets[8].append(node3)
			#G9
			if G.has_edge(node1, node3) == False and G.has_edge(node3, node1) == True and G.has_edge(node2, node1) == False and G.has_edge(node3, node2) == False and G.has_edge(node2, node3) == True:
				graphlets[9].append(node3)
			#G10
			if G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == True and G.has_edge(node2, node3) == False:
				graphlets[10].append(node3)
			#G11
			if G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == False and G.has_edge(node2, node1) == False and G.has_edge(node3, node2) == True and G.has_edge(node2, node3) == True and nodes.index(node2)>nodes.index(node3):
				graphlets[11].append(node3)
			#G12
			if G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == True and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == True and G.has_edge(node2, node3) == False and nodes.index(node1)<nodes.index(node2) and nodes.index(node1)<nodes.index(node3):
				graphlets[12].append(node3)
			#G13
			elif G.has_edge(node1, node3) == True and G.has_edge(node3, node1) == True and G.has_edge(node2, node1) == True and G.has_edge(node3, node2) == True and G.has_edge(node2, node3) == True and nodes.index(node1)<nodes.index(node2) and nodes.index(node2)<nodes.index(node3):
				graphlets[13].append(node3)		
	return [node1, node2, graphlets]
	
#function to retrive graphlet data from a pyloto object
def _graphlets(H, nproc):
	global G, nodes
	#deleting self edges
	G  = H
	nodes = list(G.nodes())
	pool = mp.Pool(processes = nproc)
	graphlets_ = pool.map(_parallel_graphlets, list(G.edges), chunksize = 1 )
	
	graphlets = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],13:[]}
	
	for data in graphlets_:
		node1 = data[0]
		node2 = data[1]
		for graphletType in data[2]:
			for node3 in data[2][graphletType]:
				graphlets[graphletType].append([node1,node2,node3])

	return graphlets


	
class pyloto:

	def __init__(self,G, weight=None, nproc=1): #weight is the name for the attribute to use as weight
		
		#self properties
		self.weight_ = weight
		
		#getting basic statics
		_stats = _basicStats(G)
		self.number_TFs = _stats[0]
		self.number_genes = _stats[1]
		self.true_edges = _stats[2]
		self.false_edges = _stats[3]

		#creating a copy of the network and removing self edges
		self.network = G
		self.network.remove_edges_from(nx.selfloop_edges(G))
		
		#computing graphlets
		self.graphlets = _graphlets(self.network,nproc)
		
		#if the network has not an associated weight
		if weight == None:
			self.weight_ = "weight"
			for node in self.network.nodes:
				self.network.nodes[node]["weight"] = 1

	"""
	function to count number of graphlets by type
	"""
	def count_graphlets(self):
		count_graphlets_ = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,13:0}
		for type_ in self.graphlets:
			count_graphlets_[type_] +=len(self.graphlets[type_])
		return count_graphlets_

	"""
	function that receive a node and return a list of
	how many times is present on determined type of graphlets
	"""
	def count(self, node):
		count_node_graphlets = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,13:0, "total":0}
		for graphletType in self.graphlets:
			count_per_graphletType = 0
			for data in self.graphlets[graphletType]:
				if node in data:
					count_per_graphletType += 1
			count_node_graphlets[graphletType] = count_per_graphletType
			count_node_graphlets["total"] += count_per_graphletType
		return count_node_graphlets
				
	"""
	function that compute if node is or not in graphlet.
	A boolean value will be returned
	"""
	def nog(self, node, nproc = 1):
		for graphletType in self.graphlets:
			count_ = _count(self.graphlets, node, graphletType, nproc)
			if count_ != 0:
				del count_
				return True
		return False
