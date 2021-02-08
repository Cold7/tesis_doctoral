from Bio import Entrez
def search(query):
	Entrez.email = 'XXXXX@XXXX.XX'
	handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='20',
                            retmode='xml',
                            term=query)
	results = Entrez.read(handle)
	return results

def fetch_details(id):
    Entrez.email = 'scontreras@dlab.cl'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=id)
    results = Entrez.read(handle)
    return results
    
if __name__ == "__main__":
	dict = {"ENSMUSG00000027556":"Car1",
			"ENSMUSG00000048878":"Hexim1",
			"ENSMUSG00000053977":"Cd8a",
			"ENSMUSG00000074272":"Ceacam1",
			"ENSMUSG00000022508":"Bcl6",
			"ENSMUSG00000039220":"Ppp1r10",
			"ENSMUSG00000050272":"Dscam",
			"ENSMUSG00000013419":"Zfp651",
			"ENSMUSG00000030364":"Clec2h",
			"ENSMUSG00000040284":"Gzmg",
			"ENSMUSG00000035441":"Myo1d",
			"ENSMUSG00000076471":"Trbv14",
			"ENSMUSG00000013653":"1810065E05Rik",
			"ENSMUSG00000022263":"Trio",
			"ENSMUSG00000035775":"Krt20",
			"ENSMUSG00000076472":"Trbv15",
			}
	print("IBD")
	for geneID in dict:
		geneName = dict[geneID]
		result = search("((IBD) OR (inflammatory bowel disease) OR (Crohn's disease) OR (Crohn disease) OR (ulcerative colitis)) AND ((T cell) OR (Th1) OR (Th2) OR (Th17) OR (Treg)) AND (("+geneName+") OR ("+geneID+"))")
		id_list = result['IdList']
		for id in id_list:
			try:
				paper = fetch_details(id)
				print(geneID, geneName,paper["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"])
			except:
				pass
	print("\n")
	print("Cancer")
	for geneID in dict:
		geneName = dict[geneID]
		result = search("((colorectal cancer) OR (colitis-associated cancer)) AND ((T cell) OR (Th1) OR (Th2) OR (Th17) OR (Treg)) AND (("+geneName+") OR ("+geneID+"))")
		id_list = result['IdList']
		for id in id_list:
			try:
				paper = fetch_details(id)
				print(geneID, geneName,paper["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"])
			except:
				pass

