from Bio import Entrez
def search(query):
	Entrez.email = 'XXXXXX@XXXXX.XX'
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
	dict = {"ENSMUSG00000098343":"Mir6240",
			"ENSMUSG00000026656":"Fcgr2b",
			"ENSMUSG00000030688":"Stard10",
			"ENSMUSG00000071291":"Zfp58",
			"ENSMUSG00000018796":"Acsl1",
			"ENSMUSG00000026193":"Fn1",
			"ENSMUSG00000037706":"Cd81",
			"ENSMUSG00000116972":"",
			"ENSMUSG00000083563":"Gm13340",
			"ENSMUSG00000032554":"Trf",
			"ENSMUSG00000040466":"Blvrb",
			"ENSMUSG00000002985":"Apoe",
			"ENSMUSG00000001249":"Hpn"}
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

