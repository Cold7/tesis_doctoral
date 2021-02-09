from Bio import Entrez
def search(query):
	Entrez.email = 'scontreras@dlab.cl'
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
			"ENSMUSG00000050272":"Dscam"
			}
	print("IBD")
	for geneID in dict:
		geneName = dict[geneID]
		result = search("((IBD) OR (inflammatory bowel disease) OR (Crohn's disease) OR (Crohn disease) OR (ulcerative colitis)) AND ((CD4 T cell) OR (CD4+ T cell) OR (Th1) OR (Th2) OR (Th17) OR (Th9) OR (T fh) OR (Tfh) or OR (T folicular helper) OR (Treg)) AND (("+geneName+") OR ("+geneID+"))")
		id_list = result['IdList']
		for id in id_list:
			try:
				paper = fetch_details(id)
				print(geneID, geneName,paper["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleTitle"])
			except:
				pass

