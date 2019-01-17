from Bio import Entrez

Entrez.email = 'yanqi219@gmail.com'

def search(query):
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='10',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

def get_abstract(pmid):
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='text', rettype='abstract')
    return handle.read()

# if __name__ == '__main__':
#     results = search('Uridine diphosphate-N-acetylglucosamine autism')
#     id_list = results['IdList']
#     papers = fetch_details(id_list)
#     for i, paper in enumerate(papers['PubmedArticle']): print(
#         "%d) %s" % (i + 1, paper['MedlineCitation']['Article']['ArticleTitle']))

with open("Controls_ExpoUnexp.txt", "r") as fd:
    metabolite_name = fd.read().splitlines()

target = '&autism'
query_name = [x + target for x in metabolite_name]

print(query_name)

f = open('AP biomarkers with autism.txt', 'w')

for i in range(1, len(query_name)):
    try:
        print(metabolite_name[i])

        print(metabolite_name[i], file=f)

        if __name__ == '__main__':
            results = search(query_name[i])
            id_list = results['IdList']
            id_length = len(id_list)
            print('Total number of papers:', file=f)
            print(id_length, file=f)
            papers = get_abstract(id_list)
            print(papers, file=f)
    except: pass

f.close()
