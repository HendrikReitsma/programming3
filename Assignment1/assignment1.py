from Bio import Entrez
import sys
import multiprocessing as mp
Entrez.api_key = "cc90e8524a1e6db189cc428e8ddb8a862208"
Entrez.email = 'h.reitsma@st.hanze.nl'

def get_citation_ids(pmid):
    """
    Input: pubmed id
    :return: references
    """
    results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                   db="pmc",
                                   LinkName="pubmed_pmc_refs",
                                   id=pmid,
                                   api_key='cc90e8524a1e6db189cc428e8ddb8a862208'))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references

def pubmed_id_to_xml(pubmed_id):
    '''
    parameter: pubmed_id
    writes the abstract of the pubmed_id to an xml file
    '''
    # handle = Entrez.esummary(db="nucleotide", id=pubmed_id, rettype="Abstract", retmode="xml")
    handle = Entrez.efetch(db='pubmed', id=pubmed_id, retmode='xml', rettype='Abstract')
    records = handle.readlines()
    handle.close()
    with open('output/'+str(pubmed_id)+'.xml', 'wb') as f:
        for line in records:
            f.write((line))
    return True

def main():
    pmid = sys.argv[1]
    # pmid = '19304878'
    cpus = mp.cpu_count()
    # record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))
    ## Get the related ID's:
    # related_ids = []
    # results = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))
    # related_ids = [link["Id"] for link in results[0]["LinkSetDb"][5]["Link"]]
    # related_ids = related_ids[:10]
    # for i in range(1,11):
    #     related_id = record[0]["LinkSetDb"][5]["Link"][i]['Id']
    #     related_ids.append(int(related_id))
    # print(related_ids)
    references = get_citation_ids(pmid)
    with mp.Pool(cpus) as pool:
        pool.map(pubmed_id_to_xml, references[:10])
    print(references)

if __name__ == '__main__':
    main()