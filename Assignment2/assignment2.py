from Bio import Entrez
import sys
import multiprocessing as mp
import pickle
import argparse as ap


Entrez.api_key = "cc90e8524a1e6db189cc428e8ddb8a862208"
Entrez.email = 'h.reitsma@st.hanze.nl'

# Your script needs to analyze the XML of each of the references further to extract all the authors of the article.
# It should save the authors in a Python tuple and use the Pickle module to save it to the disk as 
# output/PUBMED_ID.authors.pickle where PUBMEDID is of course the pubmed ID of the article in question.

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
    returns the xml info of the pubmed_id to an xml file
    '''
    handle = Entrez.efetch(db='pubmed', id=pubmed_id, retmode='xml', rettype='Abstract')
    records = Entrez.read(handle)
    # records = Entrez.parse(handle)
    handle.close()

    return records

def extract_authors(pubmed_id):
    records = pubmed_id_to_xml(pubmed_id)
    # print("authors:", records.get("AU", "?"))
    # print(records['PubmedBookArticle'])
    # print(records['PubmedArticle'])
    print(records.keys())
    return records



def main():
    # pmid = sys.argv[1]
    pmid = '30049270'
    cpus = mp.cpu_count()
    references = get_citation_ids(pmid)
    with mp.Pool(cpus) as pool:
        pool.map(extract_authors, references[:1])

if __name__ == '__main__':
    # argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    # argparser.add_argument("-n", action="store",
    #                        dest="n", required=False, type=int, default=10,
    #                        help="Number of references to download concurrently.")
    # argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    # args = argparser.parse_args()
    # print("Getting: ", args.pubmed_id)
    main()

    # with open('output/'+str(pubmed_id)+'.xml', 'wb') as f:
    #     for line in records:
    #         f.write((line))