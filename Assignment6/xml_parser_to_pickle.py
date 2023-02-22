import pandas as pd
import lxml.etree as ET
import glob
from multiprocessing import Pool, cpu_count
import os

import glob
import pandas as pd
import argparse as ap
import multiprocessing as mp
import xml.etree.ElementTree as ET
from sqlalchemy import true

def parse_pubmed_xml(file_path):
    articles = {}
    root = ET.parse(file_path).getroot()
    # count = 0  # Counter to keep track of how many articles have been processed
    
    for article in root.findall('PubmedArticle'):
        # if count >= 5:  # If we have processed 5 articles, break out of the loop
        #     break
        pmid = article.find('MedlineCitation').find('PMID').text
        
        authors = get_authors(article)
        try:
            main_author = authors[0]
        except:
            main_author = ''
        # print('got authors')
        references = get_references(article)
        # print('got references')
        publication_year = get_publication_year(article)
        # print('got publication year')
        keywords = get_keywords(article)
        # print('got keywords')
        articles[pmid] = (main_author, authors, references, publication_year, keywords)
        # count += 1  # Increment the counter
        record_df = pd.DataFrame.from_dict(articles, orient='index', columns=['main_author', 'authors', 'references', 'publication_year', 'keywords'])
    return record_df

def get_authors(article):
    """Extract the authors from the given article element and return them as a list."""
    authors = []
    for author_list in article.find('MedlineCitation').find('Article').findall('AuthorList'):
        for author in author_list.findall('Author'):
            first_name = author.find('ForeName')
            last_name = author.find('LastName')
            if first_name is not None and last_name is not None:
                author_str = f"{first_name.text} {last_name.text}"
            else:
                author_str = "[No Name]"
            authors.append(author_str)
    return authors

def get_references(article):
    """Extract the PubMed IDs of the references for the given article element and return them as a list."""
    references = []
    reference_list = article.find('PubmedData').find('ReferenceList')
    # print(reference_list)
    if reference_list is not None:  # Check whether the References element exists
        for reference in reference_list.findall('Reference'):
            # Extract the PubMed ID of the reference from the PMID element
            if reference.find('ArticleIdList') is not None:
                article_id_list = reference.find('ArticleIdList')
                if article_id_list.find('ArticleId') is not None:
                    ref = article_id_list.find('ArticleId').text
                # Add the PubMed ID to the list
                references.append(ref)

    return references

def get_publication_year(article):
    pub_date = article.find('MedlineCitation').find('Article').find('Journal').find('JournalIssue').find('PubDate')
    if pub_date is not None:
        year_element = pub_date.find('Year')
        if year_element is not None:
            year = year_element.text
            return year
    return None

def get_keywords(article):
    keywords = []
    keyword_list = article.find('MedlineCitation').find('KeywordList')
    if keyword_list is not None:
        for keyword in keyword_list.findall('Keyword'):
            keywords.append(keyword.text)
    return keywords

def run_mp(files):
    # Function that does the multiprocessing.
    cpus = args.cpu
    with mp.Pool(cpus) as pool:
        results = pool.map(process_file, files)


def process_file(xml_file):
    # Parse the XML file and extract the relevant information
    record_df = parse_pubmed_xml(xml_file)
    # Write the data to a CSV file
    csv_file = os.path.splitext(xml_file)[0] + ".pkl"
    print(csv_file)
    record_df.to_csv(f'picklefiles/{csv_file}')

if __name__ == "__main__":

    argparser = ap.ArgumentParser(description=
                                "Script that parses the xml files to create pubmedID with keyword, author and reference lists and saves them as a pickle")
    argparser.add_argument("-cpu", action="store",
                            dest="cpu", required=True, type=int,
                            help="Number of cpus to run with multiprocessing")
    args = argparser.parse_args()

    xml_files = glob.glob("/data/datasets/NCBI/PubMed/*.xml")
    
    run_mp(xml_files)