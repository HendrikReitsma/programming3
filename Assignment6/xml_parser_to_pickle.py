"""
Script that parses xml files to pickle
Author: Hendrik Reitsma
"""

import os
import glob
import argparse as ap
import multiprocessing as mp
import xml.etree.ElementTree as ET
import pandas as pd


def parse_pubmed_xml(file_path)
    """
    Parses a Pubmed XML file into a pandas DataFrame.
    
    Args:
        file_path (str): Path to the Pubmed XML file.
    
    Returns:
        pd.DataFrame: DataFrame with columns for main_author, authors, references, 
                      publication_year, and keywords.
    """
    articles = {}
    root = ET.parse(file_path).getroot()
    
    for article in root.findall('PubmedArticle'):
        pmid = article.find('MedlineCitation').find('PMID').text
        
        authors = extract_authors(article)
        try:
            main_author = authors[0]
        except:
            main_author = ''
        
        references = extract_references(article)
        publication_year = get_publication_year(article)
        keywords = get_keywords(article)
        
        articles[pmid] = (main_author, authors, references, publication_year, keywords)
    
    record_df = pd.DataFrame.from_dict(articles, orient='index', 
                                        columns=['main_author', 'authors', 'references', 
                                                 'publication_year', 'keywords'])

    return record_df

def extract_authors(article):
    """Extracts authors from the given article element.

    Args:
        article (Element): An Element object representing an article.

    Returns:
        list: A list of authors as strings.
    """
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

def extract_references(article):
    """
    Extract the PubMed IDs of the references for the given article element and return them as a list.

    Args:
        article (Element): The article element to extract the references from.

    Returns:
        list: A list of PubMed IDs of the references for the given article element.
    """
    references = []
    reference_list = article.find('PubmedData').find('ReferenceList')
    if reference_list is not None:
        for reference in reference_list.findall('Reference'):
            if reference.find('ArticleIdList') is not None:
                article_id_list = reference.find('ArticleIdList')
                if article_id_list.find('ArticleId') is not None:
                    ref = article_id_list.find('ArticleId').text
                references.append(ref)

    return references

def get_publication_year(article):
    """
    Extract the publication year from the given article element and return it as a string.

    Args:
        article (Element): An ElementTree Element representing an article in PubMed XML format.

    Returns:
        str: The publication year as a four-digit string, or None if it cannot be extracted.

    """
    pub_date = article.find('MedlineCitation').find('Article').find('Journal').find('JournalIssue').find('PubDate')
    if pub_date is not None:
        year_element = pub_date.find('Year')
        if year_element is not None:
            year = year_element.text
            return year
    return None


def get_keywords(article):
    """
    Extract the keywords from the given article element and return them as a list.

    Parameters:
        article (Element): The Element object representing the article.

    Returns:
        list: A list of keywords.
    """
    keywords = []
    keyword_list = article.find('MedlineCitation').find('KeywordList')
    if keyword_list is not None:
        for keyword in keyword_list.findall('Keyword'):
            keywords.append(keyword.text)
    return keywords


def process_file(xml_file):
    """
    Extract relevant information from a Pubmed XML file, and save it as a pickled pandas DataFrame.

    Args:
        xml_file (str): The path to the input XML file.

    Returns:
        None
    """
    # Parse the XML file and extract the relevant information
    record_df = parse_pubmed_xml(xml_file)
    # Write the data to a pickle file
    pickle_file = os.path.splitext(xml_file)[0] + ".pkl"
    record_df.to_pickle(f'picklefiles/{pickle_file}')

def run_mp(files, args):
    # Multiprocessing function
    cpus = args.cpu
    with mp.Pool(cpus) as pool:
        results = pool.map(process_file, files)

if __name__ == "__main__":
    argparser = ap.ArgumentParser(description=
                                "Script that parses the xml files")
    argparser.add_argument("-cpu", action="store",
                            dest="cpu", required=True, type=int,
                            help="Number of cpus")
    args = argparser.parse_args()
    xml_files = glob.glob("/data/datasets/NCBI/PubMed/*.xml")
    run_mp(xml_files, args)