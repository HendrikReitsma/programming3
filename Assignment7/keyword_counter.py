"""
keyword_counter.py
Script that counts keywords
Author: Hendrik Reitsma
"""
import argparse as ap
import xml.etree.ElementTree as ET
import pandas as pd

class keywordCounter:
    def __init__(self, keywords):
        self.keywords = set(keywords)
        pass

    def count_keywords(self, xml_file):
        root = ET.parse(xml_file).getroot()
        keywords = []
        for article in root.findall('PubmedArticle'):
            keywords += self.get_keywords(article)
        keyword_series = pd.Series(keywords, dtype='object')
        keyword_df = pd.DataFrame({'keyword': keyword_series.value_counts().index,
                                'count': keyword_series.value_counts().values})
        return keyword_df

    def get_keywords(self, article):
        """
        Extract the specified keywords from the given article element and return them as a list.

        Parameters:
            article (Element): The Element object representing the article.

        Returns:
            list: A list of specified keywords.
        """
        keywords = []
        keyword_list = article.find('MedlineCitation').find('KeywordList')
        if keyword_list is not None:
            keywords = [keyword.text for keyword in keyword_list.findall('Keyword') if keyword.text in self.keywords]
        return keywords

if __name__ == "__main__":
    argparser = ap.ArgumentParser(description="Script that parses the xml files")
    argparser.add_argument("-k", "--keywords", action="store", nargs='+',
                           dest="keywords", required=True,
                           help="List of keywords to count")
    argparser.add_argument("-f", "--file", action="store", nargs='+',
                           dest="file", required=True,
                           help="xml file to process")
    args = argparser.parse_args()

    xml_file = args.file
    specified_keywords = args.keywords
    keyword_counter = keywordCounter(specified_keywords)
    results = keyword_counter.count_keywords(xml_file[0])

    # total_keyword_df = pd.concat(results).groupby('keyword')['count'].sum().reset_index()
    count_keyword_df = results.groupby('keyword')['count'].sum().reset_index()

    # Filter the dataframe to include only the specified keywords
    specified_keyword_df = count_keyword_df[count_keyword_df['keyword'].isin(args.keywords)]

    # Save the counts for the specified keywords to a file
    output_file = f"/students/2021-2022/master/hreitsma_DSLS/assignment7/{xml_file[0][-17:]}_counts.csv"
    specified_keyword_df.to_csv(output_file, index=False, sep="\t")
