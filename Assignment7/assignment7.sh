#!/bin/bash

# Set the number of CPUs to use
CPUS=16

# Set the list of XML files to process
XML_FILES="$(printf '%s\n' /data/datasets/NCBI/PubMed/pubmed21n*.xml)"

# KEYWORDS="cancer therapy vaccine"

# Run the script on each file in parallel
echo "$XML_FILES" | parallel -j "$CPUS" 'python keyword_counter.py -k cancer therapy vaccine -f "{}"'
# echo "$XML_FILES" | parallel -j "$CPUS" 'python keyword_counter.py -k "$KEYWORDS" -f "{}"'

python sum_counts.py "/students/2021-2022/master/hreitsma_DSLS/assignment7"



