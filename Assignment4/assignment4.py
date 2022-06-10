import csv
from Bio import SeqIO
import sys

dirs = list(sys.argv[1:])
kmer_list = []
row = []
for dir in dirs:
    with open(dir + "/contigs.fa") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = record.seq
            kmer_list.append(len(sequence))
    kmer_list.sort(reverse=True)
    num_list = kmer_list
    average = (sum(num_list))/2
    sum_num = 0
    for kmer in kmer_list:
        sum_num = sum_num + kmer
        if sum_num >= average:
            n50 = kmer
            break
    row.append([dir, n50])

with open("output/final.csv", "a") as f:
    writer = csv.writer(f)
    writer.writerows(row)
    