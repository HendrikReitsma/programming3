#!/bin/bash
#SBATCH --time 7:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=week3Hendrik
#SBATCH --partition=assemblix

export BLASTDB=/local-fs/datasets/
mkdir -p output
export time=/usr/bin/time
export time_file=output/timings.txt

for n in {15..16}
do
$time --append -o $time_file -f "${n}\t%e" blastp -query MCRA.faa -db refseq_protein/refseq_protein -num_threads $n -outfmt 6 >> blastoutput.txt
done

## Run python script
python3 assignment3.py