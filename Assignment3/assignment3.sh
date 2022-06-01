#!/bin/bash
#SBATCH --time 7:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=week3Hendrik
#SBATCH --partition=assemblix

export BLASTDB=/local-fs/datasets/refseq_protein/
# export BLASTDB=/local-fs/datasets/
mkdir -p output
export time=/usr/bin/time
export time_file=output/timings.txt

for n in {15..16}
do
$time --append -o $time_file -f "${n}\t%e" blastp -query MCRA.faa -db refseq_protein/refseq_protein -num_threads $n
done

## Run python script
python3 assignment3.py






#!/bin/bash

# export BLASTDB=/local-fs/datasets/refseq_protein/
# mkdir -p /output/
# blastp -query MCRA.faa -db refseq_protein/refseq_protein -num_threads 1 -outfmt 6 >> blastoutput.txt

# declare -A core_time_dict
# for n in {16..16}
# do
# export start=`date +%s`
#     pycodestyle --ignore=E501 script1_martijn_vd_werff.py > script1_martijn_vd_werff.py.report
#     /usr/bin/time -o script1_martijn_vd_werff.py.times --append -f "${n}\t%e" python3 script1_martijn_vd_werff.py ./rnaseq.fastq -n $n -o script1_martijn_vd_werff.py.csv >> script1_james_gra.py.report 2>&1
# export end=`date +%s`
# export timing=`expr $end - $start`
# core_time_dict[n]=timing
# done

# ## Run python script
# python assignment3.py