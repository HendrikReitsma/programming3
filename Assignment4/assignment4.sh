# !/bin/bash
# SBATCH --time 48:00:00
# SBATCH --cpus-per-task=16
# SBATCH --job-name=Hendrikweek4
# SBATCH --nodes=1
# SBATCH --mem=1000
# SBATCH --output=outfile
# SBATCH --partition=assemblix

file1=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R1_001_BC24EVACXX.filt.fastq
file2=/data/dataprocessing/MinIONData/MG5267/MG5267_TGACCA_L008_R2_001_BC24EVACXX.filt.fastq
mkdir -p output
seq 21 2 31 | parallel mkdir -p {}
seq 21 2 31 | parallel -j16 velveth {} {} -fastq -longPaired -separate $file1 $file2
seq 21 2 31 | parallel -j16 velvetg {}
seq 21 2 31 | parallel -j16 "python assignment4.py "
python assignment4_twee.py