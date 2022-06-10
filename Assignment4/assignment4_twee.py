import csv
import numpy as np
import pandas as pd
import shutil as sh

df = pd.read_csv("output/final.csv", header = None)
df.columns = ["dir", "n50"]
df["n50"] =  pd.to_numeric(df["n50"])
i = np.max(list(df["n50"]))
best = df[df["n50"] == i]
best_n50 = best.iloc[0,1]
best_kmer = str(best.iloc[0,0])
result = f"N50: {best_n50} kmer: {best_kmer}."
x = []
x.append(result)
with open("output/final.csv", "a") as f:
    writer = csv.writer(f)
    writer.writerow(x)
source= r'{}/contigs.fa'.format(best_kmer)
destination= r"output/"
sh.copy(source,destination)
del_list = list(df['dir'])
for i in del_list:
    sh.rmtree(str(i))