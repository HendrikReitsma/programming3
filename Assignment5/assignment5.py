import pyspark
import findspark
from pyspark.sql import SparkSession, SQLContext, Row
from pyspark import SparkContext, SparkConf
from pyspark.sql.types import StructType, StructField, IntegerType
from pyspark.sql.functions import avg, col, mean, avg, countDistinct, desc, split, explode
import pyspark.sql.functions
import csv

# ---------------------------------------------------------------Read the Data -------------------------------------------------------------------------
findspark.init()
findspark.find()
conf = pyspark.SparkConf().setAppName("hendrik").setMaster('local[16]')
sc = pyspark.SparkContext(conf=conf)
spark = SparkSession(sc)
sql_c = SQLContext(sc)
df = sql_c.read.csv('/data/dataprocessing/interproscan/all_bacilli.tsv',header=False, sep="\t")
df = df.toDF('protein_accession', 
              'Sequence', 
              'Sequence_len', 
              'Analysis', 
              'Signature_accession', 
              'Signature_description', 
              'Start', 
              'Stop', 
              'score', 
              'status', 
              'date',  
              'InterPro_accession', 
              'interPro_description', 
              'GO_annotations', 
              'Pathways_annotations')

df = df.withColumn('diff_length', ( df["Stop"] - df["Start"]) )


q1 = df.filter(df.InterPro_accession != '-').groupby('InterPro_accession').count().count()
q1_ex = str(df.groupby('InterPro_accession').count()._sc._jvm.PythonSQLUtils.explainString(df.groupby('InterPro_accession').count()._jdf.queryExecution(),'simple'))


q2 = df.filter(df.InterPro_accession != '-').groupBy('protein_accession').count().summary().collect()[1][2]
explain_2 = df.filter(df.InterPro_accession != '-').groupBy('protein_accession').count().summary()
q2_ex = str(explain_2._sc._jvm.PythonSQLUtils.explainString(explain_2._jdf.queryExecution(),'simple'))

go_terms = list()
for i in df.filter(df.GO_annotations != '-').select('GO_annotations').collect():
    go_terms.extend(i[0].split("|"))

sqlContext = SQLContext(sc)
rdd = sc.parallelize(go_terms)
ppl = rdd.map(lambda x: Row(name=x))
DF_go_terms = sqlContext.createDataFrame(ppl)


q3 = list()
for i in DF_go_terms.groupby('name').count().sort(desc('count')).select('name').head(1):
    q3.append(i[0])
explain_3 = DF_go_terms.groupby('name').count().sort(desc('count'))
q3_ex = str(explain_3._sc._jvm.PythonSQLUtils.explainString(explain_3._jdf.queryExecution(),'simple'))


q4 = float(df.select('diff_length').summary().collect()[1][1])
explain_4 = df.select('diff_length').summary()
q4_ex = str(explain_4._sc._jvm.PythonSQLUtils.explainString(explain_4._jdf.queryExecution(),'simple'))


q5 = list()
for i in df.filter(df.InterPro_accession != '-').groupby('InterPro_accession').count().sort(desc('count')).select('InterPro_accession').head(10):
    q5.append(i[0])
explain_5 = df.filter(df.InterPro_accession != '-').groupby('InterPro_accession').count().sort(desc('count'))
q5_ex = str(explain_5._sc._jvm.PythonSQLUtils.explainString(explain_5._jdf.queryExecution(),'simple'))


q6 = list()
for i in df.filter(df.InterPro_accession != '-').filter(df.diff_length / df.Seq_len > 0.9).sort(desc('diff_length')).select('InterPro_accession').head(10):
    q6.append(i[0])
explain_6 =df.filter(df.InterPro_accession != '-').filter(df.diff_length / df.Seq_len > 0.9).sort(desc('diff_length'))
q6_ex = str(explain_6._sc._jvm.PythonSQLUtils.explainString(explain_6._jdf.queryExecution(),'simple'))


q7 = list()
for i in df.filter(df.interPro_description != '-').withColumn('word',explode(split(col('interPro_description'), ' '))).groupby('word').count().sort(desc('count')).head(10):
    q7.append(i[0])
explain_7 = df.filter(df.interPro_description != '-').withColumn('word',explode(split(col('interPro_description'), ' '))).groupby('word').count().sort(desc('count'))
q7_ex = str(explain_7._sc._jvm.PythonSQLUtils.explainString(explain_7._jdf.queryExecution(),'simple'))


q8 = list()
for i in df.filter(df.interPro_description != '-').withColumn('word',explode(split(col('interPro_description'), ' '))).groupby('word').count().sort('count', descending ='False').head(10):
    q8.append(i[0])
explain_8 =df.filter(df.interPro_description != '-').withColumn('word',explode(split(col('interPro_description'), ' '))).groupby('word').count()
q8_ex = str(explain_8._sc._jvm.PythonSQLUtils.explainString(explain_8._jdf.queryExecution(),'simple'))


top_10_inter = df.filter(df.interPro_description != '-').filter(df.diff_length / df.Seq_len > 0.9).sort(desc('diff_length'))
q9 = list()
for i in top_10_inter.filter(top_10_inter.interPro_description != '-').withColumn('word',explode(split(col('interPro_description'), ' '))).groupby('word').count().sort(desc('count')).head(10):
     q9.append(i[0])

explain_9 = top_10_inter.filter(top_10_inter.interPro_description != '-').withColumn('word',explode(split(col('interPro_description'), ' '))).groupby('word').count().sort(desc('count'))
q9_ex = str(explain_9._sc._jvm.PythonSQLUtils.explainString(explain_9._jdf.queryExecution(),'simple'))


q10 = df.filter(df.InterPro_accession != '-').groupby(['protein_accession', 'Seq_len']).count().corr('count', 'Seq_len')
explain_10 =  df.filter(df.InterPro_accession != '-').groupby(['protein_accession', 'Seq_len']).count()
q10_ex = str(explain_10._sc._jvm.PythonSQLUtils.explainString(explain_10._jdf.queryExecution(),'simple'))

sc.stop()

header = ['Questions', 'Answer', '.explain']
answers = [
    [1, q1, q1_ex],
    [2, q2, q2_ex],
    [3, q3, q3_ex],
    [4, q4, q4_ex],
    [5, q5, q5_ex],
    [6, q6, q6_ex],
    [7, q7, q7_ex],
    [8, q8, q8_ex],
    [9, q9, q9_ex],
    [10, q10, q10_ex]

]

with open('results.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(answers)
