import pyspark
import findspark
findspark.init()
findspark.find()
from pyspark.sql import SparkSession
from pyspark.sql import SQLContext
from pyspark import SparkContext, SparkConf
conf = pyspark.SparkConf().setAppName("hendrik_app").setMaster('local')
sc = pyspark.SparkContext(conf=conf)
spark = SparkSession(sc)
sql_c = SQLContext(sc)
from pyspark.sql.types import StructType, StructField, IntegerType
# schema = StructType([
#     StructField("protein_accession"),
#     StructField("sequence"),
#     StructField("sequence_length"),
#     StructField("analysis"),
#     StructField("signature_accession"),
#     StructField("signature_description"),
#     StructField("start"),
#     StructField("stop"),
#     StructField("score"),
#     StructField("status"),
#     StructField("date"),
#     StructField("interPro_annot_accession"),
#     StructField("interPro_annot_description"),
#     StructField("GO_annotations"),
#     StructField("Pathways_annotations"),
# ])
# df = sql_c.read.csv('/data/dataprocessing/interproscan/all_bacilli.tsv',schema=schema,header=False)

# df = sql_c.read.csv('/data/dataprocessing/interproscan/all_bacilli.tsv',header=False)
# df = df.toDF('protein_accession', 'sequence', 'sequence_length', 
# 'analysis', 'signature_accession', 'signature_description', 'start', 'stop', 
# 'score', 'status', 'date',  'interPro_annot_accession', 'interPro_annot_description', 
# 'GO_annotations', 'Pathways_annotations' 
# )

df = sql_c.read.csv('/data/dataprocessing/interproscan/all_bacilli.tsv',header=False, sep='\t')
# print(df.head())
print(len(df.columns))

