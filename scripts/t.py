import pandas as pd

cts = "/home/a/pub/ycq/output/1/01.mapping/day10/WT/4"
genes_length = "/home/a/big/ycq/db/GRCm39/ncbi_dataset/data/GCF_000001635.27/genes.length"
gtf = "/home/a/big/ycq/db/GRCm39/ncbi_dataset/data/GCF_000001635.27/genomic.gtf"
cts = pd.read_csv(cts, sep="\t", comment="#")
genes_length = pd.read_csv(genes_length, sep="\t", header=None)
genes_length.columns = ['gene', 'length']
gtf = pd.read_csv(gtf, sep="\t", header=None, comment="#")
gtf.columns = ['contig', 'db', 'feature', 'start', 'end', 'a', 'b', 'c', 'attr']
sub_gtf = gtf.query('feature == "exon"')
attr = sub_gtf['attr']

a = [x.split("; ")[0] for x in attr]
b = [x.replace("gene_id ", "").replace("\"", "") for x in a]
c = set(b)
genes1 = set(list(cts['Geneid']))
genes2 = set(list(genes_length['gene']))
num = 0
for i in c.difference(genes1):
    num += 1
    print(i)
    if num >= 100:
        break
print(len(c))
print(len(genes1))
print(len(genes2))
# print(list(c.difference(genes2))[3])