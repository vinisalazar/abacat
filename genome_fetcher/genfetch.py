from bactools import GenomeFetch, Query

genfet = GenomeFetch("/Users/viniWS/storage/neorefs/rev6/test/test_rev6.fasta")

assemblies = dict()

for acc in genfet.accessions:
    assemblies[acc.repr] = Assembly(acc.out_fasta)

for k, v in assemblies.items():
    v.load_prodigal()
