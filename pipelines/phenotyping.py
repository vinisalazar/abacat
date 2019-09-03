# import argparse
from pipelines import annotate
from bactools import CONFIG, pathways

"""
This module uses the annotate function to predict metabolic pathways using the
pathway.json file.
"""


def phenotype(input_, blast="blastx"):
    """
    Parse an annotation to get phenotype.
    :param input_: Contigs file
    :param blast: Type of blast. Default is blastx
    :return:
    """
    genome = annotate(input_, db="phenotyping", blast=blast)

    def pathway_genes(info=True):
        """
        Takes the phenotyping geneset records and checks them against the pathways object
        from the CONFIG module.
        :return: pathway genes, a dict containing pathways as keys and identified records as values.
        """
        pathway_genes = {}
        for k, v in pathways.items():
            pathway_genes[k] = []
            for gene in genome.geneset["phenotyping"]["records"]:
                desc = gene.description.split()[1].split(".")[1]
                if desc in v:
                    pathway_genes[k].append(gene)
            if info:
                print(f"Found {len(pathway_genes[k])} genes for {k};")

        return pathway_genes
    
    genome.geneset["phenotyping"]["pathways"] = pathway_genes()
    return genome


genome = phenotype("/Users/viniWS/Bio/bactools/data/synecho/GCA_900473895.1_N32/GCA_900473895.1_N32_genomic.fna")
target = "/Users/viniWS/Bio/bactools/data/synecho/"