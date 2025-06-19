# -*- coding: utf-8 -*-

import os

from multidgene import MultiDGene

###
#
# [ This code is made available for testing purposes only ]
#
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

genes_of_interest = ['BRCA1','MZT2B','TP53','DDB2','USP24','ANAPC15','PSMD14','COP1']
tissues_of_interest = ['Brain_Cortex','Brain_Hippocampus','Liver','Thyroid']#,'Stomach','Liver']

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("################################################")
print()

multidgene = MultiDGene(genes_of_interest, tissues_of_interest, verbose=True)

# Data processing
multidimensional_network, edge_colors = multidgene.construct_multidimensional_network()

# Plotting multidimensional_network, if output_path None the figure will not be saved as file.
multidgene.plot(multidimensional_network, edge_colors=edge_colors)

# Storing multidimensional_network as set of edgelist, one per layer
#multidgene.save_as_files(multidimensional_network,output_path=WORKSPACE)

print(">>> DONE <<<")
