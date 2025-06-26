# -*- coding: utf-8 -*-

import os

from app.multidgene import MultiDGene

###
#
# [ This code is made available for testing purposes only ]
#
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference, bridging Co-Expression and Expression Dynamics across tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

genes_of_interest = ['BRCA1','MZT2B','TP53','DDB2','USP24','ANAPC15','PSMD14','COP1']
tissues_of_interest = ['Brain_Cortex','Brain_Hippocampus','Liver','Thyroid']#,'Stomach','Liver']
mode = 'diff_exp'

###
# mode:str can be defined as:
# - 'equal_exp' (default) [edges refer to equal expressed genes between tissues (layers)]
# - 'diff_exp' [edges refer to differential expressed genes between tissues (layers)]
# - 'diff_co_exp' [edges refer to differential co-expressed pairs of genes between tissues (layers)]
### 

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("- Mode: " + str(mode))
print("################################################")
print()

multidgene = MultiDGene(genes_of_interest, tissues_of_interest, data_path=WORKSPACE+'data/', p_threshold=0.05, mode=mode, verbose=True)

# Data processing
multidimensional_network, edge_colors = multidgene.construct_multidimensional_network()

# Plotting multidimensional_network
multidimensional_network.plot()

multidimensional_network.write_edge_list(WORKSPACE+"test.txt", data=True)

# Plotting multidimensional_network, if output_path None the figure will not be saved as file.
# multidgene.plot(multidimensional_network, edge_colors=edge_colors)

# Storing multidimensional_network as set of edgelist, one per layer
# multidgene.save_as_files(multidimensional_network, output_path=WORKSPACE)

print(">>> DONE <<<")
