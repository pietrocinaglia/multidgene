# -*- coding: utf-8 -*-

import os

from muxgene import MUXGENE

###
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multiplex Gene Network (muxGeNe) inference from gene expression across human tissues.
# @url: https://github.com/pietrocinaglia/muxgene
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

genes_of_interest = ['BRCA1','MZT2B','TP53','DDB2','USP24','ANAPC15','PSMD14','COP1']
tissues_of_interest = ['Brain_Cortex','Brain_Hippocampus','Liver']#,'Thyroid','Stomach','Liver']

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("################################################")
print()

muxgene = MUXGENE(genes_of_interest, tissues_of_interest, verbose=True)

# Data processing
multiplex_network = muxgene.construct_multiplex_network()

# Plotting multiplex_network, if output_path None the figure will not be saved as file.
muxgene.plot(multiplex_network)

# Storing multiplex_network as set of edgelist, one per layer
#muxgene.save_as_files(multiplex_network,output_path=WORKSPACE)

print(">>> DONE <<<")
