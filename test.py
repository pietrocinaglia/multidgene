# -*- coding: utf-8 -*-

import os

# Own modules
from app.multidgene import MultiDGeNe

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
tissues_of_interest = ['Brain_Cortex','Brain_Hippocampus','Liver','Thyroid','Stomach','Pancreas']#,'Ovary']
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

multidgene = MultiDGeNe(genes_of_interest, tissues_of_interest, data_path=WORKSPACE+'app/data/', p_threshold=0.05, mode=mode, verbose=True)

# Data processing
multidimensional_network, edge_colors = multidgene.construct_multidimensional_network()

# Plotting multidimensional_network
#multidimensional_network.plot()

html = multidimensional_network.generate_3d_plotly(output_path=WORKSPACE+'multidimensional_network.html')

'''
@app.get("/network/3d", response_class=HTMLResponse)
def get_3d_network():
    html = designer.generate_3d_plotly()
    return HTMLResponse(content=html)
'''
# multidimensional_network.plot(output_path=WORKSPACE+'multidimensional_network', format='eps')
# multidimensional_network.plot(output_path=WORKSPACE+'multidimensional_network', format='png', dpi=500)

# Storing multidimensional_network as set of edgelist, one per layer
# multidimensional_network.write_edge_list(output_path=WORKSPACE+'multidimensional_network.txt', data=True)

print(">>> DONE <<<")
