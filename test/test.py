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
OUTPUT_FOLDER = WORKSPACE + "test_output/"

inherent_tissues = [
    'Brain_Cortex',
    'Brain_Substantia_nigra',
    'Brain_Hippocampus',
    'Brain_Putamen_basal_ganglia',
    'Brain_Cerebellum',
    'Brain_Caudate_basal_ganglia'
]
related_tissues = [
    'Stomach',
    'Small_Intestine_Terminal_Ileum',
    'Colon_Sigmoid',
    'Skin_Sun_Exposed_Lower_leg',
    'Adrenal_Gland',
    'Heart_Left_Ventricle'
]
non_inherent_tissues = [
    'Adipose_Subcutaneous',
    'Adipose_Visceral_Omentum',
    'Artery_Aorta',
    'Artery_Coronary',
    'Artery_Tibial',
    'Lung'
]
#
genes_of_interest = [
    'SNCA',
    'LRRK2',
    'PRKN',
    'PINK1',
    'GBA',
    'VPS35',
    'ATP13A2',
    'PLA2G6',
    'FBXO7',
    'DNAJC6',
    'SYNJ1',
    'PARK7',
    'VPS13C',
    'RAB32',
    'RAB39B'
]
#############
###
# mode:str can be defined as:
# - 'equal_exp' (default) [edges refer to equal expressed genes between tissues (layers)]
# - 'diff_exp' [edges refer to differential expressed genes between tissues (layers)]
# - 'diff_co_exp' [edges refer to differential co-expressed pairs of genes between tissues (layers)]
###
mode = 'diff_co_exp' #<--- change here
tissues_of_interest = inherent_tissues #<--- change here
co_label = 'inherent_tissues' #<--- change here (in accordance with tissues_of_interest)
#
output_filename = 'multidimensional_network_' + str(mode) + "_" + str(co_label)
#############

print("################################################")
print("> TESTING <")
print("- Genes: " + str(genes_of_interest))
print("- Tissues: " + str(tissues_of_interest))
print("- Mode: " + str(mode))
print("################################################")
print()

multidgene = MultiDGeNe(genes_of_interest, tissues_of_interest, data_path=WORKSPACE+'app/data/', p_threshold_intralayer=0.05, p_threshold_interlayer=0.05, mode=mode, verbose=True)

# Data processing
multidimensional_network, edge_colors = multidgene.construct_multidimensional_network()

# Storing multidimensional_network as set of edgelist, one per layer
multidimensional_network.write_edge_list(output_path=OUTPUT_FOLDER+output_filename+'.txt', data=True)

# Plotting multidimensional_network
# multidimensional_network.plot()
# html = multidimensional_network.generate_3d_plotly(output_path=OUTPUT_FOLDER+output_filename+'.html')
# multidimensional_network.plot(output_path=OUTPUT_FOLDER+output_filename, format='eps')
# multidimensional_network.plot(output_path=OUTPUT_FOLDER+output_filename, format='png', dpi=500)

print(">>> DONE <<<")
