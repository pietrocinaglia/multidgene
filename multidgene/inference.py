# -*- coding: utf-8 -*-

import pymnet
import pandas as pd
import scipy
import numpy as np
from itertools import combinations
import warnings

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

def inferring_intralayer_edges(gene_symbol_ids_of_interest:list, gene_expression_tpm_by_tissue:pd.DataFrame, tissues_of_interest, multidimensional_network:pymnet.MultilayerNetwork=None, p_threshold:float=10e-8, s_threshold:list=[0.5,0.9], edge_colors:set={}, verbose:bool=True) -> tuple[pymnet.MultilayerNetwork, set]:

    if verbose == True:
        print("-> Inferring intra-layer edges...")

    if multidimensional_network is None:
        multidimensional_network = pymnet.MultilayerNetwork(aspects=1, directed=False, fullyInterconnected=True)

    gene_pairs = list(combinations(gene_symbol_ids_of_interest, 2))

    for current_layer in tissues_of_interest:

        multidimensional_network.add_layer(current_layer)

        if verbose == True:
            print("-- Layer " + str(current_layer) + ") ...")

        for u, v in gene_pairs:

            u_data = gene_expression_tpm_by_tissue.loc[ (gene_expression_tpm_by_tissue['geneSymbol'] == u) & (gene_expression_tpm_by_tissue['tissueSiteDetailId'] == current_layer) ]
            v_data = gene_expression_tpm_by_tissue.loc[ (gene_expression_tpm_by_tissue['geneSymbol'] == v) & (gene_expression_tpm_by_tissue['tissueSiteDetailId'] == current_layer) ]
            
            if (u_data.empty or v_data.empty):
                continue
            
            u_gexp = u_data.data.values[0]
            v_gexp = v_data.data.values[0]

            if len(u_gexp) < 3 or len(v_gexp) < 3:
                #if verbose == True:
                    #warnings.warn("Size of samples are empty or less than 3; skipped.")
                continue

            s, p = scipy.stats.pearsonr(x=u_gexp, y=v_gexp, alternative='two-sided')

            if p < p_threshold and (s > s_threshold[0] and s < s_threshold[1]):
                multidimensional_network[u,v,current_layer] = p_threshold
                edge_colors[((u, current_layer), (v, current_layer))] = 'blue'

    if verbose == True:
        print("- [OK]")

    return multidimensional_network, edge_colors


def inferring_interlayer_edges(gene_symbol_ids_of_interest:list, gene_read_counts_by_tissue:pd.DataFrame, tissues_of_interest, multidimensional_network:pymnet.MultilayerNetwork=None, p_threshold:float=10e-8, s_threshold:list=[0.5,0.9], edge_colors:set={}, verbose:bool=True) -> tuple[pymnet.MultilayerNetwork, set]:
        
    if verbose == True:
        print("-> Inferring inter-layer edges...")
    
    layer_pairs = list(combinations(tissues_of_interest, 2))

    for gene in gene_symbol_ids_of_interest:
        for l1, l2 in layer_pairs:

            u_v_data = gene_read_counts_by_tissue.loc[ (gene_read_counts_by_tissue['tissueSiteDetailId'].isin([l1,l2]) & (gene_read_counts_by_tissue['geneSymbol'] == gene)) ]
            u_v_data = u_v_data.drop("tissueSiteDetailId", axis='columns')

            # Model expression data as array for each layer (l1,l2)
            u_gexp = u_v_data.iloc[0, 1:].dropna().values.tolist()
            v_gexp = u_v_data.iloc[1, 1:].dropna().values.tolist()

            if len(u_gexp) < 3 or len(v_gexp) < 3:
                #if verbose == True:
                    #warnings.warn("Size of samples are empty or less than 3; skipped.")
                continue

            # Normality test
            _, p_norm_u = scipy.stats.shapiro(u_gexp)
            _, p_norm_v = scipy.stats.shapiro(v_gexp)

            if p_norm_u > 0.05 and p_norm_v > 0.05:
                # u_gexp and v_gexp are normally distributed: Welch's t-test
                _, p = scipy.stats.ttest_ind(u_gexp, v_gexp, equal_var=False)
            else:
                # otherwise: Mann–Whitney
                _, p = scipy.stats.mannwhitneyu(u_gexp, v_gexp, alternative='two-sided')

            if p < p_threshold:
                multidimensional_network[gene,gene,l1,l2] = p_threshold
                edge_colors[((gene, l1), (gene, l2))] = 'green'

    return multidimensional_network, edge_colors