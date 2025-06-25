# -*- coding: utf-8 -*-

import pymnet
import pandas as pd
import scipy
import numpy as np
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import warnings

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

def fisher_z(r):
    return 0.5 * np.log((1 + r) / (1 - r))

def inferring_intralayer_edges(gene_symbol_ids_of_interest:list, gene_exp_by_tissue:pd.DataFrame, tissues_of_interest, multidimensional_network:pymnet.MultilayerNetwork=None, p_threshold:float=10e-8, s_threshold:list=[0.5,0.9], edge_colors:set={}, verbose:bool=True) -> tuple[pymnet.MultilayerNetwork, set]:

    if verbose == True:
        print("-> Inferring intra-layer edges...")

    if multidimensional_network is None:
        multidimensional_network = pymnet.MultilayerNetwork(aspects=1, directed=False, fullyInterconnected=True)

    gene_pairs = list(combinations(gene_symbol_ids_of_interest, 2))

    # Collecting all results for global correction
    all_pvals = []
    all_corrs = []
    all_keys = []

    for current_layer in tissues_of_interest:

        multidimensional_network.add_layer(current_layer)

        if verbose == True:
            print("-- Layer " + str(current_layer) + ") ...")

        for u, v in gene_pairs:

            u_data = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == u) & (gene_exp_by_tissue['tissueSiteDetailId'] == current_layer) ]
            v_data = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == v) & (gene_exp_by_tissue['tissueSiteDetailId'] == current_layer) ]
            
            if u_data.empty or v_data.empty:
                continue
            
            u_gexp = u_data.data.values[0]
            v_gexp = v_data.data.values[0]

            if len(u_gexp) < 3 or len(v_gexp) < 3:
                #if verbose == True:
                    #warnings.warn("Size of samples are empty or less than 3; skipped.")
                continue

            s, p = __statistical_analysis(u_gexp, v_gexp, mode='co-exp')

            all_pvals.append(p)
            all_corrs.append(s)
            all_keys.append((u, v, current_layer))

    # Apply global FDR correction
    _, padj_all, _, _ = multipletests(all_pvals, method='fdr_bh')

    # Filter and add edges
    for (u, v, current_layer), s, padj in zip(all_keys, all_corrs, padj_all):
        if padj < p_threshold and (s_threshold[0] < s < s_threshold[1]):
            multidimensional_network[u, v, current_layer] = padj
            edge_colors[((u, current_layer), (v, current_layer))] = 'blue'

    if verbose == True:
        print("- [OK]")

    return multidimensional_network, edge_colors


def inferring_interlayer_edges(gene_symbol_ids_of_interest:list, gene_exp_by_tissue:pd.DataFrame, tissues_of_interest, multidimensional_network:pymnet.MultilayerNetwork=None, p_threshold:float=10e-8, s_threshold:list=[0.5,0.9], mode:str="equal_exp", edge_colors:set={}, verbose:bool=True) -> tuple[pymnet.MultilayerNetwork, set]:
    ###
    # mode:str can be defined as:
    # - 'equal_exp' (default) [edges refer to equal expressed genes between tissues (layers)]
    # - 'diff_exp' [edges refer to differential expressed genes between tissues (layers)]
    # - 'diff_co_exp' [edges refer to differential co-expressed pairs of genes between tissues (layers)]
    ###

    if verbose == True:
        print("-> Inferring inter-layer edges...")
    
    gene_pairs = list(combinations(gene_symbol_ids_of_interest, 2))
    layer_pairs = list(combinations(tissues_of_interest, 2))

    all_pvals = []
    gene_pair_keys = []

    if mode == 'equal_exp' or mode == 'diff_exp':
        for gene in gene_symbol_ids_of_interest:
            for l1, l2 in layer_pairs:
                l1_l2_data = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['tissueSiteDetailId'].isin([l1,l2]) & (gene_exp_by_tissue['geneSymbol'] == gene)) ]
                l1_l2_data = l1_l2_data.drop("tissueSiteDetailId", axis='columns')

                # Modelling expression data as array for each layer (l1,l2)
                expr_l1 = l1_l2_data.iloc[0, 1:].dropna().values.tolist()
                expr_l2 = l1_l2_data.iloc[1, 1:].dropna().values.tolist()

                if len(expr_l1) < 3 or len(expr_l2) < 3:
                    #if verbose == True:
                        #warnings.warn("Size of samples are empty or less than 3; skipped.")
                    continue

                _, p = __statistical_analysis(expr_l1, expr_l2, mode='exp')

                all_pvals.append(p)
                gene_pair_keys.append((gene, l1, l2))

        # Applying global multiple testing correction
        _, padj_all, _, _ = multipletests(all_pvals, method='fdr_bh')
        for (gene, l1, l2), adj_p in zip(gene_pair_keys, padj_all):
            if mode == 'equal_exp' and adj_p > p_threshold:
                multidimensional_network[gene, gene, l1, l2] = adj_p
                edge_colors[((gene, l1), (gene, l2))] = 'green'
            elif mode == 'diff_exp' and adj_p < p_threshold:
                multidimensional_network[gene, gene, l1, l2] = adj_p
                edge_colors[((gene, l1), (gene, l2))] = 'green'
    
    # differential co-expression
    elif mode == 'diff_co_exp':
        for u, v in gene_pairs:
            for l1, l2 in layer_pairs:
                u_data_l1 = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == u) & (gene_exp_by_tissue['tissueSiteDetailId'] == l1) ]
                v_data_l1 = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == v) & (gene_exp_by_tissue['tissueSiteDetailId'] == l1) ]
                
                u_data_l2 = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == u) & (gene_exp_by_tissue['tissueSiteDetailId'] == l2) ]
                v_data_l2 = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == v) & (gene_exp_by_tissue['tissueSiteDetailId'] == l2) ]
                
                u_gexp_l1 = u_data_l1.data.values[0]
                v_gexp_l1 = v_data_l1.data.values[0]

                u_gexp_l2 = u_data_l2.data.values[0]
                v_gexp_l2 = v_data_l2.data.values[0]
                
                # Normalize size between gene expression values --start
                # l1
                if len(u_gexp_l1) > len(u_gexp_l1):
                    idx = np.random.choice(len(u_gexp_l1), len(v_gexp_l1), replace=False)
                    u_gexp_l1 = u_gexp_l1[idx]
                elif len(u_gexp_l1) < len(v_gexp_l1):
                    idx = np.random.choice(len(v_gexp_l1), len(u_gexp_l1), replace=False)
                    v_gexp_l1 = v_gexp_l1[idx]
                # l2
                if len(u_gexp_l2) > len(v_gexp_l2):
                    idx = np.random.choice(len(u_gexp_l2), len(v_gexp_l2), replace=False)
                    u_gexp_l2 = u_gexp_l2[idx]
                elif len(u_gexp_l2) < len(v_gexp_l2):
                    idx = np.random.choice(len(v_gexp_l2), len(u_gexp_l2), replace=False)
                    v_gexp_l2 = v_gexp_l2[idx]
                # Normalize size between gene expression values --end

                s1, _ = __statistical_analysis(u_gexp_l1, v_gexp_l1, mode='co-exp')
                s2, _ = __statistical_analysis(u_gexp_l2, v_gexp_l2, mode='co-exp')

                # Transform
                z1 = fisher_z(s1)
                z2 = fisher_z(s2)

                # Number of samples per tissue
                n1 = len(u_gexp_l1)
                n2 = len(u_gexp_l2)

                # Standard error
                se = np.sqrt(1/(n1 - 3) + 1/(n2 - 3))

                # z-score
                z_diff = (z1 - z2) / se

                # p-value
                p = 2 * (1 - scipy.stats.norm.cdf(abs(z_diff)))

                all_pvals.append(p)
                gene_pair_keys.append((u, v, l1, l2))

        # Applying global multiple testing correction
        _, padj_all, _, _ = multipletests(all_pvals, method='fdr_bh')
        for (u, v, l1, l2), adj_p in zip(gene_pair_keys, padj_all):
            multidimensional_network[u, v, l1, l2] = adj_p
            edge_colors[((u, l1), (v, l2))] = 'green'
    
    return multidimensional_network, edge_colors

##
# A gene is differential expressed between two tissues (layer) of interest
##
def __statistical_analysis(x, y, mode='exp'):
    s = None
    p = None
    if mode == 'exp':
        # Normality test
        _, p_norm_u = scipy.stats.shapiro(x)
        _, p_norm_v = scipy.stats.shapiro(y)

        if p_norm_u > 0.05 and p_norm_v > 0.05:
            # x and y are normally distributed: Welch's t-test
            _, p = scipy.stats.ttest_ind(x, y, equal_var=False)
        else:
            # otherwise: Mann–Whitney
            _, p = scipy.stats.mannwhitneyu(x, y, alternative='two-sided')

    elif mode == 'co-exp':
        s, p = scipy.stats.pearsonr(x=x, y=y, alternative='two-sided')

    return s, p