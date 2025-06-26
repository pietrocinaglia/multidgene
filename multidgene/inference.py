# -*- coding: utf-8 -*-

import pymnet
import pandas as pd
import scipy
import numpy as np
from itertools import combinations
from statsmodels.stats.multitest import multipletests

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

def fisher_z(r):
    return 0.5 * np.log((1 + r) / (1 - r))

def __statistical_analysis(x, y, mode='exp'):
    s = None
    p = None
    if mode == 'exp':
        # Normality test
        if scipy.stats.shapiro(x).pvalue > 0.05 and scipy.stats.shapiro(y).pvalue > 0.05:
            # x and y are normally distributed: Welch's t-test
            _, p = scipy.stats.ttest_ind(x, y, equal_var=False)
        else:
            # otherwise: Mann–Whitney
            _, p = scipy.stats.mannwhitneyu(x, y, alternative='two-sided')

    elif mode == 'co-exp':
        s, p = scipy.stats.pearsonr(x=x, y=y, alternative='two-sided')

    return s, p

def _normalize_length(x, y):
    if len(x) > len(y):
        x = np.random.choice(x, size=len(y), replace=False)
    elif len(y) > len(x):
        y = np.random.choice(y, size=len(x), replace=False)
    return x, y

def inferring_intralayer_edges(
    gene_symbol_ids_of_interest: list,
    gene_exp_by_tissue: pd.DataFrame,
    tissues_of_interest: list,
    multidimensional_network: pymnet.MultilayerNetwork = None,
    p_threshold: float = 1e-8,
    s_threshold: list = [0.5, 0.9],
    edge_colors: dict = {},
    verbose: bool = True
) -> tuple[pymnet.MultilayerNetwork, dict]:
    
    if verbose == True:
        print("-> Inferring intra-layer edges...")

    if multidimensional_network is None:
        multidimensional_network = pymnet.MultilayerNetwork(aspects=1, directed=False, fullyInterconnected=True)

    gene_pairs = list(combinations(gene_symbol_ids_of_interest, 2))

    # Collecting results for global correction
    results = []

    for tissue in tissues_of_interest:

        multidimensional_network.add_layer(tissue)

        if verbose == True:
            print("-- Layer " + str(tissue) + ") ...")

        for u, v in gene_pairs:

            u_data = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == u) & (gene_exp_by_tissue['tissueSiteDetailId'] == tissue) ]
            v_data = gene_exp_by_tissue.loc[ (gene_exp_by_tissue['geneSymbol'] == v) & (gene_exp_by_tissue['tissueSiteDetailId'] == tissue) ]
            
            if u_data.empty or v_data.empty:
                continue
            
            u_gexp = u_data.data.values[0]
            v_gexp = v_data.data.values[0]

            if len(u_gexp) < 3 or len(v_gexp) < 3:
                #if verbose == True:
                    #warnings.warn("Size of samples are empty or less than 3; skipped.")
                continue

            s, p = __statistical_analysis(u_gexp, v_gexp, mode='co-exp')

            results.append(((u, v, tissue), s, p))

    if not results:
        return multidimensional_network, edge_colors
    
    # Apply global FDR correction
    keys, scores, pvals = zip(*results)
    _, padj_all, _, _ = multipletests(pvals, method='fdr_bh')

    # Filter and add edges
    for (u, v, current_layer), s, padj in zip(keys, scores, padj_all):
        if padj < p_threshold and (s_threshold[0] < s < s_threshold[1]):
            multidimensional_network[u, v, current_layer] = padj
            edge_colors[((u, current_layer), (v, current_layer))] = 'blue'

    if verbose == True:
        print("- [OK]")

    return multidimensional_network, edge_colors


def inferring_interlayer_edges(
    gene_symbol_ids_of_interest: list,
    gene_exp_by_tissue: pd.DataFrame,
    tissues_of_interest: list,
    multidimensional_network: pymnet.MultilayerNetwork = None,
    p_threshold: float = 1e-8,
    s_threshold: list = [0.5, 0.9],
    mode: str = "equal_exp",
    edge_colors: dict = {},
    verbose: bool = True
) -> tuple[pymnet.MultilayerNetwork, dict]:
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

    results = []

    if mode in ['equal_exp', 'diff_exp']:
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
                results.append(((gene, l1, l2), None, p))

        if not results:
            return multidimensional_network, edge_colors
        
        keys, _, pvals = zip(*results)
        _, padj_all, _, _ = multipletests(pvals, method='fdr_bh')

        # Applying global multiple testing correction
        for (gene, l1, l2), adj_p in zip(keys, padj_all):
            if (mode == 'equal_exp' and adj_p > p_threshold) or (mode == 'diff_exp' and adj_p < p_threshold):
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
                
                if u_data_l1.empty or v_data_l1.empty or u_data_l2.empty or v_data_l2.empty:
                    continue

                u_gexp_l1, v_gexp_l1 = u_data_l1.data.values[0], v_data_l1.data.values[0]
                u_gexp_l2, v_gexp_l2 = u_data_l2.data.values[0], v_data_l2.data.values[0]
                
                if len(u_gexp_l1) < 3 or len(v_gexp_l1) < 3 or len(u_gexp_l2) < 3 or len(v_gexp_l2) < 3:
                    continue

                x1, y1 = _normalize_length(u_gexp_l1, v_gexp_l1)
                x2, y2 = _normalize_length(u_gexp_l2, v_gexp_l2)

                s1, _ = __statistical_analysis(x1, y1, mode='co-exp')
                s2, _ = __statistical_analysis(x2, y2, mode='co-exp')

                # Transform
                z1, z2 = fisher_z(s1), fisher_z(s2)

                # Number of samples per tissue
                n1 = len(u_gexp_l1)
                n2 = len(u_gexp_l2)

                # Standard error
                se = np.sqrt(1/(n1 - 3) + 1/(n2 - 3))

                # z-score
                z_diff = (z1 - z2) / se

                # p-value
                p = 2 * (1 - scipy.stats.norm.cdf(abs(z_diff)))

                results.append(((u, v, l1, l2), None, p))

        if not results:
            return multidimensional_network, edge_colors
        
        # Applying global multiple testing correction
        keys, _, pvals = zip(*results)
        _, padj_all, _, _ = multipletests(pvals, method='fdr_bh')

        for (u, v, l1, l2), adj_p in zip(keys, padj_all):
            multidimensional_network[u, v, l1, l2] = adj_p
            edge_colors[((u, l1), (v, l2))] = 'green'
    
    return multidimensional_network, edge_colors