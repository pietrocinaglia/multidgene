# -*- coding: utf-8 -*-

import os
import requests
import pymnet
import pandas as pd
import scipy
import statistics
from itertools import combinations
import matplotlib.pyplot as plt
import csv

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

class MultiDGene:
    WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
    #
    UNIT = 'TPM' # 'TPM' 'read count' 'log(CP10K)'
    DATASET_ID = 'gtex_v10' # 'gtex_v8' 'gtex_v10'
    LIMIT = 100000 # Max is 100000
    #
    gene_symbol_ids_of_interest = list()
    gene_code_ids_of_interest = list()
    tissues_of_interest = list()
    p_threshold = 10e-8
    s_threshold = [0.5,0.9] # lower and upper limit
    verbose = False
    #

    def __init__(self, genes_of_interest:list, tissues_of_interest:list=[], p_threshold:float=0, verbose:bool=False):
        if len(genes_of_interest) == 0:
            raise Exception("Genes of interest ('genes_of_interest:list') are mandatory; this list cannot be empty or none.")

        if self.verbose == True:
            print( "***************************************************************" )
            print( " --------------------------------------------------------------" )
            print( " >        multi-Dimensional Gene Network (multi-DGeNe)        <" )
            print( " --------------------------------------------------------------" )
            print( "  Description:\tMulti-Dimensional Gene Network (multi-DGeNe)   " )
            print( "              \tinference based on Weighted Gene Co-expression " )
            print( "              \t Network Analysis (WGCNA) across human tissues." )
            print( " --------------------------------------------------------------" )
            print( "  Version:\t0.1                                                " )
            print( "  GitHub:\t" + "https://github.com/pietrocinaglia/multidgene   " )
            print( "  Contact:\tPietro Cinaglia (cinaglia@unicz.it)                " )
            print( " --------------------------------------------------------------" )
            print( "***************************************************************" )
            print()
        
        print("[INFO] Loading metadata...")
        
        self.verbose = verbose

        self.genes_of_interest = genes_of_interest
        self.tissues_of_interest = tissues_of_interest

        if p_threshold > 0:
            self.p_threshold = p_threshold

        # reading metadata
        df = pd.read_csv(self.WORKSPACE+'metadata.csv', header=0, sep = " ")

        # filtering metadata
        df = df[ df['gene_symbol'].isin(self.genes_of_interest) ]

        # dataframe to dictionary
        metadata = dict(zip(df.gene_symbol,df.ensembl_id))
        self.gene_symbol_ids_of_interest = list(metadata.keys())
        self.gene_code_ids_of_interest = list(metadata.values())

        print("- Metadata loaded [OK] ")

        pass
  
    def __retrieve_gtex_data(self, page_i=0) -> pd.DataFrame:
        
        api = "https://gtexportal.org/api/v2/expression/geneExpression"
        #api = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        params = {'datasetId': self.DATASET_ID, 'itemsPerPage': self.LIMIT, 'page': page_i, 'gencodeId': self.gene_code_ids_of_interest, 'tissueSiteDetailId': self.tissues_of_interest, 'format':'json'}
        result =  requests.get( api, params ).json()
    
        # Checking result
        if result is None or 'data' not in result:
            raise Exception("An error occurring in data retrieving.")
        
        data = pd.DataFrame(result['data'])

        if len(self.tissues_of_interest) == 0:
            self.tissues_of_interest = list(data["tissueSiteDetailId"].unique())

        #data = data.loc[(data['unit'] == self.UNIT)]

        print(data)

        return data
    
    def construct_multidimensional_network(self) -> tuple[pymnet.MultilayerNetwork, set]:
        if self.verbose == True:
            print("[INFO] Data Retrieving...")       

        df = self.__retrieve_gtex_data()

        gene_pairs = list(combinations(self.gene_symbol_ids_of_interest, 2))

        multilayer_network = pymnet.MultilayerNetwork(aspects=1, directed=False, fullyInterconnected=True)
        edge_colors = {}

        if self.verbose == True:
            print("- [OK] ")
            print()
            print("[INFO] Multi-Dimensional Network Construction...")
        
        multilayer_network, edge_colors = self.__inferring_intralayer_edges(df, gene_pairs, multilayer_network, edge_colors)

        # layer_pairs = list(combinations(self.tissues_of_interest, 2))
        multilayer_network, edge_colors = self.__inferring_interlayer_edges(df, self.gene_symbol_ids_of_interest, self.tissues_of_interest, multilayer_network, edge_colors)

        '''
        # Only for testing (1/2)
        with open(self.WORKSPACE+'statistics_result.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            field = ['gene1', 'gene2', 'tissue', 's', 'p', 'is_significative']
            writer.writerow(field)
        # Only for testing (2/2)
        with open(self.WORKSPACE+'statistics_result.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            # if p < self.p_threshold and (s > self.s_threshold[0] and s < self.s_threshold[1]):
            writer.writerow([u, v, current_tissue, s, p, ('yes' if p < self.p_threshold and (s > self.s_threshold[0] and s < self.s_threshold[1]) else 'no')])
        '''

        if self.verbose == True:
            print("- Multilayer Network modelled.")
        
        return multilayer_network, edge_colors

    def __inferring_intralayer_edges(self, df:pd.DataFrame, gene_pairs:list, multilayer_network:pymnet.MultilayerNetwork, edge_colors:set) -> tuple[pymnet.MultilayerNetwork, set]:

        if self.verbose == True:
            print("- Inferring intra-layer edges...")

        for current_layer in self.tissues_of_interest:

            multilayer_network.add_layer(current_layer)

            if self.verbose == True:
                print("-- Layer " + str(current_layer) + ") ...")

            for u, v in gene_pairs:

                u_data = df.loc[ (df['geneSymbol'] == u) & (df['tissueSiteDetailId'] == current_layer) ]
                v_data = df.loc[ (df['geneSymbol'] == v) & (df['tissueSiteDetailId'] == current_layer) ]
                
                if (u_data.empty or v_data.empty):
                    continue
                
                u_gexp = u_data.data.values[0]
                v_gexp = v_data.data.values[0]
    
                if ( len(u_gexp) < 3 or len(v_gexp) < 3):
                    if self.verbose == True:
                        raise Warning("Size of sample less than 3; skipped.")
                    continue

                s, p = scipy.stats.pearsonr(x=u_gexp, y=v_gexp, alternative='two-sided')

                if p < self.p_threshold and (s > self.s_threshold[0] and s < self.s_threshold[1]):
                    multilayer_network[u,v,current_layer] = self.p_threshold
                    edge_colors[((u, current_layer), (v, current_layer))] = 'blue'

            if self.verbose == True:
                print("-- [OK]")

        return multilayer_network, edge_colors


    def __inferring_interlayer_edges(self, df:pd.DataFrame, genes:list, layers:list, multilayer_network:pymnet.MultilayerNetwork, edge_colors:set) -> tuple[pymnet.MultilayerNetwork, set]:
        
        if self.verbose == True:
            print("- Inferring inter-layer edges...")

        for i in range(len(layers)):
            j = i + 1
            if j > (len(layers)-1):
                break

            for g in genes:
                multilayer_network[g,g,layers[i],layers[j]] = self.p_threshold
                edge_colors[((g, layers[i]), (g, layers[j]))] = 'green'

        if self.verbose == True:
            print("--- [OK]")

        return multilayer_network, edge_colors
    
    def plot(self, multilayer_network:pymnet.MultilayerNetwork, edge_colors:set={}, output_path:str=None):
        
        fig = pymnet.draw(
            multilayer_network,
            layershape='rectangle',
            layerPadding=0.2,
            defaultLayerColor = "#BDBBBB",
            defaultLayerLabelLoc = (1.1,1.1),
            defaultLayerLabelColor = "#6E6E6E",
            defaultNodeLabelColor = "#000000",
            defaultNodeColor = "#000000",
            nodeSizeRule = {'rule': 'scaled', 'scalecoeff': 0.2},
            defaultEdgeColor = 'black',
            edgeColorDict = edge_colors,
            defaultNodeLabelSize = 6,
            alignedNodes = True,
            figsize = (7,4),
            autoscale = True,
            )
        
        if output_path is not None:
            filepath = output_path + "multilayer_network"
            fig.savefig(filepath + ".pdf")
            #pymnet.webplot(multilayer_network, outputfile=filepath + ".html")
            print("- Plot saved as file [OK]")

        plt.show()

    def save_as_files(self, multilayer_network:pymnet.MultilayerNetwork, output_path:str):
        
        if output_path is None:
            raise Exception("You must provide an output path for storing data.")
        
        filepath = output_path + "multilayer_network.txt"

        pymnet.write_edge_file(multilayer_network, filepath, sep='\t', weights=True, numericNodes=False)
        print("- multilayer Network saved as file [OK]")