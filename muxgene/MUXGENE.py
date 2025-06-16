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
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multiplex Gene Network (muxGeNe) inference from gene expression across human tissues.
# @url: https://github.com/pietrocinaglia/muxGeNe
###

class MUXGENE:
    WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
    #
    UNIT = 'TPM'
    DATASETID = "gtex_v10"
    LIMIT = 100000 # Max is 100000
    #
    gene_symbol_ids_of_interest = list()
    gene_code_ids_of_interest = list()
    tissues_of_interest = list()
    p_threshold = 10e-8
    s_threshold = [0.5,0.9]
    verbose = False
    #

    def __init__(self, genes_of_interest:list, tissues_of_interest:list=[], p_threshold:float=0, verbose:bool=False):
        if len(genes_of_interest) == 0:
            raise Exception("Genes of interest ('genes_of_interest:list') are mandatory; this list cannot be empty or none.")

        if self.verbose == True:
            print( "***********************************************************" )
            print( " ------------------------------------" )
            print( " > Multiplex Gene Network (muxGeNe) < " )
            print( " ------------------------------------" )
            print( "  Version:\t0.1" )
            print( "  GitHub:\t" + "https://github.com/pietrocinaglia/muxgene" )
            print( "  Contact:\tPietro Cinaglia (cinaglia@unicz.it)" )
            print( "***********************************************************" )
            print()
        
        print("[INFO] Loading metadata...")
        
        self.genes_of_interest = genes_of_interest
        self.tissues_of_interest = tissues_of_interest
        if p_threshold > 0:
            self.p_threshold = p_threshold
        self.verbose = verbose

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
  
    def __retrieve_gtex_data(self) -> pd.DataFrame:
        result = None

        api = "https://gtexportal.org/api/v2/expression/geneExpression"
        #api = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        params = {'itemsPerPage': self.LIMIT, 'datasetId': self.DATASETID, 'gencodeId': self.gene_code_ids_of_interest, 'tissueSiteDetailId': self.tissues_of_interest, 'format':'json'}
        result =  requests.get( api, params ).json()
    
        # Checking result
        if result is None:
            raise Exception("An error occurring in data retrieving.")
        
        data = pd.DataFrame(result['data'])

        if len(self.tissues_of_interest) == 0:
            self.tissues_of_interest = list(data["tissueSiteDetailId"].unique())

        return data
    
    def construct_multiplex_network(self) -> list:
        if self.verbose == True:
            print("[INFO] Data Retrieving...")       

        df = self.__retrieve_gtex_data()

        print(df) ###################

        gene_pairs = list(combinations(self.gene_symbol_ids_of_interest, 2))

        if self.verbose == True:
            print("- Data concerning gene expression across human tissues [OK] ")
            print()
            print("[INFO] Multiplex Network Construction...")
        
        # Only for testing (1/2)
        with open(self.WORKSPACE+'statistics_result.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            field = ['gene1', 'gene2', 'tissue', 's', 'p', 'is_significative']
            writer.writerow(field)

        multiplex_network = pymnet.MultiplexNetwork()

        for i in range(len(self.tissues_of_interest)):

            current_tissue = self.tissues_of_interest[i]

            multiplex_network.add_layer(current_tissue)

            if self.verbose == True:
                print("- Layer " + str(i+1) + " (" + str(current_tissue) + ") ...")

            ### Work in progress
            for u, v in gene_pairs:
                
                u_data = df.loc[ (df['geneSymbol'] == u) & (df['unit'] == self.UNIT) & (df['tissueSiteDetailId'] == current_tissue) ]
                v_data = df.loc[ (df['geneSymbol'] == v) & (df['unit'] == self.UNIT) & (df['tissueSiteDetailId'] == current_tissue) ]
                
                if (u_data.empty or v_data.empty):
                    continue
                
                u_gexp = u_data.data.values[0]
                v_gexp = v_data.data.values[0]
    
                if ( len(u_gexp) < 3 or len(v_gexp) < 3):
                    if self.verbose == True:
                        raise Warning("Size of sample less than 3; skipped.")
                    continue

                s, p = scipy.stats.pearsonr(u_gexp, v_gexp)
                
                p = round(p, 8)
                if p < self.p_threshold and (s > self.s_threshold[0] and s < self.s_threshold[1]):
                    #current_layer.add_edge(u,v, pvalue=p)
                    multiplex_network[u,v,current_tissue,current_tissue] = self.p_threshold
                
                # Only for testing (2/2)
                with open(self.WORKSPACE+'statistics_result.csv', 'a', newline='') as file:
                    writer = csv.writer(file)
                    # if p < self.p_threshold and (s > self.s_threshold[0] and s < self.s_threshold[1]):
                    writer.writerow([u, v, current_tissue, s, p, ('yes' if p < self.p_threshold and (s > self.s_threshold[0] and s < self.s_threshold[1]) else 'no')])

            if self.verbose == True:
                print("-- [OK]")
        
        if self.verbose == True:
            print("- Multiplex Network built.")

        return multiplex_network

    def plot(self, multiplex_network:pymnet.MultiplexNetwork, output_path:str=None):
        
        fig = pymnet.draw(
            multiplex_network,
            autoscale = True,
            layershape='rectangle',
            layerPadding=0.2,
            #defaultLayerColor = "#DDDDDD",
            defaultLayerLabelLoc = (1.1,1.1),
            defaultLayerLabelColor = 'blue',
            defaultNodeLabelColor = 'black',
            defaultNodeColor = 'black',
            nodeSizeRule = {'rule': 'scaled', 'scalecoeff': 0.1},
            defaultEdgeColor = 'green',
            defaultNodeLabelSize = 6,
            alignedNodes = True,
            figsize = (7,4)
            )
        
        if output_path is not None:
            filepath = output_path + "multiplex_network"
            fig.savefig(filepath + ".pdf")
            #pymnet.webplot(multiplex_network, outputfile=filepath + ".html")
            print("- Plot saved as file [OK]")

        plt.show()

    def save_as_files(self, multiplex_network:pymnet.MultiplexNetwork, output_path:str):
        
        if output_path is None:
            raise Exception("You must provide an output path for storing data.")
        
        filepath = output_path + "multiplex_network.txt"

        pymnet.write_edge_file(multiplex_network, filepath, sep='\t', weights=True, numericNodes=False)
        print("- Multiplex Network saved as file [OK]")