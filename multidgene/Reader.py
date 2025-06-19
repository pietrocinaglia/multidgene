# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import requests

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

class GTExDataReader:

    # GTEx API v2
    DATASET_ID = 'gtex_v10' # 'gtex_v8' 'gtex_v10' 'gtex_snrnaseq_pilot'
    LIMIT = 100000 # Max is 100000
    #
    data_path = None

    def __init__(self, data_path:str=None):

        # default data_path
        if data_path is None: 
            self.data_path = os.path.dirname(os.path.realpath(__file__)) + '/data/'
        else:
            self.data_path = data_path

        pass
    
    def load_metadata_from_file(self, genes_of_interest):

        source_path = self.data_path + 'metadata.csv'

        metadata = pd.read_csv(source_path, header=0, sep = " ")
        metadata = metadata[ metadata['gene_symbol'].isin(genes_of_interest) ]
        metadata = dict(zip(metadata.gene_symbol,metadata.ensembl_id))
        
        gene_symbol_ids_of_interest = list(metadata.keys())
        gene_code_ids_of_interest = list(metadata.values())

        return gene_symbol_ids_of_interest, gene_code_ids_of_interest
    
    def load_gene_expression_read_counts_from_file(self) -> pd.DataFrame:

        source_path = self.data_path + 'GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz'
        # index_col => [ 0 is Ensembl-Id, 1 is Gene Name ]
        df = pd.read_csv(source_path, delimiter='\t', header=2, index_col=1, compression='infer')  

        return df

    def retrieve_gene_expression_tpm_from_api(self, gene_code_ids_of_interest:list, tissues_of_interest:list=[], page_i=0) -> pd.DataFrame:
        
        api = "https://gtexportal.org/api/v2/expression/geneExpression"
        #api = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        params = {'datasetId': self.DATASET_ID, 'itemsPerPage': self.LIMIT, 'page': page_i, 'gencodeId': gene_code_ids_of_interest, 'tissueSiteDetailId': tissues_of_interest, 'format':'json'}
        result =  requests.get( api, params ).json()
    
        # Checking result
        if result is None or 'data' not in result:
            raise Exception("An error occurring in data retrieving.")
        
        data = pd.DataFrame(result['data'])

        if len(tissues_of_interest) == 0:
            tissues_of_interest = list(data["tissueSiteDetailId"].unique())

        # GTEx API v2 only provides data in TPM (as unit), we only used this filter as "double-check" in case of unexpected API updates.
        data = data.loc[(data['unit'] == 'TPM')]

        return data, tissues_of_interest