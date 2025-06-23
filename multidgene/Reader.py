# -*- coding: utf-8 -*-

import os
import pandas as pd
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
    verbose = True

    def __init__(self, data_path:str=None, verbose:bool=True):

        # default data_path
        if data_path is None: 
            self.data_path = os.path.dirname(os.path.realpath(__file__)) + '/data/'
        else:
            self.data_path = data_path

        self.verbose = verbose

        pass
    
    def load_metadata_from_file(self, genes_of_interest):

        if self.verbose == True:
            print("- Loading Metadata from file...")

        source_path = self.data_path + 'metadata.csv'

        metadata = pd.read_csv(source_path, header=0, sep = " ")
        metadata = metadata[ metadata['gene_symbol'].isin(genes_of_interest) ]
        metadata = dict(zip(metadata.gene_symbol,metadata.ensembl_id))
        
        gene_symbol_ids_of_interest = list(metadata.keys())
        gene_code_ids_of_interest = list(metadata.values())

        return gene_symbol_ids_of_interest, gene_code_ids_of_interest
    
    def load_gene_expression_read_counts_from_file(self, gene_code_ids_of_interest:list, tissues_of_interest:list=None) -> pd.DataFrame:

        if self.verbose == True:
            print("- Retrieving Gene Expression Read Counts by tissue (from file)...")

        if tissues_of_interest is None:
            tissues_of_interest = pd.read_csv(self.data_path+'tissues.csv')
            tissues_of_interest = list(tissue_names.iloc[:, 1])

        data = []
        for tissue in tissues_of_interest:
            file_path = os.path.join(self.data_path, "gene_read_counts_by_tissue/gene_reads_v10_"+tissue.lower()+".gct.gz")
            print("\t - " + tissue)
            try:
                tmp = pd.read_csv(file_path, delimiter='\t', header=2, compression='infer')
                tmp = tmp.loc[ tmp['Name'].isin(gene_code_ids_of_interest) ] # Name is Gene-ID, Description is Gene-Name
                tmp['tissueSiteDetailId'] = tissue
                data.extend(tmp.to_dict(orient='records'))
                tmp = None
                print("\t -- [OK]" )
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

        gene_expression_read_counts_by_tissue = pd.DataFrame(data)
        data = None
        gene_expression_read_counts_by_tissue.rename({'Description': 'geneSymbol'}, axis=1, inplace=True)
        gene_expression_read_counts_by_tissue.drop("Name", axis='columns', inplace=True) # unused identifier
        #gene_expression_read_counts_by_tissue.set_index(['geneSymbol', 'tissueSiteDetailId'], inplace=True) # Name is Gene-ID, Description is Gene-Name

        return gene_expression_read_counts_by_tissue


    def retrieve_gene_expression_tpm_from_api(self, gene_code_ids_of_interest:list, tissues_of_interest:list=[], page_i=0) -> pd.DataFrame:
        
        if self.verbose == True:
            print("- Retrieving Gene Expression TPM by tissue (from API)...")

        api = "https://gtexportal.org/api/v2/expression/geneExpression"
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