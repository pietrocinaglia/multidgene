# -*- coding: utf-8 -*-
import os
import pandas as pd
import requests


class Reader:
    DATASET_ID = 'gtex_v10'
    LIMIT = 100000

    def __init__(self, data_path: str = None, verbose: bool = True):
        self.data_path = data_path or os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/')
        self.verbose = verbose

    def load_metadata_from_file(self, genes_of_interest: list):
        if self.verbose:
            print("- Loading Metadata from file...")

        metadata_path = os.path.join(self.data_path, 'metadata.csv')
        metadata = pd.read_csv(metadata_path, header=0, sep=" ")
        metadata = metadata[metadata['gene_symbol'].isin(genes_of_interest)]

        symbol_to_id = dict(zip(metadata.gene_symbol, metadata.ensembl_id))
        return list(symbol_to_id.keys()), list(symbol_to_id.values())

    def load_gene_expression_read_counts_from_file(self, gene_code_ids_of_interest: list, tissues_of_interest: list = None) -> pd.DataFrame:
        if self.verbose:
            print("- Retrieving Gene Expression Read Counts by tissue (from file)...")

        if tissues_of_interest is None:
            tissue_file = os.path.join(self.data_path, 'tissues.csv')
            tissues_df = pd.read_csv(tissue_file)
            tissues_of_interest = list(tissues_df.iloc[:, 1])

        records = []
        for tissue in tissues_of_interest:
            file_path = os.path.join(self.data_path, f"gene_read_counts_by_tissue/gene_reads_v10_{tissue.lower()}.gct.gz")
            if self.verbose:
                print(f"\t - {tissue}")
            try:
                df = pd.read_csv(file_path, delimiter='\t', header=2, compression='infer')
                df = df[df['Name'].isin(gene_code_ids_of_interest)]
                df['tissueSiteDetailId'] = tissue
                records.append(df)
                if self.verbose:
                    print("\t -- [OK]")
            except Exception as e:
                print(f"\t -- Error reading {file_path}: {e}")

        if not records:
            raise RuntimeError("No gene expression data was loaded. Check gene IDs or tissue names.")

        result = pd.concat(records, ignore_index=True)
        result.rename(columns={'Description': 'geneSymbol'}, inplace=True)
        result.drop(columns=['Name'], inplace=True)

        return result

    def retrieve_gene_expression_tpm_from_api(self, gene_code_ids_of_interest: list, tissues_of_interest: list = None, page_i: int = 0) -> tuple[pd.DataFrame, list]:
        if self.verbose:
            print("- Retrieving Gene Expression TPM by tissue (from API)...")

        url = "https://gtexportal.org/api/v2/expression/geneExpression"
        params = {
            'datasetId': self.DATASET_ID,
            'itemsPerPage': self.LIMIT,
            'page': page_i,
            'gencodeId': gene_code_ids_of_interest,
            'tissueSiteDetailId': tissues_of_interest or [],
            'format': 'json'
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            result = response.json()
        except Exception as e:
            raise RuntimeError(f"API request failed: {e}")

        if not result or 'data' not in result:
            raise ValueError("Invalid response from GTEx API: 'data' field not found.")

        data = pd.DataFrame(result['data'])

        if tissues_of_interest is None or len(tissues_of_interest) == 0:
            tissues_of_interest = list(data["tissueSiteDetailId"].unique())

        data = data[data['unit'] == 'TPM']

        return data, tissues_of_interest
