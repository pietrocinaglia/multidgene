'''
' NOTE: 'pip install pyarrow' for parquet support
'''

import os
import pandas as pd
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

###
# @project: multi-DGeNe (PRE-PROCESSING: DATASET CONSTRUCTION)
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues. This is used for preprocessing in order to build a (whole) dataset from the fragmented GTEx's data.
# @url: https://github.com/pietrocinaglia/multidgene
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"

FILE_LIST = WORKSPACE + 'tissues.csv'
DATA_DIR = WORKSPACE + 'gene_read_counts_by_tissue/' # put here all .gct.gz for the read counts data by tissue (see "Gene read counts by tissue" at https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)
OUTPUT_FILE_PATH = WORKSPACE + 'gene_reads_v10_by_tissue.parquet'
COMBINING_OUTPUT = False

tissue_names = pd.read_csv(FILE_LIST)

print("Loading files:")
data = []
for tissue in tissue_names.iloc[:, 1]:
    file_path = os.path.join(DATA_DIR, "gene_reads_v10_"+tissue.lower()+".gct.gz")
    print("\t - " + tissue)
    try:
        df = pd.read_csv(file_path, delimiter='\t', header=2, compression='infer')
        df['tissue'] = tissue

        data.extend(df.to_dict(orient='records')) # duct optimize memory usage in combining pd.DataFrame

        print("\t -- [OK]" )
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

print("Combining data...")
gene_reads_v10_by_tissue = pd.DataFrame(data)
data = None

# Preview
print(gene_reads_v10_by_tissue.head(2))

print("Storing data as file (compressed .parquet)...")
gene_reads_v10_by_tissue.to_parquet(OUTPUT_FILE_PATH, index=False, compression='zstd') # , partition_cols=['tissue']) # compression='snappy'