'''
' NOTE: 'pip install pyarrow' for parquet support
'''

import os

import pandas as pd
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

import numpy as np
import pickle as pk

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
OUTPUT_FILE_PATH = WORKSPACE + 'gene_reads_v10_by_tissue'
COMBINING_OUTPUT = False

print("Loading tissue identifiers...")
tissue_names = pd.read_csv(FILE_LIST)

print("Collecting samples identifiers...")
sample_ids = set()
for tissue in tissue_names.iloc[:, 1]:  # tissue names are in the second column
    file_path = os.path.join(DATA_DIR, "gene_reads_v10_"+tissue.lower()+".gct.gz")
    try:
        df = pd.read_csv(file_path, delimiter='\t', header=2, nrows=1, compression='infer')
        sample_ids.update(df.columns)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

print("Defining resulting data model...")
df_model = pd.DataFrame(columns=['tissue']+list(sample_ids))

print("Loading files:")
data = []
for tissue in tissue_names.iloc[:, 1]: # names are located in (index) 1
    file_path = os.path.join(DATA_DIR, "gene_reads_v10_"+tissue.lower()+".gct.gz")
    print("\t - " + tissue)
    try:
        tmp = pd.read_csv(file_path, delimiter='\t', header=2, compression='infer')
        tmp = pd.concat([df_model, tmp]) # port to reference model
        tmp['tissue'] = tissue
        #data.extend(df.to_dict(orient='records')) # dict optimize memory usage in combining pd.DataFrame
        data.append(df.to_numpy())
        print("\t -- [OK]" )
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

print("Concatenating data as file...")
output = np.concatenate(data)

print("Storing data as file...")
with open(OUTPUT_FILE_PATH+'.pkl', 'wb') as fp:
    pk.dump(output, fp)

'''
print("Loading data from file...")
data = None
with open(OUTPUT_FILE_PATH+'.pkl', 'rb') as fp:
    data = pk.load(fp)
    
print("Converting to dataframe...")
gene_reads_v10_by_tissue = pd.DataFrame(data)
data = None

print("Storing data as file...")
gene_reads_v10_by_tissue.to_parquet(OUTPUT_FILE_PATH+'.parquet', index=False, compression='zstd') # , partition_cols=['tissue']) # compression='snappy'
'''

print("[ DONE ]")