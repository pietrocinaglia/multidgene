import os
import pandas as pd

from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

###
# @project: multi-DGeNe (PRE-PROCESSING: DATASET CONSTRUCTION)
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference based on Weighted Gene Co-expression Network Analysis (WGCNA) across human tissues. This is used for preprocessing in order to build a (whole) dataset from the fragmented GTEx's data.
# @url: https://github.com/pietrocinaglia/multidgene
###

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
OUTPUT_FILE_PATH = WORKSPACE + 'gene_reads_v10_by_tissue.parquet.zstd'

FILE_LIST = WORKSPACE + 'tissues.csv'
DATA_DIR = WORKSPACE + 'gene_read_counts_by_tissue/' # put here all .gct.gz for the read counts data by tissue (see "Gene read counts by tissue" at https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression)

all_data = []
all_columns = set()

print("Collecting tissue identifiers...")
file_names = pd.read_csv(FILE_LIST)
#
for file in file_names.iloc[:, 1]:  # tissue names are in the second column
    file_path = os.path.join(DATA_DIR, "gene_reads_v10_"+file.lower()+".gct.gz")
    try:
        df = pd.read_csv(file_path, delimiter='\t', header=2, compression='infer')
        all_columns.update(df.columns)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

print("Loading files:")
for file in file_names.iloc[:, 1]:
    file_path = os.path.join(DATA_DIR, "gene_reads_v10_"+file.lower()+".gct.gz")
    print("\t - " + file)
    try:
        df = pd.read_csv(file_path, delimiter='\t', header=2, compression='infer')
        missing_cols = all_columns - set(df.columns)
        for col in missing_cols:
            df[col] = 0
        # Reorder columns consistently
        df = df[list(all_columns)] 
        # Add a column with the tissue name
        df['source_file'] = file  
        all_data.append(df)
        print("\t -- [OK]" )
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

print("Data concatenation...")
combined_df = pd.concat(all_data, ignore_index=True)

print("Storing data as file (.parquet with ZSTD compression)...")
combined_df.to_parquet(OUTPUT_FILE_PATH, index=False, compression='zstd')

# Preview
print(combined_df.head(2))