import os
import re
import sys
import pandas as pd
import argparse
import numpy as np

maf = pd.read_csv(snakemake.input["tsv"], sep="\t")
# vcf_path = '/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/test'
# csv_path = '/Users/asherpreskasteinberg/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/lab_notebook/APS017.1_3x3_SV_analysis/test'
# maf.to_csv(csv_path, sep='\t', index=False)
# split into many mafs for gene annotation step
os.makedirs(snakemake.params["split_out"], exist_ok=True)
out_paths = snakemake.output["split_tsv"]
num_splits = len(out_paths)
# Split the DataFrame
split_dfs = np.array_split(maf, num_splits)
# Convert split arrays back into DataFrames
split_dfs = [pd.DataFrame(split) for split in split_dfs]
for i in np.arange(0, len(split_dfs)):
    split_df = split_dfs[i]
    split_path = out_paths[i]
    split_df.to_csv(split_path, sep='\t', index=False)
