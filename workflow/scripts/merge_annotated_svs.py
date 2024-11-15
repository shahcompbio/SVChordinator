import pandas as pd
import numpy as np

split_tsvs = snakemake.input["split_tsvs"]
out_tsv = snakemake.output["all_SVs"]
#######
all_SVs = pd.DataFrame()

for split_tsv in split_tsvs:
    SVs = pd.read_csv(split_tsv, sep="\t")
    all_SVs = pd.concat([all_SVs, SVs])

all_SVs.to_csv(out_tsv, sep='\t', index=False)
