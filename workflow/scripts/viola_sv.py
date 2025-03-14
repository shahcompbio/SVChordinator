import viola
import gzip
import shutil
import os

"""
make variant tables for illumina SVs using viola-sv
"""

# functions

def decompress_vcf(compressed_vcf, decompressed_vcf):
    # Decompress vcf and keep original file
    with gzip.open(compressed_vcf, 'rb') as f_in:
        with open(decompressed_vcf, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return

# inputs

vcf = snakemake.input["vcf"]
temp_vcf = snakemake.params["vcf"]
output_table = snakemake.output["tsv"]

# check if the vcf is compressed
if vcf.endswith(".gz"):
    # path to temp vcf
    decompress_vcf(vcf, temp_vcf)
else:
    os.system(f"cp {vcf} {temp_vcf}")

# run viola-sv; gripss format works for gripss/manta/svaba
vcf = viola.read_vcf(temp_vcf,
                     variant_caller="gripss",
                     patient_name="sample")
df = vcf.get_table("positions")
cols = ["chrom1", "pos1", "id", "svtype"]
smalldf = df[cols].copy()
# Combine `strand1` and `strand2` into a new column
smalldf['strands'] = df['strand1'] + df['strand2']
# blank columns to conform to annotation rule
smalldf["STRAND"] = "."
smalldf["BP_notation"] = "."
smalldf.to_csv(output_table, sep="\t", index=None)

