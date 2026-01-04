import os
import numpy as np
import pandas as pd

### paths ###
out_dir = config["out_dir"]
sample_metadata = config["samples"]
#sample_name = config["sample_name"]
oncokb = config["annotate"]["oncokb"]
gene_annotations = config["annotate"]["gene_annotation"]
sv_type = config["sv_type"]
ideo = config["reference"]
min_callers = config["min_callers"]
### helper functions ###
def extract_sample_list(sample_metadata):
    """
    extract sample names from sample metadata file
    :return: list of sample names
    """
    df = pd.read_csv(sample_metadata, sep="\t")
    sample_names = list(df["sample"].unique())
    return sample_names

def extract_individual_calls(sample_metadata):
    """
    extract a dict of info for each individual set of SV calls
    :param sample_metadata: metadata sheet
    :return: pandas dataframe of vcf paths, names of individual callers
    """
    df = pd.read_csv(sample_metadata, sep="\t")
    callers = list(df["caller"])
    return callers

def define_caller_table_targets(sample_metadata, sample=None):
    """
    define targets for raw SV calls for ONT and ILL
    :param sample_metadata: metadata sheet
    :return: list of target paths for raw SV calls
    """
    df = pd.read_csv(sample_metadata, sep="\t")
    print(sample)
    if sample is not None:
        df = df[df["sample"] == sample]
    for row in df.itertuples():
        tech = row.tech
        caller = row.caller
        sample = row.sample
        assert tech in ["ONT", "ILL"], f"technology {tech} not supported for {sample}/{caller}"
        target = os.path.join(out_dir, "raw_SVs", 
        sample, tech, sample + f".{caller}.tsv")
        targets.append(target)
    return targets
# extract list of samples
samples = extract_sample_list(sample_metadata)
print(samples)
def get_output():
    """
    returns targets for rule all
    :return:
    """
    output = []
    target1 = expand(os.path.join(out_dir, "minda", "{sample}", "{sample}_minda_ensemble.vcf"), sample=samples)
    output.extend(target1)
    target2 = expand(os.path.join(out_dir,f"{sv_type}_SVs","{sample}_filtered_ensemble.vcf"), sample=samples)
    output.extend(target2)
    if config["annotate"]["activate"]:
        target7 = expand(os.path.join(out_dir,f"{sv_type}_SVs","{sample}.filtered_ensemble.tsv"), sample=samples)
        target8 = expand(os.path.join(out_dir,f"{sv_type}_SVs","split_out",
            "{sample}", "output.filtered.annotated.{split}.tsv"), split=np.arange(0, 20), sample=samples)
        target9 = expand(os.path.join(out_dir,f"{sv_type}_SVs",
             "{sample}.filtered_ensemble.annotated.tsv"), sample=samples)
        target10 = define_caller_table_targets(sample_metadata)
        output.extend(target7)
        output.extend(target8)
        output.extend(target9)
        output.extend(target10)
    if config["visualize"]["activate"]:
        target10  = expand(os.path.join(out_dir,f"{sv_type}_SVs",
             "{sample}.circos.pdf"), sample=samples)
        output.extend(target10)
    return output

