import os
import numpy as np
import pandas as pd

### paths ###
out_dir = config["out_dir"]
minda_tsv = config["samples"]
sample_name = config["sample_name"]
ONT = config["ONT"]
oncokb = config["annotate"]["oncokb"]
gene_annotations = config["annotate"]["gene_annotation"]
### helper functions ###
def extract_individual_calls(minda_tsv):
    """
    extract a dict of info for each individual set of SV calls
    :param minda_tsv: input file for minda
    :return: pandas dataframe of vcf paths, names of individual callers
    """
    df = pd.read_csv(minda_tsv, sep="\t", header=None)
    df.columns = ["vcf_path", "caller", "nickname"]
    callers = list(df["caller"])
    return df, callers

# extract individual caller paths
caller_df, callers = extract_individual_calls(minda_tsv)
print(callers)

def get_output():
    """
    returns targets for rule all
    :return:
    """
    output = []
    target1 = os.path.join(out_dir, "minda", sample_name+"_minda_ensemble.vcf")
    output.append(target1)
    if config["genotype"]["activate"]:
        target3 = os.path.join(out_dir, "sniffles", sample_name + "_normal_genotypes.vcf")
        target4 = os.path.join(out_dir,"sniffles",sample_name + "_tumor_genotypes.vcf")
        target5 = os.path.join(out_dir, "sniffles", sample_name + "_filtered_IDs.tsv")
        target6 = os.path.join(out_dir, "somatic_SVs", sample_name + "_filtered_ensemble.vcf")
        output.append(target3)
        output.append(target4)
        output.append(target5)
        output.append(target6)
    # just reformat minda if genotyping hasn't been run
    else:
        target1 = os.path.join(out_dir,"somatic_SVs",sample_name + "_filtered_ensemble.vcf")
        output.append(target1)
    if config["annotate"]["activate"]:
        target7 = os.path.join(out_dir,"somatic_SVs",sample_name+ ".filtered_ensemble.tsv")
        target8 = expand(os.path.join(out_dir,"somatic_SVs","split_out",
            sample_name, "output.filtered.annotated.{split}.tsv"), split=np.arange(0, 20))
        target9 = os.path.join(out_dir,"somatic_SVs",
             sample_name+ ".filtered_ensemble.annotated.tsv")
        target10 = expand(os.path.join(out_dir, "raw_SVs", sample_name, sample_name+ ".{caller}.tsv"),
            caller=callers)
        output.append(target7)
        output.extend(target8)
        output.append(target9)
        output.append(target10)
    if config["visualize"]["activate"]:
        target10  = os.path.join(out_dir,"somatic_SVs",
             sample_name + ".circos.pdf")
        output.append(target10)
    return output

