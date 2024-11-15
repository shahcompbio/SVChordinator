import os
import numpy as np
# import config
out_dir = config["out_dir"]
minda_tsv = config["samples"]
sample_name = config["sample_name"]

def get_output():
    output = []
    target1 = os.path.join(out_dir, "minda", sample_name+"_minda_ensemble.vcf")
    output.append(target1)
    if config["genotype"]["activate"]:
        target3 = os.path.join(out_dir, "sniffles", sample_name + "_normal_genotypes.vcf")
        target4 = os.path.join(out_dir,"sniffles",sample_name + "_tumor_genotypes.vcf")
        target5 = os.path.join(out_dir, "sniffles", sample_name + "filtered_IDs.tsv")
        target6 = os.path.join(out_dir, "somatic_SVs", sample_name + "_filtered_ensemble.vcf")
        output.append(target3)
        output.append(target4)
        output.append(target5)
        output.append(target6)
    if config["annotate"]["activate"]:
        target7 = os.path.join(out_dir,"somatic_SVs",sample_name+ ".filtered_ensemble.tsv")
        target8 = expand(os.path.join(out_dir,"somatic_SVs","split_out",
            sample_name, "output.filtered.annotated.{split}.tsv"), split=np.arange(0, 20))
        target9 = os.path.join(out_dir,"somatic_SVs",
             sample_name+ ".filtered_ensemble.annotated.tsv")
        output.append(target7)
        output.extend(target8)
        output.append(target9)
    return output
