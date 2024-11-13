import os
out_dir = config["out_dir"]
minda_tsv = config["samples"]
sample_name = config["sample_name"]

def get_output():
    output = []
    target1 = os.path.join(out_dir, "minda", sample_name+"_minda_ensemble.vcf")
    output.append(target1)
    if config["genotype"]["activate"]:
        target3 = os.path.join(out_dir, "sniffles", sample_name + "_normal_genotypes.vcf")
        output.append(target3)
    return output
