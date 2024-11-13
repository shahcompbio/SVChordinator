"""
reformat minda vcf to format that appeals to sniffles2
"""
minda_vcf = snakemake.inputs["minda_vcf"]
out_vcf = snakemake.outputs["out_vcf"]

# read in vcf
f = open(minda_vcf, "r")
input_vcf = f.readlines()

# remove cmd line
with open(out_vcf, "w+") as f:
    for line in input_vcf:
        # skip offensive line
        if line.startswith("cmd"):
            continue
        else:
            f.write(line)