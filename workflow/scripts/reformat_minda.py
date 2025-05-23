"""
reformat minda vcf to format that appeals to sniffles2
"""
minda_vcf = snakemake.input["minda_vcf"]
out_vcf = snakemake.output["out_vcf"]

# read in vcf
f = open(minda_vcf, "r")
input_vcf = f.readlines()

# remove cmd line
with open(out_vcf, "w+") as f:
    for line in input_vcf:
        # skip offensive line
        if line.startswith("cmd"):
            continue
        # replace quality with dummy number for sniffles
        elif not line.startswith("#"):
            terms = line.split("\t")
            terms[5] = "30"
            newline = ("\t").join(terms)
            f.write(newline)
        else:
            f.write(line)