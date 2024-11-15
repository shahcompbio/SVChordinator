# reformat minda output to conform to sniffles specifications
# we also need to give them proper SV types as well to properly genotype
rule reformat_minda:
    input:
        minda_vcf = os.path.join(out_dir,"minda",sample_name + "_minda_ensemble.vcf"),
    output:
        out_vcf = os.path.join(out_dir, "sniffles", sample_name + "_sniffles_ensemble.vcf")
    threads: 1
    resources:
        mem_mb = 4000,
        time = 20,
        retries = 0
    container:
        "docker://quay.io/preskaa/biopython:v241011a"
    script:
        "../scripts/reformat_minda.py"
# genotype the consensus SVs in the normal to extract read support
rule genotype_normal:
    input:
        merged_vcf=os.path.join(out_dir, "sniffles", sample_name + "_sniffles_ensemble.vcf"),
        normal_bam=config["genotype"]["normal_bam"]
    output:
        vcf = os.path.join(out_dir, "sniffles", sample_name + "_normal_genotypes.vcf")
    threads: 10,
    container:
        "docker://quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0",
    shell:
        """
        sniffles --input {input.normal_bam} --genotype-vcf {input.merged_vcf} \
        --vcf {output.vcf} --threads {threads}
        """

# genotype the consensus SVs in the tumor
# the issue with this rule is that the SV types need to be correct to genotype properly
rule genotype_tumor:
    input:
        merged_vcf=os.path.join(out_dir, "sniffles", sample_name + "_sniffles_ensemble.vcf"),
        tumor_bam=config["genotype"]["tumor_bam"]
    output:
        vcf = os.path.join(out_dir, "sniffles", sample_name + "_tumor_genotypes.vcf")
    threads: 10,
    container:
        "docker://quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0",
    shell:
        """
        sniffles --input {input.tumor_bam} --genotype-vcf {input.merged_vcf} \
        --vcf {output.vcf} --threads {threads}
        """
# filter SVs post genotyping
rule filter_SV:
    input:
        tumor_genotypes = os.path.join(out_dir,"sniffles",sample_name + "_tumor_genotypes.vcf"),
        norm_genotypes = os.path.join(out_dir, "sniffles", sample_name + "_normal_genotypes.vcf")
    output:
        read_support=os.path.join(out_dir, "sniffles", sample_name + "_read_support.tsv"),
        filtered_IDs = os.path.join(out_dir, "sniffles", "filtered_IDs.tsv")
    threads: 1,
    resources:
        mem_mb = 4000,
        time = 20,
        retries = 0
    container:
        "docker://quay.io/biocontainers/r-tidyverse:1.2.1"
    script:
        "../scripts/filter_SV.R"