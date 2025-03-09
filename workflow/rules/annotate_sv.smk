include: "common.smk"
# reformat minda vcf if it hasn't been done already
rule reformat_minda_no_genotype:
    input:
        minda_vcf = os.path.join(out_dir,"minda",sample_name + "_minda_ensemble.vcf"),
    output:
        out_vcf = os.path.join(out_dir, "somatic_SVs", sample_name + "_filtered_ensemble.vcf")
    threads: 1
    resources:
        mem_mb = 4000,
        time = 20,
        retries = 0
    container:
        "docker://quay.io/preskaa/biopython:v241011a"
    script:
        "../scripts/reformat_minda.py"
# convert sniffles-format vcf to a tsv
rule convert_ensemble_vcf:
    input:
        vcf = os.path.join(out_dir, "somatic_SVs", sample_name + "_filtered_ensemble.vcf"),
    output:
        tsv = os.path.join(out_dir,"somatic_SVs",sample_name+ ".filtered_ensemble.tsv"),
        split_tsv = expand(os.path.join(out_dir,"somatic_SVs", "split_out",
            sample_name, "output.filtered.{split}.tsv"), split=np.arange(0, 20))
    params:
        split_out = os.path.join(out_dir,"somatic_SVs","split_out", sample_name),
        ref = ideo
    container:
        "docker://quay.io/preskaa/annotate_genes:v240817"
    script:
        "../scripts/vcf2csv.py"

# annotate those genes ... marginally faster ...
rule annotate_genes:
    input:
        split_tsv = os.path.join(out_dir,"somatic_SVs", "split_out",
            sample_name, "output.filtered.{split}.tsv")
    output:
        split_tsv = os.path.join(out_dir,"somatic_SVs","split_out",
            sample_name, "output.filtered.annotated.{split}.tsv")
    params:
        oncokb = oncokb,
        annotations = gene_annotations
    resources:
        time = 360,
        retries = 1,
    singularity:
        "docker://quay.io/preskaa/annotate_genes:v240817"
    script:
        "../scripts/annotate_genes_faster.py"

# merge gene annotation outputs into a unified dataframe ....
rule merge_annotated_SVs:
    input:
        split_tsvs = expand(os.path.join(out_dir,"somatic_SVs", "split_out",
            sample_name, "output.filtered.annotated.{split}.tsv"), split=np.arange(0, 20))
    output:
        all_SVs = temp(os.path.join(out_dir,"somatic_SVs",
             sample_name+ ".filtered_ensemble.temp.annotated.tsv"))
    container:
        "docker://quay.io/preskaa/annotate_genes:v240817"
    threads: 1,
    resources:
        time = 20,
        mem_mb = 10000,
        retries = 0,
    script:
        "../scripts/merge_annotated_svs.py"

# attempt to annotate any "unresolved" breakpoints (i.e., BND variants)
# https://gatk.broadinstitute.org/hc/en-us/articles/9022476791323-Structural-Variants

# convert vcfs to tsvs from each caller

def _fetch_vcf(wildcards):
    df = caller_df[caller_df["caller"] == wildcards.caller]
    assert len(df) == 1, f"{len(df)} vcfs for {wildcards.caller}"
    return list(df["vcf_path"])[0]



rule variants2table:
    input:
        vcf = _fetch_vcf
    output:
        tsv = os.path.join(out_dir, "raw_SVs",
            sample_name, "ONT", sample_name+ ".{caller}.tsv")
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    threads: 1
    resources:
        time = 20,
        mem_mb = 4000,
        retries = 0,
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/STRANDS\t%INFO/BP_NOTATION\n' {input.vcf} -u -H -o {output.tsv}
        """

# convert illumina variants to table
rule ILL_variants2table:
    input:
        vcf = _fetch_vcf
    params:
        vcf = temp(os.path.join(out_dir,"raw_SVs",
            sample_name, "ILL", sample_name + ".{caller}.vcf")),
    output:
        tsv = os.path.join(out_dir, "raw_SVs",
            sample_name, "ILL", sample_name+ ".{caller}.tsv")
    container:
        "docker://quay.io/preskaa/viola-sv:1.0.2"
    threads: 1
    resources:
        time = 20,
        mem_mb = 4000,
        retries = 0,
    script:
        "../scripts/viola_sv.py"


rule annotate_svtypes:
    input:
        all_SVs = os.path.join(out_dir,"somatic_SVs",
            sample_name + ".filtered_ensemble.temp.annotated.tsv"),
        caller_tables = define_caller_table_target(callers),
    output:
        all_SVs = os.path.join(out_dir,"somatic_SVs",
            sample_name + ".filtered_ensemble.annotated.tsv")
    container:
        "docker://quay.io/preskaa/annotate_genes:v240817"
    threads: 1,
    resources:
        time = 60,
        mem_mb = 10000,
        retries = 0,
    script:
        "../scripts/annotate_svtypes.py"


