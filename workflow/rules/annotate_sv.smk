# convert sniffles-format vcf to a tsv
rule convert_ensemble_vcf:
    input:
        vcf = os.path.join(out_dir, "somatic_SVs", sample_name + "_filtered_ensemble.vcf"),
    output:
        tsv = os.path.join(out_dir,"somatic_SVs",sample_name+ ".filtered_ensemble.tsv"),
        split_tsv = expand(os.path.join(out_dir,"somatic_SVs", "split_out",
            sample_name, "output.filtered.{split}.tsv"), split=np.arange(0, 20))
    params:
        split_out = os.path.join(out_dir,"somatic_SVs","split_out", sample_name)
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
        oncokb = config["annotate"]["oncokb"],
        annotations = config["annotate"]["gene_annotation"]
    resources:
        time = 360,
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
        time = 60,
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
        tsv = os.path.join(out_dir, "raw_SVs", sample_name, sample_name+ ".{caller}.tsv")
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

rule annotate_svtypes:
    input:
        all_SVs = os.path.join(out_dir,"somatic_SVs",
            sample_name + ".filtered_ensemble.temp.annotated.tsv"),
        caller_tables = expand(os.path.join(out_dir, "raw_SVs",
            sample_name, sample_name+ ".{caller}.tsv"), caller=callers),
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


