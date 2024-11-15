## convert sniffles vcf to a tsv
rule convert_vcf:
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
        all_SVs = os.path.join(out_dir,"somatic_SVs",
             sample_name+ ".filtered_ensemble.annotated.tsv")
    container:
        "docker://quay.io/preskaa/annotate_genes:v240817"
    threads: 1,
    resources:
        time = 60,
        mem_mb = 10000,
        retries = 0,
    script:
        "../scripts/merge_annotated_svs.py"