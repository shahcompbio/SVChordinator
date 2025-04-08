# merge SVs
rule sv_merge:
    input:
        tsv=minda_tsv
    output:
        merged_vcf=os.path.join(out_dir, "minda", sample_name + "_minda_ensemble.vcf")
    params:
        sample_name=config["sample_name"],
        out_dir=os.path.join(out_dir, "minda"),
        filter_bed=config["filter_bed"],
        min_support=min_callers,
        tolerance=100,
        min_size=50
    resources:
        mem_mb = 20000,
        time = 360,
        retries = 0
    threads: 1,
    container:
        "docker://quay.io/preskaa/minda:v250408",
    shell:
        """
        /minda/minda.py ensemble --tsv {input.tsv} --out_dir {params.out_dir} \
        --sample_name {params.sample_name} --min_support {params.min_support} \
        --tolerance {params.tolerance} --min_size {params.min_size} \
        --bed {params.filter_bed}
        """


