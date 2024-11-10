import os
out_dir = config["out_dir"]
# merge SVs
rule sv_merge:
    input:
        tsv=config["samples"]
    output:
        merged_vcf=os.path.join(out_dir, "minda", config["sample_name"]+"_minda_ensemble.vcf")
    params:
        sample_name=config["sample_name"],
        out_dir=os.path.join(out_dir, "minda"),
        min_support=2,
        tolerance=100,
        min_size=50
    resources:
        mem_mb = 20000,
        time = 360,
        retries = 0
    threads: 1,
    container:
        "docker://quay.io/preskaa/minda:v241109",
    shell:
        """
        ./minda.py ensemble --tsv {input.tsv} --outdir {params.out_dir} \
        --min_support {params.min_support} --tolerance {params.tolerance} \
        --min_size {params.min_size}
        """


