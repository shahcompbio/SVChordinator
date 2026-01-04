# merge SVs
def _write_minda_input(wildcards):
    """
    write minda input tsv for a given sample
    :param wildcards:
    :return: path to minda input tsv
    """
    sample = wildcards.sample
    minda_tsv_path = os.path.join(out_dir, "minda", sample, sample + "_minda_input.tsv")
    os.makedirs(os.path.dirname(minda_tsv_path), exist_ok=True)
    df = pd.read_csv(sample_metadata, sep="\t")
    sample_df = df[df["sample"] == sample]
    minda_df = sample_df[["vcf_path", "caller", "tech"]]
    minda_df.to_csv(minda_tsv_path, sep="\t", index=False, header=False)
    return minda_tsv_path

rule sv_merge:
    input:
        tsv=_write_minda_input
    output:
        merged_vcf=os.path.join(out_dir, "minda", "{sample}","{sample}_minda_ensemble.vcf")
    params:
        out_dir=os.path.join(out_dir, "minda", "{sample}"),
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
        --sample_name {wildcards.sample} --min_support {params.min_support} \
        --tolerance {params.tolerance} --min_size {params.min_size} \
        --bed {params.filter_bed}
        """


