rule plot_circos:
    input:
        annotated_SVs = os.path.join(out_dir,"somatic_SVs",
             sample_name + ".filtered_ensemble.annotated.tsv")
    output:
        circos_plot = os.path.join(out_dir,"somatic_SVs",
             sample_name + ".circos.pdf")
    container:
        "docker://quay.io/biocontainers/r-rcircos:1.1.3--r3.2.2_0"
    threads: 1,
    resources:
        time = 60,
        mem_mb = 20000,
        retries = 0,
    script:
        "../scripts/plot_circos.R"