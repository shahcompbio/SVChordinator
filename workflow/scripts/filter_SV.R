library(tidyverse)
library(data.table)


# Severus VCF
severus_vcf <- read.table(snakemake@input$minda_vcf, comment.char = "#")
# severus_vcf <- read.table("results/SPECTRUM-OV-004_S1_LEFT_ADNEXA_R1/severus/SPECTRUM-OV-004_S1_LEFT_ADNEXA_R1_PASS.vcf", comment.char = "#")

# Sniffles genotyped normal VCF
sniffles_vcf <- read.table(snakemake@input$sniffles_normal_vcf, comment.char = "#")
# sniffles_vcf <- read.table("results/SPECTRUM-OV-004_S1_LEFT_ADNEXA_R1/severus/SPECTRUM-OV-004_S1_LEFT_ADNEXA_R1_PASS_normal_genotypes.vcf", comment.char = "#")

# Sniffles genotyped tumor VCF
sniffles_tumor_vcf <- read.table(snakemake@input$sniffles_tumor_vcf, comment.char = "#")

# Parse Sniffles normal VCF for normal reads and VAF
sniffles_vcf_parsed <- sniffles_vcf %>%
 mutate(norm_ref_reads = as.numeric(str_split(V10, ":", simplify = TRUE)[, 3]),
     norm_variant_reads = as.numeric(str_split(V10, ":", simplify = TRUE)[, 4])) %>%
 filter(!(norm_ref_reads == 0 & norm_variant_reads == 0)) %>%
 mutate(norm_vaf = norm_variant_reads/(norm_ref_reads + norm_variant_reads))

# also compute these statistics for the genotyped tumor vcf
sniffles_tumor_vcf_parsed <- sniffles_vcf %>%
 mutate(tumor_ref_reads = as.numeric(str_split(V10, ":", simplify = TRUE)[, 3]),
     tumor_variant_reads = as.numeric(str_split(V10, ":", simplify = TRUE)[, 4])) %>%
 filter(!(tumor_ref_reads == 0 & tumor_variant_reads == 0)) %>%
 mutate(tumor_vaf = tumor_variant_reads/(tumor_ref_reads + tumor_variant_reads))

# Filter for IDs with 0 normal reads and in Severus VCF
sniffles_somatic_id <- sniffles_vcf_parsed %>%
 filter(V3 %in% severus_vcf$V3,
     norm_variant_reads == 0) %>%
 select(V3)

# Write out to file
write.table(sniffles_somatic_id, file = snakemake@output$somatic, sep = "\t", col.names = F, row.names = F, quote = F)