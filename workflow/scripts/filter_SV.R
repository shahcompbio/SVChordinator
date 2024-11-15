library(tidyr)
library(dplyr)
library(data.table)

# paths
sniffles_tumor_path <- snakemake@input[["tumor_genotypes"]]
sniffles_norm_path <- snakemake@input[["norm_genotypes"]]
somatic_ids <- snakemake@output[["filtered_IDs"]]
read_support <- snakemake@output[["read_support"]]

# normal genotyping
norm_vcf <- read.table(sniffles_norm_path, comment.char = "#")
# tumor genotyping
tumor_vcf <- read.table(sniffles_tumor_path, comment.char = "#")

# Parse Sniffles normal VCF for normal reads
fields <- c("norm_genotype", "norm_genotype_quality", "norm_ref_reads", "norm_variant_reads")

norm_genotypes <- norm_vcf %>%
separate(V10, into = fields, sep = ":", convert = TRUE) %>%
mutate(norm_coverage = norm_ref_reads + norm_variant_reads) %>%
select(V3, norm_ref_reads, norm_variant_reads, norm_coverage)

# parse genotyped tumor for tumor reads
fields <- c("tumor_genotype", "tumor_genotype_quality", "tumor_ref_reads", "tumor_variant_reads")

tumor_genotypes <- tumor_vcf %>%
separate(V10, into = fields, sep = ":", convert = TRUE) %>%
mutate(tumor_coverage = tumor_ref_reads + tumor_variant_reads) %>%
select(V3, tumor_ref_reads, tumor_variant_reads, tumor_coverage)

# Merge on V3, then rename the column in the combined data frame
read_support <- tumor_genotypes %>%
  inner_join(norm_genotypes, by = "V3") %>%
  rename(ID = V3)

"
Let's filter out those SVs with the following:

1. < 5 read coverage in tumor or normal
2. > 0 reads in the normal
3. 0 supporting reads in the tumor
"
somatic_ids <- read_support %>%
  filter(ID %in% read_support$ID,
         norm_variant_reads == 0,
         tumor_variant_reads > 0,
         tumor_coverage >= 5,
         norm_coverage >=5) %>%
  select(ID)

# Write out to file
write.table(somatic_ids, file = "somatic_ids", sep = "\t", col.names = F, row.names = F, quote = F)

# write read support table
write.table(read_support, file = "sniffles_read_support.tsv", sep = "\t", col.names = F, row.names = F, quote = F)