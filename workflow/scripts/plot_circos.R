library(RCircos)
library(dplyr)

# paths
sv.data.path <- snakemake@input[["annotated_SVs"]]
circos.path <- snakemake@output[["circos_plot"]]

sv.data <- read.csv(sv.data.path, sep="\t")

# make links database
sv.link.data <- sv.data[, c("chrom1", "base1", "base1", "chrom2", "base2", "base2")]
colnames(sv.link.data) <- c("Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1")

# load in gene info
sv.gene.data <- sv.data[, c("chrom1", "base1", "base1", "gene_name_1")]
new_cols <- c("Chromosome", "chromStart", "chromEnd", "Gene")
colnames(sv.gene.data) <- new_cols
sv.gene.data.1 <- sv.data[, c("chrom2", "base2", "base2", "gene_name_2")]
colnames(sv.gene.data.1) <- new_cols
sv.gene.data <- dplyr::bind_rows(sv.gene.data, sv.gene.data.1)
sv.gene.data <- sv.gene.data %>%
  filter(!(sv.gene.data$Gene == "")) %>%
  distinct(Gene, .keep_all = TRUE)

# give colors based on SV type ... let's do a color-blind friendly palette
pal <- palette.colors(palette = "Okabe-Ito")
link.colors <- c()
sv_types <- sv.data$SV_Type
link.widths <- rep(1, nrow(sv.data))
for (i in 1:length(sv_types)){
  sv_type <- sv_types[i]
  if (sv_type == "INS"){
    link.color <- pal[1]
  } else if (sv_type == "DEL") {
     link.color <- pal[2]
    }
     else if (sv_type == "BND"){
      link.color <- pal[3]
    } else if (sv_type == "INV"){
      link.color <- pal[4]
    } else if (sv_type == "DUP"){
      link.color <- pal[5]
  }
  link.colors <- c(link.colors, link.color)
  link.widths[i] <- 1.5
}

# plot
out.file <- circos.path;
pdf(file=out.file, height=8, width=8, compress=TRUE);
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info,
    tracks.inside=10, tracks.outside=0 )
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
# add SVs
track.num <- 4;
# set link colors
sv.link.data["PlotColor"] <- link.colors
RCircos.Link.Plot(sv.link.data, track.num, FALSE, lineWidth=link.widths);
# add gene annotations ....
name.col <- 4;
side <- "in";
track.num <- 1;
RCircos.Gene.Connector.Plot(sv.gene.data,track.num, side);
track.num <- 2;
RCircos.Gene.Name.Plot(sv.gene.data,name.col,track.num, side);

# Add a legend for the colors
legend("topright", legend = c("Insertion", "Deletion", "Translocation", "Inversion", "Duplication"),
       fill = pal[1:5],
       title = "Band Colors", cex = 0.8)

