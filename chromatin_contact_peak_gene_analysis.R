library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


excel_sheets("hicar_table_s2_maps_loops.xlsx")

peek <- read_excel("hicar_table_s2_maps_loops.xlsx",
                   sheet = "GM12878 interactions",
                   n_max = 5)
peek

gm <- read_excel("hicar_table_s2_maps_loops.xlsx",
                 sheet = "GM12878 interactions")
dim(gm)
names(gm)

gm$`interaction type`

table(gm$end1 - gm$start1)

sum(!is.na(gm$anchor1_gene))
sum(!is.na(gm$anchor2_gene))
sum(!is.na(gm$anchor1_gene) | !is.na(gm$anchor2_gene))

gm$span <- ((gm$start2 + gm$end2) / 2) - ((gm$start1 + gm$end1) / 2)

summary(gm$span)

mean(gm$span > 100000)

ggplot(gm, aes(x = span)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  labs(title = "How far do GM12878 HiCAR loops reach?",
       x = "Distance between loop anchors (bp, log scale)",
       y = "Number of loops")


tapply(gm$span, gm$`interaction type`, median)

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")


all_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
length(all_genes)
all_genes[1:3]

main_chr <- paste0("chr", c(1:22, "X", "Y"))
genes_main <- all_genes[seqnames(all_genes) %in% main_chr]

length(genes_main)


gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                    keys = names(genes_main),
                    keytype = "ENTREZID",
                    columns = c("SYMBOL", "GENETYPE"))

head(gene_info)
table(gene_info$GENETYPE)

coding_ids <- gene_info$ENTREZID[gene_info$GENETYPE == "protein-coding"]

genes_coding <- genes_main[names(genes_main) %in% coding_ids]

length(genes_coding)

symbol_lookup <- gene_info$SYMBOL[match(names(genes_coding), gene_info$ENTREZID)]
genes_coding$symbol <- symbol_lookup

head(genes_coding)

testable <- gm[!is.na(gm$anchor2_gene) & gm$span < 2000000, ]
nrow(testable)

summary(testable$span)
table(testable$`interaction type`)

gaps <- GRanges(
  seqnames = testable$chr1,
  ranges = IRanges(start = testable$end1, end = testable$start2)
)

length(gaps)
gaps[1:3]

genome(gaps) <- "hg38"
gaps[1:3]

n_genes_between <- countOverlaps(gaps, genes_coding)

table(n_genes_between)
mean(n_genes_between > 0)

summary(n_genes_between)

plot(testable$span, n_genes_between,
     xlab = "Loop span (bp)", ylab = "Protein-coding genes skipped",
     pch = ".", log = "x")


fig1 <- ggplot(gm, aes(x = span)) +
  geom_histogram(bins = 60, fill = "grey30") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "GM12878 HiCAR loops typically span ~145 kb",
       subtitle = "48,516 MAPS interactions, FDR < 0.01 (Wei et al. 2022, Table S2)",
       x = "Distance between loop anchors (bp, log scale)",
       y = "Number of loops") +
  theme_minimal()

ggsave("fig1_loop_spans.png", fig1, width = 7, height = 4.5, dpi = 300)

skip_counts <- data.frame(genes_skipped = n_genes_between)

fig2 <- ggplot(skip_counts, aes(x = genes_skipped)) +
  geom_bar(fill = "grey30") +
  coord_cartesian(xlim = c(-0.5, 15)) +
  labs(title = "90.5% of loops reach past at least one protein-coding gene",
       subtitle = "8,344 testable GM12878 loops (gene annotated at anchor 2, span < 2 Mb); median 2 genes skipped",
       x = "Protein-coding genes between loop anchors",
       y = "Number of loops") +
  theme_minimal()


ggsave("fig2_genes_skipped.png", fig2, width = 7, height = 4.5, dpi = 300)

scatter_data <- data.frame(span = testable$span,
                           genes_skipped = n_genes_between)

fig3 <- ggplot(scatter_data, aes(x = span, y = genes_skipped)) +
  geom_point(alpha = 0.1, size = 0.8) +
  scale_x_log10(labels = scales::comma) +
  labs(title = "Longer loops skip more genes, but gene density varies widely",
       subtitle = "At any given span, loops skip anywhere from 0 to 20+ genes",
       x = "Loop span (bp, log scale)",
       y = "Protein-coding genes skipped") +
  theme_minimal()


ggsave("fig3_span_vs_skipped.png", fig3, width = 7, height = 4.5, dpi = 300)





