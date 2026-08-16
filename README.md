# Do chromatin loops reach the nearest gene?

A re-analysis of published HiCAR chromatin interaction data from GM12878 cells, testing how often the standard "nearest-gene" assumption matches measured DNA contacts.

*Analysis run 15 August 2026*

---

## Result

**90.5% of testable chromatin loops reach past at least one protein-coding gene.** Only 9.5% connect to the nearest gene. The median loop skips 2 protein-coding genes.

---

## Background

Genes are switched on and off partly by **enhancers** — short stretches of regulatory DNA that
can sit a long way from the gene they control. An enhancer works by physically touching its
target gene's promoter, with the DNA in between bulging out into a **loop**.

When researchers find a regulatory element, they usually cannot see which gene it controls, so
they assign it to **the nearest gene**. This shortcut is used throughout genomics because there is
often no alternative. This analysis asks how well it performs when checked against direct
measurements of DNA contact.

### Motivation

In an earlier project (`pbmc_multiome_peak_gene_links`) I found that nearest-gene assignment
agreed with multiome-derived peak-to-gene links about 63% of the time in PBMCs. That compared
one prediction method against another. This analysis compares a prediction against a physical
measurement of the genome contacting itself.

**The two numbers are not directly comparable** — different cell systems, different definitions
of agreement, different classes of evidence. The PBMC result motivated the question; it is not
combined with these data.

---

## Data

Wei X, Xiang Y, Peters DT, *et al.* HiCAR is a robust and sensitive method to analyze
open-chromatin-associated genome organization. *Molecular Cell* 82(6):1225-1238 (2022).
DOI: 10.1016/j.molcel.2022.01.023

**Supplementary Table S2** (`mmc3.xlsx`), sheet `GM12878 interactions` — the full list of MAPS-called interactions. Note that the preprint and published versions number their supplementary tables differently; S2 is correct for the published paper.

| Property | Value | How verified |
|---|---|---|
| Cell type | GM12878 (B-cell-derived lymphoblastoid line) | sheet name |
| Interactions | 48,516 | `nrow()` |
| Resolution | 5 kb bins | `end1 - start1` constant at 5000 |
| Orientation | all *cis* (same chromosome) | `chr1 == chr2` all TRUE |
| Significance | pre-filtered, FDR < 0.01 | max FDR in file = 0.0099 |
| Genome build | hg38 | stated in paper's STAR Methods |

Gene annotation: `TxDb.Hsapiens.UCSC.hg38.knownGene`, filtered to the 24 main chromosomes and to protein-coding genes only (18,193 genes). Both datasets are hg38, confirmed independently on each side.

---

## Method

1. **Loop span** — distance between the midpoints of the two 5 kb anchor bins.
2. **Filter to testable loops** — required a gene annotated at anchor 2, and a span below 2 Mb. The 2 Mb cap removes chromosome-scale contacts (max span in the raw data was 92.9 Mb), which reflect large-scale genome organisation rather than enhancer-promoter regulation. This threshold is my own judgement call, not one specified by the paper.
3. **Count intervening genes** — built the region between the inner edges of the two anchors (so genes sitting *on* an anchor are not counted as being *between* anchors), then used `countOverlaps()` against the protein-coding gene set.

A loop skipping zero genes reached the nearest protein-coding gene, meaning nearest-gene
assignment would have been correct. A loop skipping one or more reached past closer genes,
meaning nearest-gene assignment would have been wrong.

---

## Results

### Loop spans

- Median span: **145,000 bp**
- **65.1%** of loops span more than 100 kb
- By chromatin state: active **135,000 bp**, poised **175,000 bp**

The paper reports the same direction in H1 hESCs (poised 145 kb vs active 125 kb). Recovering that pattern in a different cell type is a partial replication of a published finding.

### Genes skipped

- Testable loops: **8,344**
- Median span of testable subset: 150,000 bp (close to the full set, so filtering did not bias the sample)
- **90.5%** skip at least one protein-coding gene
- Median genes skipped: **2**
- Most common outcome: **1 gene skipped**

### Robustness

Restricting to the well-defined `active` category alone gives **90.2%**, essentially unchanged from the overall figure. The result is therefore not driven by the `other` category, whose definition I could not locate in the paper's Methods.

Longer loops skip more genes, as expected. However, at any given span the number of genes
skipped varies from 0 to over 20, because gene density differs greatly across the genome.
This is why distance alone cannot answer the question and the gene count is required.

---

## Limitations

These constrain what the result can be claimed to show.

1. **No null model.** I have not run a distance-matched shuffle, so I cannot say how much of the
   90.5% exceeds what random region pairs of the same span distribution would give.
   Genes are common, so some of this figure reflects arithmetic rather than biology. This is the
   most important open question and the obvious next step.

2. **Contact is not proof of regulation.** A loop shows two regions are in physical proximity. It
   does not establish that one controls the other. Functional validation would be required.

3. **GM12878 is not PBMCs.** It is an immortalised, EBV-transformed B-lymphoblastoid line. It
   shares a lineage with one PBMC population but is not primary blood cells. No claim is made
   about the PBMC dataset from these results.

4. **Published loop calls, not raw processing.** I used the authors' MAPS interaction calls. I did
   not process raw HiCAR sequencing data or run the nf-core/hicar pipeline. The loop calls inherit
   all of the original study's processing decisions.

5. **Undefined category.** 26% of loops are classed as `other`, which I could not find defined in
   the Methods. I report but do not interpret this group.

6. **Anchor-gene assignment is the authors'.** The `anchor2_gene` annotation came with the
   supplementary table. The rule used to assign genes to anchors is not stated, which matters
   given that this analysis is itself about how positions get linked to genes.

7. **5 kb resolution.** Anchor positions are known only to within a 5 kb bin.

---

## Files

| File | Contents |
|---|---|
| `hicar_loop_analysis.R` | Complete analysis, runs top to bottom |
| `fig1_loop_spans.png` | Distribution of loop spans |
| `fig2_genes_skipped.png` | Protein-coding genes between anchors (main result) |
| `fig3_span_vs_skipped.png` | Span against genes skipped |

## Reproducing

Download Supplementary Table S2 (`mmc3.xlsx`) from the paper, rename it to
`hicar_table_s2_maps_loops.xlsx`, place it beside the script, and run the script.

```r
install.packages(c("readxl", "ggplot2", "dplyr", "BiocManager"))
BiocManager::install(c("GenomicRanges",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "org.Hs.eg.db"))
```

Written under R 4.6.0.

