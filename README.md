# Do chromatin loops reach the nearest gene?

This is a re-analysis of published HiCAR chromatin interaction data from GM12878 cells, asking a scientific question: When two regions of DNA physically contact each other, does the contact usually involve
the nearest protein-coding gene?

## Result

**90.5% of testable chromatin loops pass at least one protein-coding gene before reaching the gene at the other anchor.**

Only 9.5% reach the nearest protein-coding gene without passing another one first. The median testable loop passes 2 protein-coding genes.

That does not mean that 90.5% of enhancers regulate a non-nearest gene. A physical contact is not
the same thing as a regulatory relationship. What it does show is that chromatin contacts
frequently cross over genes that would look closer by genomic distance alone.

## Why look at this?

Enhancers are regulatory regions that can influence gene expression from far away. The enhancer and
promoter can be brought together through a chromatin loop, allowing regions that are separated
along the linear genome to physically interact.

The problem is that when researchers identify a regulatory region, they often do not immediately
know which gene it affects. One common shortcut is to assign the region to the nearest gene.

It is a useful shortcut, but it is still a shortcut.

I wanted to see how well that assumption holds up against a direct measurement of chromatin contact.

This started from an earlier analysis I did with PBMC multiome data, where nearest-gene assignment
agreed with multiome-derived peak-to-gene links about 63% of the time. That result compared two
computational approaches. Here the question is different: does genomic proximity agree with measured
physical interaction?

The two percentages should not be compared directly. They come from different datasets, cell systems,
definitions, and types of evidence. The PBMC result is what led me to this question, not something I
am combining with the HiCAR result.

## Data

I used the published HiCAR interaction calls from:

Wei X, Xiang Y, Peters DT, et al. HiCAR is a robust and sensitive method to analyze
open-chromatin-associated genome organization. *Molecular Cell* 82(6):1225-1238 (2022).
DOI: 10.1016/j.molcel.2022.01.023

The main dataset was Supplementary Table S2, `mmc3.xlsx`, specifically the `GM12878 interactions`
sheet. This contains the MAPS-called chromatin interactions for GM12878 cells.

| Property | Value |
|---|---|
| Cell type | GM12878 |
| Interactions | 48,516 |
| Resolution | 5 kb |
| Contact orientation | all cis (both anchors on the same chromosome) |
| FDR | < 0.01 |
| Genome build | hg38 |

I checked the table directly for the interaction count, bin size, chromosome pairing, and FDR. The
genome build is stated in the paper's STAR Methods.

For genes, I used `TxDb.Hsapiens.UCSC.hg38.knownGene`, restricted to the 24 main chromosomes and to protein-coding genes. That left 18,193 genes.

Both the interaction data and the gene annotation are on hg38.

## How I measured it

I wanted the question to be as simple as possible:

**How many protein-coding genes lie between the two anchors of a chromatin loop?**

If there are no genes between the anchors, the loop has not passed a closer protein-coding gene.
Under this definition, nearest-gene assignment would be correct.

If there is one or more, the contact reaches past at least one closer protein-coding gene. If the
gene at the far anchor is the actual target, nearest-gene assignment would point to the wrong one.

**1. Measure the loop**

Each anchor is a 5 kb genomic bin. I calculated the distance between the midpoints of the two anchors.

**2. Remove loops that cannot answer the question**

I kept loops that:

- had a gene annotated at anchor 2
- were shorter than 2 Mb

The 2 Mb cutoff is my own choice. It is not a threshold used by the original paper. I used it to remove very large-scale contacts that are less useful for asking an enhancer-promoter-style question. The original interaction table contains contacts as large as 92.9 Mb.

**3. Count genes between the anchors**

For each remaining loop, I looked only at the region between the inner edges of the two anchor bins.
This matters because a gene sitting directly on an anchor should not be counted as an intervening gene.

I then counted overlaps with the protein-coding gene annotation using `countOverlaps()`.

That gives a straightforward measure:

- 0 genes skipped = the loop reached the nearest protein-coding gene
- 1 or more = the contact passes at least one closer protein-coding gene

## What I found

### Loop size

Across the interaction set:

- Median loop span: 145 kb
- 65.1% span more than 100 kb
- Active loops: 135 kb median
- Poised loops: 175 kb median

The paper reports the same general pattern in H1 hESCs, where poised interactions were longer than
active ones (145 kb vs 125 kb). Seeing the same direction in GM12878 was a good sign that the data
and the basic span calculation were behaving sensibly.

### Genes skipped

After filtering, 8,344 loops were testable. The median span of this subset was 150 kb, very close to the 145 kb median of the full interaction set.

The main result:

- 90.5% passed at least one protein-coding gene
- 9.5% did not
- Median number of genes skipped: 2
- Most common outcome: 1 gene skipped

So in this dataset, it was far more common for a chromatin contact to pass one or more protein-coding
genes than to reach the nearest one first.

### The relationship with distance

Longer loops generally passed more genes, which is not surprising. There is simply more genomic
sequence available to contain genes.

But distance alone does not determine how many genes are crossed. At the same loop span, some
interactions passed 0 genes while others passed more than 20. That is because gene density varies
substantially across the genome.

This is one reason I think the simple distance-based intuition is worth testing, rather than assuming
that a longer distance automatically tells us how likely nearest-gene assignment is to be correct.

### A robustness check

I repeated the main calculation using only loops classified as active. The result was 90.2%, almost identical to the 90.5% overall.

That makes it unlikely that the main result is being driven simply by the `other` category, which makes up about 26% of the loops and whose definition I could not find in the paper's Methods. I report that category but do not try to interpret what `other` represents biologically.

## What this does and does not show

The 90.5% number is interesting, but it is easy to overstate what it means.

**1. There is no null model yet**

This is the biggest limitation.

I have not compared these loops against randomly generated pairs of regions with the same distance
distribution. Genes are common across the genome. If two regions are 150 kb apart, it is already
fairly likely that some gene will happen to lie between them.

So I cannot yet say that 90.5% is evidence of biological targeting away from the nearest gene. As it stands, it tells us what the observed HiCAR interactions look like. A distance-matched null model is the obvious next test.

**2. A loop is not proof of regulation**

HiCAR measures physical proximity. A contact between an enhancer and a region containing a gene does
not prove that the enhancer regulates that gene. Functional experiments would be needed to establish
that relationship.

So the result is about chromatin contact, not definitive enhancer-target relationships.

**3. GM12878 is not PBMC**

GM12878 is an immortalized, EBV-transformed B-lymphoblastoid cell line. It sits in the broader B-cell
context of PBMC datasets, but it is not primary blood cells.

So I am not using this result to make any claim about PBMCs.

**4. I used the published interaction calls**

I did not start from the raw HiCAR sequencing data. I used the MAPS interaction calls provided by the
authors, which means the result inherits the original study's processing and filtering decisions.

**5. The anchor annotation is not completely transparent**

The `anchor2_gene` field was already present in the supplementary table. The paper does not make the
exact rule used to assign that annotation clear enough for me to reproduce it independently.

That matters here, because gene assignment is essentially the question being tested.

**6. The resolution is 5 kb**

The anchors are 5 kb bins, so the exact position of an interaction is not known at single-base
resolution. Whether a gene falls just inside or just outside an anchor can therefore come down to
where the bin boundaries happen to fall.

**7. The 2 Mb cutoff is mine**

The original data contain much larger interactions. I chose 2 Mb because I wanted to focus on contacts
that are more plausibly relevant to local regulatory relationships rather than chromosome-scale
organization.

That choice should be treated as an analysis decision, not as a biological boundary.

## What I would do next

The next step is not to collect another dataset just to get another percentage. It is to answer the
question the current analysis cannot:

**Is 90.5% actually higher than what we would expect from random contacts at the same genomic distances?**

I would generate distance-matched random anchor pairs, keeping the chromosome and ideally the local
genomic context, and count intervening genes the same way.

If random pairs also produce something close to 90%, then the current result is mostly a consequence
of gene density and loop length. If the observed interactions consistently skip substantially more
genes than the matched null, that becomes much more interesting.

From there, the analysis could also ask whether the effect changes with:

- enhancer/promoter annotation
- chromatin state
- gene density
- loop length
- expression level
- directionality of the interaction
- active versus poised chromatin

That would move the project from "nearest-gene assignment often disagrees with chromatin contacts"
toward a more meaningful question about what actually predicts which gene a regulatory region contacts.

## Reproducing the analysis

| File | Purpose |
|---|---|
| `hicar_loop_analysis.R` | Complete analysis script |
| `fig1_loop_spans.png` | Distribution of loop spans |
| `fig2_genes_skipped.png` | Main result: genes between anchors |
| `fig3_span_vs_skipped.png` | Loop span vs genes skipped |

To reproduce it:

1. Download Supplementary Table S2 (`mmc3.xlsx`) from the HiCAR paper.
2. Rename it to `hicar_table_s2_maps_loops.xlsx`.
3. Put it in the same directory as `hicar_loop_analysis.R`.
4. Run the script from top to bottom.

The script runs the full analysis end to end, rather than relying on manually generated intermediate
results.

