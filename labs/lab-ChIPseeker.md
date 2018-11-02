---
layout: default
title:  'ChIP-seq down-stream analysis'
---

# ChIP-seq down-stream analysis: ChIPseeker

### Learning outcomes
Using `ChIPseeker` package
- to profile ChIP signal by genomics location and by ChIP binding to TSS regions
- to annotate peaks, visualise and compare annotations
- to run and compare functional enrichment

## Content
- [Introduction](#Introduction)
- [Data & Methods](#DataMethods)
- [Setting-up](#Setting-up)

- [ChIP profiling](#Profile)
- [Peaks Annotations](#Annotations)
- [Functional analysis](#Functional)

- [Concluding remarks and next steps](#Next)
- [Appendix: figures](#Next)

## Introduction <a name="Introduction"></a>
In this tutorial we use another package, `ChIPseeker`, to have a look at the ChIP profiles, annotate peaks and visualise annotations and run functional enrichment. In a way, `ChIPseeker` can be seen as an alternative and newer workflow to `ChIPpeakAnno`. It also offers additional functionality, e.g. when it comes to visualising ChIP profiles and comparing functional annotations.

_It supports annotating ChIP peaks and provides functions to visualize ChIP peaks coverage over chromosomes and profiles of peaks binding to TSS regions. Comparison of ChIP peak profiles and annotation are also supported. Moreover, it supports evaluating significant overlap among ChIP-seq datasets. Currently, ChIPseeker contains 17,000 bed file information from GEO database. These datasets can be downloaded and compare with user’s own data to explore significant overlap datasets for inferring co-regulation or transcription factor complex for further investigation._

## Data & Methods <a name="DataMethods">
We will build upon the main labs, using the same dataset and results from `DiffBind` analyses saved under `DiffBind.RData`. The tutorial is based on the `ChIPseeker` package tutorial [html](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)

## Setting-up  <a name="Setting-up">
If you have not done it already, install R and R-Studio. Refer back to [pre-course](../precourse) preparations for instructions.

You can continue working in the `diffBind` directory. We need access to `diffBind.RData` object & some libraries

```bash

# Load libraries (install if needed)
library(DiffBind)
library(ChIPseeker)
library(ReactomePA)
library(clusterProfiler)

library(org.Hs.eg.db)  
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

```

## ChIP profile <a name="Profile"></a>
### ChIP peaks coverage plot <a name="Profile_coverage"></a>

After peak calling one may want to visualise distribution of peaks locations over the whole genome. Function `covplot` calculates coverage of peaks regions over chromosomes.

Let's use data saved in `DiffBind.RData` objects. From this object we can easily extract peaks called for all our libraries as well as consensus peakset. In principle, we could also use `ChIPseeker` on raw .BED files.

```bash

# Let's start fresh removing all objects from R environment
rm(list = ls())

# Loading diffBind.RData
load("diffBind.RData")

# Do you remember what was have we saved in the diffBind.RData
ls()

# res.cnt3 object was the final one containing consensus peaks and differential binding results

# to view all samples
dba.show(res.cnt3)

# this should show you our 8 libraries
> dba.show(res.cnt3)
          ID Tissue Factor Replicate Caller Intervals FRiP
1 REST_chip1   HeLa   REST         1 counts      6518 0.11
2 REST_chip2   HeLa   REST         2 counts      6518 0.08
3 REST_chip3 neural   REST         1 counts      6518 0.07
4 REST_chip4 neural   REST         2 counts      6518 0.09
5 REST_chip5  HepG2   REST         1 counts      6518 0.08
6 REST_chip6  HepG2   REST         2 counts      6518 0.06
7 REST_chip7  sknsh   REST         1 counts      6518 0.10
8 REST_chip8  sknsh   REST         2 counts      6518 0.06

```

To plots peaks over genomic locations we need to extract from `res.cnt3` peaks of interest, e.g. consensus peaks or present in a single replicate etc. Here, we will focus on peaks present in HeLa replicates.

```bash

# extracting consensus peaks set with 6518 peaks
peaks.consensus <- dba.peakset(res.cnt3, bRetrieve = T)

# extracting HeLA peaks
peaks.HeLa_rep1 <- peaks.consensus[res.cnt3$called[,1]==1] # peaks called in rep 1
peaks.HeLa_rep2 <- peaks.consensus[res.cnt3$called[,2]==1] # peaks called in rep 2

# adding an unified affinity scores column (re-formatting data)
peaks.HeLa_rep1$Score <- peaks.HeLa_rep1$REST_chip1
peaks.HeLa_rep2$Score <- peaks.HeLa_rep2$REST_chip2

# plotting coverage for replicate 1, using affinity scores as a weight for peaks height
covplot(peaks.HeLa_rep1, weightCol = "Score")

```

We can also compare peaks across replicates. This should give us visual assessment of variability between replicates: peaks locations and strength should match in an ideal scenario

```bash

# creating genomicRangesList object holding replicates 1 and 2
grL.HeLa <- GenomicRangesList(HeLa_rep1=peaks.HeLa_rep1, HeLa_rep2=peaks.HeLa_rep2)

# plotting using affinity scores as a weight for peaks height
covplot(grL.HeLa, weightCol = "Score")

```

What do you think?
- are these peaks reproducible?
- which pair of replicates is most consistent, HeLA, neural, HepG2 or sknsh?
- why is it good to always look at the data instead of simply trusting the output of the summary statistics, after all, we do rely on diffBind to call peaks being consistent


### Profile of ChIP peaks binding to TSS regions
For calculating the profile of ChIP peaks binding to TSS regions, we need to prepare the TSS regions, which are defined as the flanking sequence of the TSS sites. Then we can align the peaks that are mapping to these regions, and generate the tagMatrix used for plotting.

Here, we will select peaks present per cell type, i.e. found in two replicates. We will also create tagMatrix list to enable groups comparisons across cell lines.

```bash

# extracting peaks for each cell line present across replicates
peaks.HeLa <- peaks.consensus[res.cnt3$called[,1]==1 & res.cnt3$called[,2]==1]
peaks.neural <- peaks.consensus[res.cnt3$called[,3]==1 & res.cnt3$called[,4]==1]
peaks.HepG2 <- peaks.consensus[res.cnt3$called[,5]==1 & res.cnt3$called[,6]==1]
peaks.sknsh <- peaks.consensus[res.cnt3$called[,7]==1 & res.cnt3$called[,8]==1]

# TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# calculating tagMatrix
tagMatrix.1 <- getTagMatrix(peaks.HeLa, windows=promoter)
tagMatrix.2 <- getTagMatrix(peaks.neural, windows=promoter)
tagMatrix.3 <- getTagMatrix(peaks.HepG2, windows=promoter)
tagMatrix.4 <- getTagMatrix(peaks.sknsh, windows=promoter)

# preparing tagMatrix list to enable cell lines comparisions
tagMatrixList <- list(HeLa=tagMatrix.1, neural=tagMatrix.2, HepG2=tagMatrix.3, sknsh=tagMatrix.4)

# plotting tagMatrix heatmaps for each cell line
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

# plotting average profile of ChIP peaks among different cell lines
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

```

## Peaks Annotations

## Functional analysis

## Concluding remarks and next steps <a name="Next">


----

## Appendix: figures <a name="Appendix">

![fig-covplot](../figures/lab-chipseeker/fig-covplot-1.pdf)

Fig: Coverage plot for HeLa replicate 1 peaks

----

![fig-covplot-cmp](../figures/lab-chipseeker/fig-covplot-cmp-1.pdf)

Fig: Coverage plot for HeLa replicate 1 and 2 peaks

----

![fig-TSS-heatmap](../figures/lab-chipseeker/fig-heatmap-1.pdf)

Fig: Heatmap of ChIP peaks among different cell lines

----

![fig-TSS-avg-profile](../figures/lab-chipseeker/fig-TSS-avg-profile-1.pdf)

Fig: Average profile of ChIP peaks among different cell lines


----
