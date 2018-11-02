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
We will build upon the main labs, using the same dataset and results from `DiffBind` analyses saved under `DiffBind.RData`

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

# res.cnt3 object was the final one containing consensus peakset and differential binding results

# to view all samples
dba.show(res.cnt3)

# should give you our 8 libraries
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
- are the peaks reproducible?
- which pair of replicates is most consistent, HeLA, neural, HepG2 or sknsh?
- why is good to always look at the data instead of simply trusting the output of the summary statistics, after all, we do rely on diffBind to call peaks being consistent





###

### Setting-up `DiffBind`  <a name="DB_install">
Create a separate directory on your local computer where you would like to work and name it e.g. `diffBind`.

To be able to run `DiffBind` locally, we will need access to a set of files
- BAM files
- BED files with called peaks regions
- sample sheet information `.txt` file

Download these from Uppmax with `_scp_`command:

```bash

scp -r <username>@rackham.uppmax.uu.se:~/chipseq/data/bam/* .
scp -r <username>@rackham.uppmax.uu.se:~/chipseq/results/peaks_bed/* .
scp -r <username>@rackham.uppmax.uu.se:~/chipseq/analysis/R/samples_REST.txt .

```

You may want to place the downloaded files in the `diffBind` directory or at least keep a track of their location.

Also we need to modify `samples_REST.txt` so the pathways are pointing to the BAM and BED files on your local computer. Adjust the pathways in any editor your like.

Now, we can open R-Studio and set working directory to working folder e.g. `diffBind` folder by `Session -> Set Working Directory -> Choose Directory`. Now, all R commands will be in respect to this directory.

You can type commands directly in the Console window. A bit smarter way is to open a new R script under `File -> New File -> R Script` and type commands there, saving it from time to time. This way if you want to go back and repeat commands you can. To execute commands written in script, copy and paste commands to Console window and press Enter, press `Run` button in R-Studio or ask for a demo.

To use `DiffBind` package we need to install it first. To do so:
```bash

source("https://bioconductor.org/biocLite.R")
biocLite("DiffBind")

```

If the above worked, we should be able to load DiffBind library:
```bash

library(DiffBind)

```

### Running `DiffBind`  <a name="DB_run">

We will now follow `DiffBind` example to obtain differentially bound sites, given our samples. You may want to open `DiffBind` tutorial and read section [3 Example: Obtaining differentially bound sites](http://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) while typing the command to get more information about each step.

```bash

# reading in the sample information (metadata)
samples = read.csv("samples_REST.txt", sep="\t")

#	inspecting the metadata
samples

#	creating an object containing data
res=dba(sampleSheet=samples, config=data.frame(RunParallel=FALSE))

# inspecting the object: how many peaks are identified given the default settings?
res

# counting reads mapping to intervals (peaks)
# at this step a normalisation is applied by the default set to: score=DBA_SCORE_TMM_MINUS_FULL
res.cnt = dba.count(res, minOverlap=2, score=DBA_SCORE_TMM_MINUS_FULL, fragmentSize=130)

# inspecting the object: notice the FRiP values!
res.cnt

# plotting the correlation of libraries based on normalised counts of reads in peaks
pdf("correlation_libraries_normalised.pdf")
plot(res.cnt)
dev.off()

# PCA scores plot: data overview
pdf("PCA_normalised_libraries.pdf")
dba.plotPCA(res.cnt,DBA_TISSUE,label=DBA_TISSUE)
dev.off()

# setting the contrast
res.cnt2 = dba.contrast(res.cnt, categories=DBA_TISSUE, minMembers=2)

# inspecting the object: how many contrasts were set in the previous step
res.cnt2

# performing analysis of differential binding
res.cnt3 = dba.analyze(res.cnt2)

# inspecting the object: which condition are most alike, which are most different, is this in line with part one of the tutorial?
dba.show(res.cnt3, bContrasts = T)

# correlation heatmap  using only significantly differentially bound sites
# choose the contrast of interest e.g. HeLa vs. neuronal (#1)
pdf("correlation_HeLa_vs_neuronal.pdf")
plot(res.cnt3, contrast=1)
dev.off()

# boxplots to view how read distributions differ between classes of binding sites
# are reads distributed evenly between those that increase binding affinity HeLa vs. in neuronal?
pdf("Boxplot_HeLa_vs_neuronal.pdf")
pvals <- dba.plotBox(res.cnt3, contrast=1)
dev.off()

# extracting differentially binding sites in GRanges
res.db1 = dba.report(res.cnt3, contrast=1)
head(res.db1)

# plotting overlaps of sites bound by REST in different cell types
pdf("binding_site_overlap.pdf")
dba.plotVenn(res.cnt3, 1:4, label1="HeLa",label2="neuron",label3="HepG2",label4="sknsh")
dev.off()

# finally, let's save our R session including the generated data. We will need everything in the next section
save.image("diffBind.RData")
```

## Functional analysis <a name="FA">

So now we have list of differentially bound sites for comparisons of interest but we do not know much about them besides the genomic location. It is time to them in a biological context. To do so, we will use another `Bioconductor` package [ChIPpeakAnno](http://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html).

ChIPpeakAnno _"is for facilitating the downstream analysis for ChIP-seq experiments. It includes functions to find the nearest gene, exon, miRNA or custom features such as the most conserved elements and other transcription factor binding sites supplied by users, retrieve the sequences around the peak, obtain enriched Gene Ontology (GO) terms or pathways. Starting 2.0.5, new functions have been added for finding the peaks with bi-directional promoters with summary statistics (peaksNearBDP), for summarizing the occurrence of motifs in peaks (summarizePatternInPeaks) and for adding other IDs to annotated peaks or enrichedGO (addGeneIDs). Starting 3.4, permutation test has been added to determine whether there is a significant overlap between two sets of peaks. In addition, binding patterns of multiple transcription factors (TFs) or distributions of multiple epigenetic markers around genomic features could be visualized and compared easily using a side-by-side heatmap and density plot._

Here, we will annotate deferentially bound sites, summarise them in a genomic feature context and obtain enriched GO terms and pathways.


### Setting-up `ChIPpeakAnno`  <a name="FA_install">

We will continue our R-Studio session. If you have logged-out or lost connection or simply want to start fresh follow setting up instructions for running DiffBind locally.

To install ChIPpeakAnno
```bash

source("https://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")

```

We will also need to load DiffBind results saved in the differential binding session. We will build on them.
```bash

load("diffBind.RData")

```

### Running `ChIPpeakAnno`  <a name="FA_run">

Like with DiffBind package there is a nice [ChIPpeakAnno tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html#annotate-peaks) that you can view along this exercise to read more about the various steps.

```bash

# Loading DiffBind library
# we will need it to extract interesting peaks for down-stream analysis
library(DiffBind)

# Loading ChIPpeakAnno library
library(ChIPpeakAnno)

# Loading TSS Annotation For Human Sapiens (GRCh37) Obtained From BiomaRt
data(TSS.human.GRCh37)

# Choosing the peaks for the interesting comparison, e.g.
data.peaks = dba.report(res.cnt3, contrast=1)
head(data.peaks)

# Annotate peaks with information on closest TSS using precompiled annotation data
data.peaksAnno=annotatePeakInBatch(data.peaks, AnnotationData=TSS.human.GRCh37)

# View annotated peaks: can you see the added information in comparsition to data.peaks?
head(data.peaksAnno)

# Saving results
write.table(data.peaksAnno, file="peaks_HeLa_vs_neuronal.txt", sep="\t", row.names=F)
```


### Loading GO and REACTOME database <a name="FA_go_run">
Locally, we can install few more R libraries and annotation data to inspect our peaks a bit more. We will need libraries `org.Hs.eg.db`, `TxDb.Hsapiens.UCSC.hg19.knownGene` and `reactome.db`. To install:


```bash

source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

source("https://bioconductor.org/biocLite.R")
biocLite("reactome.db")

source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

```

### Enriched GO / REACTOME terms  <a name="FA_go_run">

```bash

library(org.Hs.eg.db)
library(reactome.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Peak distribution over genomic features
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peaks.featuresDist<-assignChromosomeRegion(data.peaksAnno, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs","Exons", "Introns"), TxDb=txdb)

pdf("peaks_featuresDistr_HeLa_vs_neuronal.pdf")
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(peaks.featuresDist$percentage, las=1, horiz=T)
dev.off()

# GO ontologies
peaks.go <- getEnrichedGO(data.peaksAnno, orgAnn="org.Hs.eg.db", maxP=.1, minGOterm=10, multiAdjMethod="BH", condense=TRUE)

# Preview GO ontologies results
head(peaks.go$bp[, 1:2])
head(peaks.go$mf[, 1:2])
head(peaks.go$cc[, 1:2])

# REACTOME pathways
peaks.pathways <- getEnrichedPATH(data.peaksAnno, "org.Hs.eg.db", "reactome.db", maxP=.05)

# REACTOME pathways: preview data
head(peaks.pathways)

# REACTOME pathways: list all pathways
print(unique(peaks.pathways$PATH))

```

Feel free to build more on the exercises. Follow the [ChIPpeakAnno tutorial](http://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html#annotate-peaks) for ideas.


## Concluding remarks and next steps <a name="Next">

The workflow presented in the tutorials is quite common and it includes recommended steps for analysis of ChIP-seq data. Naturally, there may be different tools or ways to preform similar tasks. New tools are being developed all the time and no single tool can do it all.

In the extra labs we have prepared you can find for instance an alternative way of quality control of ChIP-seq data with R package called `ChIPQC` as well as alternative differential binding workflow with a packaged called `csaw`. Note, these labs were not extensively tested so you may need to experiment and draw from the knowledge gained in the main labs.

Also, there are more types of analyses one can do beyond the one presented here. A common further analysis, for instance, includes identification of short sequence motifs enriched in regions bound by the assayed factor (peaks). There are several tools available here and we recommend you test one or two with on the tutorial data: [Homer](http://homer.salk.edu/homer/), [GEM](http://groups.csail.mit.edu/cgs/gem/), [RSAT](http://floresta.eead.csic.es/rsat/peak-motifs_form.cgi)m [MEME](http://meme-suite.org/)

Above all, we recommend that you keep trying to analyze your own data. Practice makes perfect :)

----

## Appendix: figures <a name="Appendix">

![covplot](../figures/lab-chipseeker/fig-covplot-1.pdf)

Fig: Coverage plot for HeLa replicate 1 peaks

----

![covplot-cmp](../figures/lab-chipseeker/fig-covplot-cmp-1.pdf)

Fig: Coverage plot for HeLa replicate 1 and 2 peaks

----
