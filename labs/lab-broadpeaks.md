# Detection of broad peaks from ChIP-seq data


## Requirements

* MACS 2.1.2 (this version is not available as a module on Rackham)


**MACS 2.1.2 installation (on Uppmax):**

in your home folder on Rackham:

```
## clone pyenv repository from github
git clone git://github.com/yyuu/pyenv.git ~/.pyenv

## make pyenv start at a login
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile
echo 'eval "$(pyenv init -)"' >> ~/.bash_profile

source ~/.bash_profile

## install required version of python and set it to use
pyenv install 2.7.9
pyenv global 2.7.9

## create directory for macs2
mkdir macs2
cd macs2

## clone macs repository from github
git clone https://github.com/taoliu/MACS.git

## install MACS2 with its dependencies
pip install MACS2

## export the location of /macs/MACS/bin folder to $PATH
## to obtain the path you can use `pwd`

in my case:
export PATH=$PATH:/home/agata/soft/macs/MACS/bin

```

The detailed instructions on how to install and use pyenv on Uppmax are at

https://www.uppmax.uu.se/support/user-guides/python-modules-guide/


## Data

We will use ChIP-seq of H3K79me2 from Orlando et al, 2014 ("Quantitative ChIP-Seq Normalization Reveals Global Modulation of the Epigenome"). H3K79me2 is enriched at active promoters and linked to transcriptional activation. This is a SE data set, which admitedly is not the best design for broad marks. To use this procedure with PE data, please follow modifications listed on https://github.com/taoliu/MACS.


GEO accession is `GSE60104`
ENA accession is ` PRJNA257491`

files in the dataset:

sample | GEO accession | SRA accession
--- | --- | --- 
Jurkat_K79_100%_R1 | GSM1465008 | SRR1536561
Jurkat_K79_100%_R2 | GSM1464998 | SRR1536551
Jurkat_K79_50%_R1 | GSM1465006 | SRR1536559
Jurkat_K79_0%_R1 | GSM1465004 | SRR1536557
Jurkat_WCE_100%_R1 | GSM1511469 | SRR1584493
Jurkat_WCE_100%_R2 | GSM1511474 |SRR1584498
Jurkat_WCE_50%_R1 | GSM1511467 | SRR1584491
Jurkat_WCE_0%_R1 | GSM1511465 | SRR1584489


We will call peaks from one sample only, and compare the results to other samples processed earlier.


Data have been processed in the same way as for the TF ChIP-seq, i.e. the duplicated reads were removed, as were the reads mapped to blacklisted regions. In this case the reads were mapped to `hg38` assembly of human genome.


## Quality control

As always, one should start the analysis from assesment of data quality. This is already performed, and the plots and metrics are below.

### Cross-correlation and related metrics

The files discussed in this section can be accessed at 
`/sw/share/compstore/courses/ngsintro/chipseq/broad_peaks/results_pre/fingerprint`
and
`/sw/share/compstore/courses/ngsintro/chipseq/broad_peaks/results_pre/xcor`
.

These metrics have been developed with application to TF ChIP-seq in mind, and you can see that the results for broad peaks are not as easy to interpret as for point-source factors. Below are cross correlation plots for the IP and input you are going to use for the exercise. Already from these plots alone it is evident that the data has some quality issues. At this point you should be able to identify them.


ChIP:

<img src="../figures/lab-broadpeaks/SRR1536557_xcor.png" alt="" style="width: 100px;"/><br>

input:

<img src="../figures/lab-broadpeaks/SRR1584489_xcor.png" alt="" style="width: 100px;"/><br>


As for the ChIP, the cross correlation profile of factors with broad occupancy patterns is not going to be as sharp as for TFs, and the values of NSC and RSC tend to be lower, which does not mean that the ChIP failed. In fact, the developers of the tool do not recommend using the same NSC / RSC values as quality cutoffs for broad marks. However, input samples should not display signs of enrichment, as is the case here.

### Cumulative enrichment

Another plot worth examining is cumulative enrichment (aka fingerprint from deepTools):

<img src="../figures/lab-broadpeaks/cmplGSE60104fingerprint.png" alt="" style="width: 100px;"/><br>

You can see that even though the cross correlation metrics don't look great, to put it mildly, some enrichment can be observed for the ChIP samples, and not for the input samples. As this data is data from very shallow sequencing, the fraction of the genome covered by reads is smaller than expected (0.3 for the best sample). Thus we do not expect to detect all occupancy sites, only the ones which give the strongest signal (this is actually an advantage for this class, as it reduces the running time).


## Peak calling

You will call peaks using sample Jurkat_K79_50_R1 (`SRR1536557`) and its matching input `SRR1584489`.
Effective genome size for hg38 is `3.0e9`.
The estimated fragment size is `180 bps` (`phantompeakqualtools`).

```
mkdir -p results/macs
cd results/macs

ln -s /sw/share/compstore/courses/ngsintro/chipseq/broad_peaks/bam/SRR1536557.bwt.hg38_dm6.sorted.hg38.BLfilt.bam
ln -s /sw/share/compstore/courses/ngsintro/chipseq/broad_peaks/bam/SRR1584489.bwt.hg38_dm6.sorted.hg38.BLfilt.bam

#if it is a different session than when installing pyenv and macs2:
pyenv global 2.7.9

macs2 callpeak -t SRR1536557.bwt.hg38_dm6.sorted.hg38.BLfilt.bam -c SRR1584489.bwt.hg38_dm6.sorted.hg38.BLfilt.bam -n 50_R1 --outdir 50_R1 -f BAM --gsize 3.0e9 -q 0.1 --nomodel --extsize 180 --broad --broad-cutoff 0.1

```

You can now inspect the results in the output folder `50_R1`. The structure is alike the output for calling narrow peaks. The file `*.broadPeak` is in `BED6+3` format which is similar to `narrowPeak` file used for point-source factors, except for missing the 10th column for annotating peak summits. Look here (https://github.com/taoliu/MACS) for details.

How many peaks were identified?

```
[agata@r483 50_R1]$ wc -l *Peak
  46664 50_R1_peaks.broadPeak
```

This is a preliminary peak list, and in case of broad peaks, it almost always needs some processing or filtering.

## Visual inspection of the peaks

You will use IGV for this step, and it is recommended that you run it locally on your own computer. Please load `hg38` reference genome.

Required files are:

* SRR1536557.bwt.hg38_dm6.sorted.hg38.BLfilt.bam and bai
* SRR1584489.bwt.hg38_dm6.sorted.hg38.BLfilt.bam and bai
* 50_r1_peaks.broadPeak

You can access the bam and bai files from
`/sw/share/compstore/courses/ngsintro/chipseq/broad_peaks/bam/`.

You can look at the locations of interest. Some peaks with low FDR (q value) or high fold enrichment may be worth checking out. Or check your favourite gene.

Some ideas:

```
chr1:230,145,433-230,171,784
chr1:235,283,256-235,296,431
chr1:244,857,626-244,864,213
chr1:45,664,079-45,690,431

chr1:45,664,079-45,690,431
```

The first two locations visualise peaks longer than 2kb. The third and the fourth are a 4 kb-long peaks with fold erichment over background >15.


An example (two upper tracks are ChIP samples, the bottom track is input; the annotation is refseq genes and peaks called for sample 100_r1):

<img src="../figures/lab-broadpeaks/broad3.png" alt="" style="width: 400px;"/><br>


All the above but, perhaps the fifth location most of all, demonstrate one of the common caveats of calling broad peaks: regions obviously enriched in a mark of interest are represented as a series of adjoining peaks which in fact should be merged into one long enrichment domain. You may leave it as is, or merge the peaks into longer ones, depending on the downstream application.

## Postprocessing of peak candidates

Please note that this step is only an example, as ***any postprocessing of peak calling results is highly project specific***.

Normally, you would work with replicated data. As in the case of TFs earlier, it is recommended to continue working with peaks reproducible between replicates.

The peak candidate lists can and should be further filtered, based on fold enrichment and pileup value, to remove peaks which could have a high fold enrichment but low signal, as these are likely non-informative. Any filtering, however has to be performed having in mind the biological characteristics of the signal.

You can merge peaks which are close to one another using bedtools (https://bedtools.readthedocs.io/en/latest/). You will control the distance of features to be merged using option `-d`. Here we arbitrarily choose 1 kb.

```
cp 50_r1_peaks.broadPeak 50_r1.bed

module load bioinfo-tools
module load BEDTools/2.27.1

bedtools merge -d 1000 -i 50_r1.bed > 50_r1.merged.bed

#how many peaks?
wc -l 50_r1.merged.bed 

#11732 50_r1.merged.bed
```
