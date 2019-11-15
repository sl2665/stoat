# stoat
**St**ere**o**scopic **a**nalysis of the **t**ranscriptome using PRO-seq and TED-seq

## Summary <img src="img/stoat.png" alt="drawing" width="100" align="left"/>
Stereoscopic transcriptome analysis that adds dynamics-perception using binocular perspectives of the static transcriptome; transcription rate and polyadenylation status. Transcription rate, measured by nascent RNA sequencing (PRO-seq), reflects the RNA synthesis level. Poly(A) tail length (PAL), measured by TED-seq, reflects the quality of the RNA, and is associated with translation and decay status. The starts and the ends of transcripts can be mapped using PRO-seq and TED-seq to re-define high-confidence annotation including novel transcripts. Transcriptional and post-transcriptional expression analysis identifies dynamically regulated transcripts.

## Installation

### Supported OS
Linux, Mac OS X

### Prerequisites
* samtools (http://www.htslib.org/)
* bedtools (https://bedtools.readthedocs.io/en/latest/)
* cutadapt (https://cutadapt.readthedocs.io/en/stable/)
* STAR or bowtie aligner (https://github.com/alexdobin/STAR/releases or http://bowtie-bio.sourceforge.net/index.shtml)
* dREG (https://github.com/Danko-Lab/dREG)
* bedGraphToBigWig, gtfToGenePred, genePredToBed (http://hgdownload.soe.ucsc.edu/admin/exe/)

These softwares should be accessibile in your PATH

### Install stoat
* Download to your installation directory
* Go to your installation directory and run make
```
cd /(your)/(installation)/(directory)
make
```
* Copy the stoat file to one of your PATH accessible directory
```
cp /(your)/(installation)/(directory)/bin/stoat ~/.local/bin/ 
```

## Flowchart
<img src="img/stoat.jpg" alt="drawing" width="800"/>

## Quickstart
In the example directory,
```
cd example
gunzip *.gz
```
PRO-seq pipeline
```
proseq-align -f PROseq.chr22.fastq -r (human reference genome)
proseq-make-dREG -p proseq.out -s (dREG SVM RData) > proseq.dREG.bed
proseq-hmm -p proseq.out -mp mappability.chr22.bedgraph > proseq.hmm.bed
proseq-getexpr -p proseq.out -g GencodeComprehensiveV26-hg38.chr22.bed > proseq.expr.txt
```
TED-seq pipeline
```
tedseq-align -f TEDseq.chr22.fastq -r (human reference genome)
tedseq-find3cps -t tedseq.out > tedseq.cps.bed
tedseq-makepal -t tedseq.out.bam -g GencodeComprehensiveV26-hg38.chr22.bed > tedseq.pal.txt
tedseq-getexpr -t tedseq.out -g GencodeComprehensiveV26-hg38.chr22.bed > tedseq.expr.txt
```
Combined analysis pipeline
```
stoat-getannot -a proseq.dREG.bed -b proseq.hmm.bed -c tedseq.3cps.bed -g GencodeComprehensiveV26-hg38.chr22.bed > stoat.annot.txt
stoat-getdge -pro proseq.expr.txt -ted tedseq.expr.txt -pal tedseq.pal.txt -ins 300 > stoat.dge.txt
```

## Usage

### proseq-align
```
Usage:   proseq-align [options] -f <fastq> -r <reference genome>
Options:
        -a      alinger (STAR/BOWTIE; default = STAR)
        -b      output filename base (default = proseq.out)
```
Output: generates 3 files
  * \<output filename base>.bam : aligned bam file with unique molecular identifiers collapsed
  * \<output filename base>.pl.bedgraph : (+) strand bedgraph file of PRO-seq raw read counts
  * \<output filename base>.mn.bedgraph : (-) strand bedgraph file of PRO-seq raw read counts
         

## Documentation

## How to cite
Lee S.A., Kwak H. Stereoscopic transcriptome analysis depicts transcriptional and post-transcriptional RNA regulation
