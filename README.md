# stoat
**St**ere**o**ptic **a**nalysis of the **t**ranscriptome using PRO-seq and TED-seq

## Summary <img src="img/stoat.png" alt="drawing" width="100" align="left"/>
Stereoptic transcriptome analysis that adds dynamics-perception using binocular perspectives of the static transcriptome; transcription rate and polyadenylation status. Transcription rate, measured by nascent RNA sequencing (PRO-seq), reflects the RNA synthesis level. Poly(A) tail length (PAL), measured by TED-seq, reflects the quality of the RNA, and is associated with translation and decay status. The starts and the ends of transcripts can be mapped using PRO-seq and TED-seq to re-define high-confidence annotation including novel transcripts. Transcriptional and post-transcriptional expression analysis identifies dynamically regulated transcripts.

## Installation

### Supported OS
Linux, Mac OS X

### Prerequisite
* Samtools (http://www.htslib.org/)
* Bedtools (https://bedtools.readthedocs.io/en/latest/)
* STAR or bowtie aligner (https://github.com/alexdobin/STAR/releases or http://bowtie-bio.sourceforge.net/index.shtml)
* dREG (https://github.com/Danko-Lab/dREG)

## Flowchart
<img src="img/STOAT-FLOWCHART.png" alt="drawing" width="800" />

## Usage

### proseq-align
```
Usage:   proseq-align [options] -f <fastq> -r <reference genome>
Options:
        -a      alinger (STAR/BOWTIE; default = STAR)
        -b      output filename base (default = proseq.out)
```
Output:
  * (output filename base).bam : aligned bam file with unique molecular identifiers collapsed
  * (output filename base).pl.bedgraph : (+) strand bedgraph file of PRO-seq raw read counts
  * (output filename base).mn.bedgraph : (-) strand bedgraph file of PRO-seq raw read counts
         
### proseq-makedREG
```
Usage:   proseq-makedREG [options] -p <PRO-seq plus bedgraph> -m <PRO-seq minus bedgraph>
Options:
        -r      dREG options
```
Output:  [output filename base].bed
### proseq-HMM
```
Usage:  proseq-HMM      [options] -p <PRO-seq plus bedgraph> -m <PRO-seq minus bedgraph>
Options:
        -r      ChromHMM options 
```
Output: [TAR].bed 
### proseq-getexpr
```
Usage:  proseq-makedREG [options] -p <PRO-seq plus bedgraph> -m <PRO-seq minus bedgraph>
Options: 
        -r      Hg19 Promoters Options
```
Output: [Transcript activity].bed 
### tedseq-align
```
Usage:   tedseq-align [options] -f <fastq> -r <reference genome>
Options:
        -a      alinger (STAR/BOWTIE; default = STAR)
        -b      output filename base (default = tedseq.out)
```
Output:
  * (output filename base).bam : aligned bam file with unique molecular identifiers collapsed
  * (output filename base).pl.bedgraph : (+) strand bedgraph file of TED-seq raw read counts
  * (output filename base).mn.bedgraph : (-) strand bedgraph file of TED-seq raw read counts
```
### tedseq-find3cps   

        Usage: tedseq-find3cps [options] -p <TED-seq plus bedgraph> -m <TED-seq minus bedgraph>
        Options:
                -r      3CPS options
 ```
Output:[tedseq.3CPS].bed
 
### tedseq-makepal 
```
        Usage: tedseq-makepal [options] -a <bam> -b <bed>
        Options:
        -bin    bin size (default = 1)
        -win    window size (default = 500)
```
Output:[palmatrix.3CPS].bed
### tedseq-getexpr
```
        Usage:  TEDseq.bam [options] -p <TED-seq plus bedgraph> -m <TED-seq minus bedgraph>
        Options: 
                -r      Hg19 Promoters Options
```
Output: [tedseq.expressed].bed
### stoat-getannot
```
        Usage:  tedseq.3CPS.bed[options] -r <reference genome>
        Options:
                -r      Hg19 Reference genome  
```
Output:[STOAT.annotated].bed 
### stoat-getdge
```
        Usage: Transcript activity.bed, RNA-exp.bed, palmatrix.3CPS.bed 
        Options:
                -r      
```
Output:[STOAT.expressed].bed
## Documentation

## How to cite
Lee S.H., Woo Y.M., Kwak H. (2018). Stereoptic transcriptome analysis refines functional gene annotation and identifies polyadenylation of atypical transcripts.
