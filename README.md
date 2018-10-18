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
Output: generates 3 files
  * \<output filename base>.bam : aligned bam file with unique molecular identifiers collapsed
  * \<output filename base>.pl.bedgraph : (+) strand bedgraph file of PRO-seq raw read counts
  * \<output filename base>.mn.bedgraph : (-) strand bedgraph file of PRO-seq raw read counts
         
### proseq-makedREG
```
Usage:   proseq-makedREG [options] -p <PRO-seq filename base> -s <dREG SVM RData>
Options:
        -c      number of CPU cores (default = 4)
        -g      GPU ID if multiple GPUs are available(default = NA)
```
Output: bedgraph format of dREG scores (transcriptional activity)

### proseq-HMM
```
Usage:  proseq-HMM      [options] -p <PRO-seq filename base>
Options:
        -w      window size (default = 500 bp) 
```
Output: bed format of transcription active regions (TAR)
1. chromosome
2. start
3. end
4. TAR ID
5. (empty)
6. strand

### proseq-getexpr
```
Usage:  proseq-getexpr [options] -b <PRO-seq filename base> -g <gene annotation bed12>
Options: 
        -w      Promoter range (default = 500 bp)
        -wu     Promoter upstream range (default = 500 bp)
        -wd     Promoter downstream range (default = 500 bp)
```
Output: reports
1. Gene ID/name in annotated bed file
2. Promoter raw read counts
3. Gene body read coverage
4. Exon read coverage

### tedseq-align
```
Usage:   tedseq-align [options] -f <fastq> -r <reference genome>
Options:
        -a      alinger (STAR/BOWTIE; default = STAR)
        -b      output filename base (default = tedseq.out)
```
Output:
  * \<output filename base>.bam : aligned bam file with unique molecular identifiers collapsed
  * \<output filename base>.pl.bedgraph : (+) strand bedgraph file of TED-seq raw read counts
  * \<output filename base>.mn.bedgraph : (-) strand bedgraph file of TED-seq raw read counts

### tedseq-find3cps   
```
Usage:  tedseq-find3cps [options] -b <TED-seq filename base>
Options:
        -c      Raw read count cut-off (default = 5)
 ```
Output: bed format of 3\` cleavage poly-adenylation sites
1. chromosome
2. start
3. end
4. 3\`CPS ID
5. read counts supporting the 3\`CPS
6. strand
 
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
