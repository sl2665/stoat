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
cd /(your)/(installation)/(directory)/bin/stoat/example
gunzip *.gz
```
PRO-seq pipeline
```
stoat make-pro -f PROseq.chr22.fastq -g gencode.v26.annotation.chr22.gtf -r (human reference genome)
```
TED-seq pipeline
```
stoat make-ted -f TEDseq.chr22.fastq -g gencode.v26.annotation.chr22.gtf -r (human reference genome)
```
Combined analysis pipeline
```

```

## Usage

### stoat
```
stoat:   stereoscopic analysis of the transcriptome using PRO-seq and TED-seq
version: 0.1.190924

usage:   stoat <command> [options]

command: make-pro   process PRO-seq data
         make-ted   process TED-seq data
         redef3     redefine 3CPS in high resolution
         pap        generate polyA length profiles
         prop       generate PRO-seq profiles
         elongHMM   calculate PRO-seq elongation rates
```

### make-pro
```
tool:    stoat make-pro
version: 0.1.191014

usage:   stoat make-pro [options] -f <fastq> -g <gtf> -r <reference genome>

options:
         -a   aligner (STAR/BOWTIE; default = STAR)
         -o   output directory (default = proseq.out)
```
Output: generates a directory structure
  (output directory)/
  ├── alignment
  │   ├── a.bam
  │   ├── a.mn.bedgraph
  │   └── a.pl.bedgraph
  ├── annotation
  │   └── transcripts.bed13
  └── sample_info.txt
  * /alignment/a.bam : aligned bam file with unique molecular identifiers collapsed
  * /alignment/a.pl.bedgraph : (+) strand bedgraph file of PRO-seq raw read counts
  * /alignment/a.mn.bedgraph : (-) strand bedgraph file of PRO-seq raw read counts
  

## Documentation

## How to cite
Lee S.A., Kwak H. Stereoscopic transcriptome analysis depicts transcriptional and post-transcriptional RNA regulation
