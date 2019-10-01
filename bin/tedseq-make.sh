#!/bin/bash
# process TED-seq files to the data directory structure

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat make-ted"
	echo -e "version: 0.1.190924"
	echo -e ""
	echo -e "usage:   stoat make-ted [options] -f <fastq> -g <gtf> -r <reference genome>"
	echo -e ""
	echo -e "options:"
	echo -e "         -a   aligner (STAR/BOWTIE; default = STAR)"
	echo -e "         -o   output directory (default = tedseq.out)"
	echo -e ""
	exit
fi

ODIR="tedseq.out"
ALIGNER="STAR"

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-f)
		FASTQ="$2"
		shift
		;;
		-g)
		GTF="$2"
		shift
		;;
		-r)
		REF="$2"
		shift
		;;
		-a)
		ALIGNER="$2"
		shift
		;;
		-o)
		ODIR="$2"
		shift
		;;
		--sdir)
		SDIR="$2"
		shift
		;;
		--default)
		;;
		*)
		;;
	esac
	shift
done

echo generating output directory >&2
# Generate output directory structure
if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi
if [ ! -d $ODIR/alignment ]; then
	mkdir $ODIR/alignment
fi
if [ ! -d $ODIR/annotation ]; then
	mkdir $ODIR/annotation
fi
if [ ! -d $ODIR/table ]; then
	mkdir $ODIR/table
fi
if [ ! -d $ODIR/Rdata ]; then
	mkdir $ODIR/Rdata
fi
if [ ! -d $ODIR/browser ]; then
	mkdir $ODIR/browser
fi
if [ ! -d $ODIR/_tmp ]; then
	mkdir $ODIR/_tmp
fi

if false; then
# Align fastq files to bam and bedgraph files
echo performing alignment >&2
$SDIR/tedseq-align -f $FASTQ -r $REF -a $ALIGNER -b $ODIR/alignment/a --sdir $SDIR
# Generate annotation
echo generating gene annotations >&2
$SDIR/tedseq-make-annotation.sh $GTF $ODIR
# Make expression table
echo calculating gene expression levels >&2
$SDIR/tedseq-getexpr -t $ODIR/alignment/a -g $ODIR/annotation/transcripts.bed13 \
	--sdir $SDIR > $ODIR/table/expression.txt
# Make PAL matrix
echo generating polyA length matrix >&2
$SDIR/tedseq-makepal -a $ODIR/alignment/a.bam -b $ODIR/annotation/transcripts.bed13 \
	--sdir $SDIR > $ODIR/table/palmatrix.txt
# Filter PAL matrix to expressed genes and calculate median poly A lengths
# Covert PAL matrix to Rdata format
$SDIR/process_palmatrix.sh $ODIR 200 250
fi
# Generate 3'CPS region
echo annotating 3\'CPS regions >&2
$SDIR/tedseq-find3cps 
# Make transcriptome-wide heatmap
