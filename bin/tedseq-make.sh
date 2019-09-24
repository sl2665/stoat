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
	echo -e "         -d   output directory (default = tedseq.out)"
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
		-r)
		REF="$2"
		shift
		;;
		-a)
		ALIGNER="$2"
		shift
		;;
		-d)
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

# Align fastq files to bam and bedgraph files
$SDIR/tedseq-align -f $FASTQ -r $REF -a $ALIGNER -b $ODIR/alignment/a --sdir $SDIR

# Generate annotation  
