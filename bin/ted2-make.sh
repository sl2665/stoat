#!/bin/bash
# process TED2 bam files to the data directory structure

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat make-ted2"
	echo -e "version: 0.1.190924"
	echo -e ""
	echo -e "usage:   stoat make-ted2 [options] -b <bam> -g <gtf>"
	echo -e ""
	echo -e "options:"
	echo -e "         -o   output directory (default = tedseq.out)"
	echo -e "         --s  library size (default = 425 bp)"
	echo -e ""
	exit
fi

ODIR="tedseq.out"
SIZE=425
while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-b)
		BAM="$2"
		shift
		;;
		-g)
		GTF="$2"
		shift
		;;
		-o)
		ODIR="$2"
		shift
		;;
		--s)
		SIZE="$2"
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

# Copy bam file
cp $BAM $ODIR/alignment/a.bam 

# Make bedgraph files
bedtools genomecov -ibam $BAM -bg -strand + -5 \
	> $ODIR/alignment/a.pl.bedgraph
bedtools genomecov -ibam $BAM -bg -strand - -5 \
	| awk '{print $1"\t"$2"\t"$3"\t"$4*-1}' \
	> $ODIR/alignment/a.mn.bedgraph

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
	--sdir $SDIR > $ODIR/table/palmatrix_raw.txt

echo processing polyA length matrix >&2
# Process PAL matrix to remove outlier peaks
RDIR=${SDIR%*/bin}/rscript
Rscript --vanilla --quiet $RDIR/ted2_make.R ${ODIR}

# Filter PAL matrix to expressed genes and calculate median poly A lengths
# Covert PAL matrix to Rdata format
SPOS=$(( 500 - $SIZE + 136 ))
$SDIR/process_palmatrix.sh $ODIR $SPOS 250

# Generate 3'CPS region
# echo annotating 3\'CPS regions >&2
# $SDIR/tedseq-find3cps 
# Make transcriptome-wide heatmap
