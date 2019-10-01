#!/bin/bash
# Redefine new 3'CPS from TED-seq junction data files

if [ $# -lt 2 ]; then
	echo -e "Usage:\tnew3cps.sh [options] -d <TED-seq dir>"
	echo -e "Options:"
	echo -e "\t-l\tdirectory to TED-seq junction library (default = HEK)"
	exit
fi

A=10

if [ ! -d _tmp ]; then
	mkdir _tmp
fi

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-d)
		TED="$2"
		shift
		;;
		-l)
		LIB="$2"
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

# Use tedseq-find3cps to find approximate 3'CPS sites from TED-seq read clusters
tedseq-find3cps -t $TED/alignment/a > _tmp/cpsregion.bed

# Split 3'CPS sites between plus and minus strands, and take 200 bp region from the 3'ends
awk '$6=="+"{print $1"\t"$3-100"\t"$3+100 > "_tmp/pl.bed";next}
	$2>100{print $1"\t"$2-100"\t"$2+100 > "_tmp/mn.bed"}' _tmp/cpsregion.bed

# Extract TED-seq read peaks that mapped to the 3'CPS clusters from the 100bp short insert TED-seq library
# Get read counts for all positions within blocks
bedtools map -a _tmp/pl.bed -b data/3CPS.HEK/alignment/a.pl.bedgraph -o collapse,collapse -c 2,4 \
	> _tmp/pl.col
bedtools map -a _tmp/mn.bed -b data/3CPS.HEK/alignment/a.mn.bedgraph -o collapse,collapse -c 2,4 \
	> _tmp/mn.col

# Find the position with the maximum read count within blocks
awk '{n = split($4, pos, ","); split($5,readcount,","); \
	maxpos = 1; for(i = 1; i<=n; ++i) if(readcount[ i ]>readcount[ maxpos ]) maxpos=i; \
	print $1"\t"pos[ maxpos ]-300"\t"pos[ maxpos ]"\t"$1":"pos[ maxpos ]"\t"readcount[ maxpos ]"\t+"}' \
	_tmp/pl.col > _tmp/pl.max

awk '{n = split($4, pos, ","); split($5,readcount,","); \
	maxpos = 1; for(i = 1; i<=n; ++i) if(readcount[ i ]>readcount[ maxpos ]) maxpos=i; \
	print $1"\t"pos[ maxpos ]"\t"pos[ maxpos ]+300"\t"$1":"pos[ maxpos ]"\t"readcount[ maxpos ]"\t-"}' \
	_tmp/mn.col > _tmp/mn.max

# Merge the plus and minus strand max peak results into bed files
cat _tmp/pl.max _tmp/mn.max | \
	awk '$2!="."' | sort -k1,1 -k2,2n > $TED/annotation/3cps.ref.bed

