#!/bin/bash
# calibrate TED2 bam files to the data directory structure

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    ted2-calibrate.sh"
	echo -e "version: 0.1.220315"
	echo -e ""
	echo -e "usage:   ted2-calibrate.sh [options] -g <genome.fa> -n <conversion.net>"
	echo -e ""
	echo -e "options:"
  echo -e "         -t   tedseq directory (default = tedseq.out)"
	echo -e ""
	exit
fi

TD="tedseq.out"

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-t)
		TD="$2"
		shift
		;;
		-g)
		FA="$2"
		shift
		;;
		-n)
		NET="$2"
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

# Read PAL table

# Make UTR sequence fasta of PAL table coordinates, separate strands

# Construct input sequence lists for all PAL table coordinates with >= 1 reads

# Run neural network on sequences to get conversion factor

# Multiply conversion factors to PAL table

