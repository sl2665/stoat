#!/bin/bash
# Find elongation rates from the clearance wave data
# as in 2014 eLife paper  DOI: 10.7554/eLife.02407

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat elongHMM"
	echo -e "version: 0.1.191114"
	echo -e ""
	echo -e "usage:   stoat elongHMM [options] -f <Clearance PRO-seq dir> -b <Baseline PRO-seq dir>"
	echo -e ""
	echo -e "options:"
	echo -e "         -bs  Bin size (default = 1000)"
	echo -e "         -bc  Bin count (default = 50)"
	echo -e "         -out Output file (default = out.elong.txt)"
	echo -e ""
	exit
fi

BS=1000
BC=50
OUT="out.elong.txt"

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-f) F="$2"; shift ;;
		-b) B="$2"; shift ;;
		--sdir) SDIR="$2"; shift;;
		--default) ;;
		*) ;;
	esac
	shift
done

${SDIR}/hmm2 -p ${F}/alignment/a.pl.bedgraph -m ${F}/alignment/a.mn.bedgraph \
	-p0 ${B}/alignment/a.pl.bedgraph -m0 ${B}/alignment/a.mn.bedgraph \
	-g ${B}/annotation/transcripts.bed13 -bs ${BS} -bc ${BC} \
	> ${OUT}
