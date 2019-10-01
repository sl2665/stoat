#!/bin/bash
# Generate polyA profile 

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat pap - poly(A) profile"
	echo -e "version: 0.1.190928"
	echo -e ""
	echo -e "usage:   stoat pap [options] -d <TEDseq dirs> -g <gene lists txt>"
	echo -e ""
	echo -e "options:"
	echo -e "         -b   3'CPS bed files instead of gene list (default = none)"
	echo -e "         -o   output plot file name (default = palout.pdf)"
	echo -e "         -w   width of the output pdf (default = 6)"
	echo -e "         -h   height of the output pdf (default = 4)"
	echo -e "         --od output data (default = none)"
	echo -e "         --de sample/gene group descriptions (default = none)"
	echo -e "         --os overlay TED-seq data (default = false)"
	echo -e "         --og overlay gene groups (default = false)"
	echo -e ""
	exit
fi

BED=false
OUT=false
GN=false
DSC=false
OF="palout.pdf"
WIDTH=5
HEIGHT=4
OS=false
OG=false

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-d) TED=("$2"); while [[ ${3:0:1} != "-" ]]; do TED+=("$3"); shift; done; shift ;;
		-g) GN="$2"; shift ;;
		-o) OF="$2"; shift ;;
		-w) WIDTH="$2"; shift ;;
		-h) HEIGHT="$2"; shift ;;
		-b) BED=("$2"); while [[ ${3:0:1} != "-" ]]; do BED+=("$3"); shift; done; shift ;;
		--od) OUT="$2"; shift ;;
		--de) DSC=("$2"); while [[ ${3:0:1} != "-" ]]; do DSC+=("$3"); shift; done; shift ;;
		--sdir) SDIR="$2"; shift ;;
		--os) OS=true ;;
		--og) OG=true ;;
		--default) ;;
		*) ;;
	esac
	shift
done

# Directories for plotting poly(A) profiles
RDIR=${SDIR%*/bin}/rscript
if [ ! -d _tmp ]; then mkdir _tmp; fi
if [ -f _tmp/bed ]; then rm _tmp/bed; fi
if [ -f _tmp/pl ]; then rm _tmp/pl; fi
if [ -f _tmp/mn ]; then rm _tmp/mn; fi

# For producing PAL from new bed files (fast assuming no soft clipping)
if [ "$BED" != false ] ; then
	if [ "$DSC" == false ] ; then DSC=("${BED[@]}"); fi
	# process multiple bed files to one adding index
	for ((i=0;i<${#BED[@]};++i)); do
		awk -v i=${DSC[i]} '{print $0"\t"i}' ${BED[i]} >> _tmp/bed
	done
	# extend bed file to -500 from the CPS to +0
	awk '$6=="+"{if($3>500) print $1"\t"$3-500"\t"$3"\t"$4";"$7"\t"$5"\t+">"_tmp/pl";next}
		{print $1"\t"$2"\t"$2+500"\t"$4";"$7"\t"$5"\t-">"_tmp/mn"}' _tmp/bed
	
	# Generate poly(A) matrix
	if [ -f _tmp/pl ]; then
		sort -k1,1 -k2,2n _tmp/pl > _tmp/pls
		$SDIR/bgmatrix -a _tmp/pls -b ${TED[0]}/alignment/a.pl.bedgraph -bin 1 > _tmp/plmat
	else
		echo '' > _tmp/pls
		echo '' > _tmp/plmat
	fi
	if [ -f _tmp/mn ]; then
		sort -k1,1 -k2,2n _tmp/mn > _tmp/mns
		$SDIR/bgmatrix -a _tmp/mns -b ${TED[0]}/alignment/a.mn.bedgraph -bin 1 | \
			awk '{printf int(-1)*$NF; for(i=NF-1;i>=1;--i) printf "\t"int(-1)*$i; printf "\n"}' > _tmp/mnmat
	else
		echo '' > _tmp/mns
		echo '' > _tmp/mnmat
	fi
	# Merge the two strands
	cat <(paste <(cut -f1-4 _tmp/pls) <(cat _tmp/plmat)) <(paste <(cut -f1-4 _tmp/mns) <(cat _tmp/mnmat)) | \
		sort -k1,1 -k2,2n -k3,3n | cut -f4- > _tmp/palmat
	# Run R script to read palmat and ouput results
	Rscript --vanilla --quiet $RDIR/palmat2profile.R \
		${OF} ${WIDTH} ${HEIGHT} &> _tmp/Routput
	rm -rf _tmp
	exit
else
	if [ "$GN" != false ] ; then
	# Producing PAL from gene lists (support multiple TED-seq samples)
		# Make TED-seq sample list
		if [ "$DSC" == false ] ; then DSC=("${TED[@]}"); fi
		for ((i=0;i<${#TED[@]};++i)); do
			echo ${TED[i]},${DSC[i]} >> _tmp/sample; done
		# Generate gene id list from gene annotation tables
		for ((i=0;i<${#TED[@]};++i)); do
			cut -f1 ${TED[i]}/table/palmatrix.expr.txt >> _tmp/expGenes
		done
		sort -k1,1 -u _tmp/expGenes > _tmp/expGenesUniq
		awk '$6=="+"{print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t+\t"$13;next}
		{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t-\t"$13}' \
			${TED[0]}/annotation/transcripts.bed13 | sort -k1,1 -k2,2n > _tmp/geneAnnot
		# Filter for the expressed 3'CPS genes
		awk 'NR==FNR{id[$1]=1;next}
			id[$4]{print $0}' _tmp/expGenesUniq _tmp/geneAnnot > \
			_tmp/expAnnot
		sort -k1,1 -k2,2n -u _tmp/expAnnot > _tmp/uniqAnnot
		# Filter out any 3' annotation loacted within 300 bp of each other
		awk 'NR==1{prev=$2;fl[NR]=$0;next}
				{cur=$2;if(cur-prev<300) { fl[NR-1]=0; fl[NR]=0} else fl[NR]=$0; prev=cur}
				END{for(i in fl) if(fl[i]) print fl[i]}' _tmp/uniqAnnot > _tmp/sepAnnot
		# Find gene ids in the separated annotation
		awk 'NR==FNR{pos[$1":"$2]=$0;next}
			pos[$1":"$2]{print pos[$1":"$2]}' _tmp/geneAnnot _tmp/sepAnnot | \
			cut -f4,7 > _tmp/annot
		# to plot individual genes without grouping, do not filter out any annotations
		awk 'NR==FNR{if(NF==1) p=1;next}
			p{print $4"\t"$7 > "_tmp/annot"}' ${GN} _tmp/geneAnnot
		# Run Rscript with the 
		Rscript --vanilla --quiet $RDIR/genelist2profile.R \
			_tmp/sample ${GN} ${WIDTH} ${HEIGHT} ${OS} ${OG} ${OF} &> _tmp/Routput
		rm -rf _tmp
		exit
	fi
fi
