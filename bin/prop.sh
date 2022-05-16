#!/bin/bash
# Generate PRO-seq profile 

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat prop - PRO-seq profile"
	echo -e "version: 0.1.191001"
	echo -e ""
	echo -e "usage:   stoat prop [options] -d <PROseq dirs> -g <gene lists txt>"
	echo -e ""
	echo -e "options:"
	echo -e "         -b   TSS bed files instead of gene list (default = none)"
	echo -e "         -o   output plot file name (default = proout.pdf)"
	echo -e "         -w   width of the output pdf (default = 6)"
	echo -e "         -h   height of the output pdf (default = 4)"
	echo -e "         --sc scaled gene body profile (default = TSS proximal)"
	echo -e "         --de sample/gene group descriptions (default = none)"
	echo -e "         --ol overlay plots (default = separate)"
	echo -e "         --lg y axis in log scale (default = linear)"
	echo -e ""
	exit
fi

BED=false
OUT=false
GN=false
DSC=false
OF="pro.out.pdf"
WIDTH=6
HEIGHT=4
SC=false
LOG=false
OL=false

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-d) PRO=("$2"); while [[ ${3:0:1} != "-" ]]; do PRO+=("$3"); shift; done; shift ;;
		-g) GN="$2"; shift ;;
		-o) OF="$2"; shift ;;
		-w) WIDTH="$2"; shift ;;
		-h) HEIGHT="$2"; shift ;;
		-b) BED=("$2"); while [[ ${3:0:1} != "-" ]]; do BED+=("$3"); shift; done; shift ;;
		--de) DSC=("$2"); while [[ ${3:0:1} != "-" ]]; do DSC+=("$3"); shift; done; shift ;;
		--sdir) SDIR="$2"; shift ;;
		--sc) SC=true ;;
		--lg) LOG=true ;;
		--ol) OL=true ;;
		*) ;;
	esac
	shift
done

# Directories for plotting PRO-seq profiles
RDIR=${SDIR%*/bin}/rscript
if [ ! -d _tmp ]; then mkdir _tmp; 
else rm -rf _tmp; mkdir _tmp;
fi

# PROseq sample informations
for ((i=0;i<${#PRO[@]};++i)); do
	awk 'NR==1{print $2}' ${PRO[i]}/sample_info.txt >> _tmp/sample
done

# For producing PRO-seq profiles from new bed files
if [ "$BED" != false ] ; then
	if [ "$DSC" == false ] ; then DSC=("${BED[@]}"); fi
	# process multiple bed files to one adding index
	for ((i=0;i<${#BED[@]};++i)); do
		awk -v i=${DSC[i]} '{print $0"\t"i}' ${BED[i]} >> _tmp/bed
	done
	sort -k1,1 -k2,2n -k3,3n -u _tmp/bed > _tmp/tmp
	mv _tmp/tmp _tmp/bed
else
# For producing bed file from gene list txt file (with or without grouping)
	# Bed file as the reference
	cut -f1-6,13 ${PRO[0]}/annotation/transcripts.bed13 > _tmp/annot
	# Expression level from PRO-seq data
	cut -f1,6 ${PRO[0]}/table/expression.txt > _tmp/expre
	# Select isoforms with the highest expression level
	awk 'NR==FNR{name[$4]=$7;bed[$4]=$0;next}
		{print $1"\t"name[$1]"\t"$2"\t"bed[$1]}' _tmp/annot _tmp/expre | \
		sort -k2,2 -k2,2nr | sort -k2,2 -u | cut -f4- > _tmp/uannot  
	
	# Find gene ids from the reference annotation
	awk 'NR==FNR{bed4id[$4]=$0;ids4name[$7]=$4","ids4name[$7];next}
		bed4id[$1]{print bed4id[$1];next}
		ids4name[$1]{n=split(ids4name[$1],a,",");
			if(NF==1) for(i=1;i<n;++i) print bed4id[a[i]];
			else for(i=1;i<n;++i) print bed4id[a[i]]"\t"$2;
		}' _tmp/uannot ${GN} | awk '{print $0"\t"NR}' \
		| sort -k1,1 -k2,2n -k3,3n -u > _tmp/bed
fi

# From the bed file to generate proseq matrix 
if [ "$SC" == true ] ; then
	# For scaled metagene, generate the whole GB matrix in 100 bp resolution
	for ((i=0;i<${#PRO[@]};++i)); do
		$SDIR/proseq_matrix.sh -a _tmp/bed \
			-p ${PRO[i]}/alignment/a.pl.bedgraph -m ${PRO[i]}/alignment/a.mn.bedgraph \
			-bin 100 --sdir $SDIR > _tmp/mat_$i
	done
else
	# For non-scaled metagene, generate high res promoter proximal matrix
	awk '$6=="+"{if($2>1000) { printf $1"\t"$2-1000"\t"$2+1000"\t"$4"\t"$5"\t+";
			for(i=7;i<=NF;++i) printf "\t"$i; printf "\n"} next}
		$3>1000{printf $1"\t"$3-1000"\t"$3+1000"\t"$4"\t"$5"\t-";
			for(i=7;i<=NF;++i) printf "\t"$i; printf "\n"; next}' \
		_tmp/bed > _tmp/prmbed2
	# create proseq matrix
	for ((i=0;i<${#PRO[@]};++i)); do
		$SDIR/proseq_matrix.sh -a _tmp/prmbed2 -3 0 -5 0\
			-p ${PRO[i]}/alignment/a.pl.bedgraph -m ${PRO[i]}/alignment/a.mn.bedgraph \
			-bin 5 --sdir $SDIR > _tmp/mat_$i
	done
	if [ "$GN" != false ] ; then
		sort -k8,8n _tmp/prmbed2 > _tmp/_tmp
		cut -f1-7 _tmp/_tmp > _tmp/prmbed
	fi
fi
# Calculate normalization factors (total reads in million)
for ((i=0;i<${#PRO[@]};++i)); do
	awk '$4>0{s+=($4);next}{s-=$4}END{print s/1000000}' ${PRO[i]}/alignment/a.pl.bedgraph \
		${PRO[i]}/alignment/a.mn.bedgraph >> _tmp/trc
done

# From the matrix generate plot using R
if [ "$SC" == false ] ; then
	# promoter proximal profile
	Rscript --vanilla --quiet $RDIR/proseq_profile.R \
		${OF} ${WIDTH} ${HEIGHT} ${#PRO[@]} &> _tmp/Routput
else
	# scaled profile
	Rscript --vanilla --quiet $RDIR/proseq_meta.R \
		${OF} ${WIDTH} ${HEIGHT} ${LOG} ${OL} &> _tmp/Routput
fi

rm -rf _tmp
