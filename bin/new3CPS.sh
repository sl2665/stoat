#!/bin/bash
# Redefine new 3'CPS from TED-seq junction bam file

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat redef3"
	echo -e "version: 0.1.191018"
	echo -e ""
	echo -e "usage:   stoat redef3 [options] -d <TEDseq dir>"
	echo -e ""
	echo -e "options:"
	echo -e "         -b   3'CPS junction bam (default = stoat/refdata/HEK)"
	echo -e "         -S   antisense strand bam (default = sense strand)"
	echo -e "         -fa  reference genome fasta (default = none)"
	echo -e ""
	exit
fi

BAM=false
STR="ss"
FA=false

if [ ! -d _r3tmp ]; then
	mkdir _r3tmp
fi

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-d) TED="$2"; shift ;;
		-b) BAM="$2"; shift ;;
		-S) STR="as"; ;;
		-fa) FA="$2"; shift ;;
		--sdir) SDIR="$2"; shift;;
		--default) ;;
		*) ;;
	esac
	shift
done

if [ "$BAM" == false ]; then
	CPS=${SDIR%*/bin}/refdata/HEK.CPS
else
	if [ "$STR" == "ss" ]; then
		bedtools genomecov -ibam ${BAM} -bg -strand + -3 | \
			sort -k1,1 -k2,2n > _r3tmp/a.pl.bedgraph
		bedtools genomecov -ibam ${BAM} -bg -strand - -3 | \
			sort -k1,1 -k2,2n > _r3tmp/a.mn.bedgraph
	else
		bedtools genomecov -ibam ${BAM} -bg -strand - -5 | \
			sort -k1,1 -k2,2n > _r3tmp/a.pl.bedgraph
		bedtools genomecov -ibam ${BAM} -bg -strand + -5 | \
			sort -k1,1 -k2,2n > _r3tmp/a.mn.bedgraph
	fi
	CPS="_r3tmp/a"
fi

# Use tedseq-find3cps to find approximate 3'CPS sites from TED-seq read clusters
$SDIR/tedseq-find3cps -t $TED/alignment/a --sdir $SDIR > _r3tmp/cpsregion.bed
cp _r3tmp/cpsregion.bed $TED/annotation/3cps.bed

# Split 3'CPS sites between plus and minus strands, and take 200 bp region from the 3'ends
awk '$6=="+"{print $1"\t"$3-100"\t"$3+100 > "_r3tmp/pl.bed";next}
	$2>100{print $1"\t"$2-100"\t"$2+100 > "_r3tmp/mn.bed"}' _r3tmp/cpsregion.bed

# Extract TED-seq read peaks that mapped to the 3'CPS clusters from the 100bp short insert TED-seq library
# Get read counts for all positions within blocks
bedtools map -a _r3tmp/pl.bed -b ${CPS}.pl.bedgraph -o collapse,collapse -c 2,4 \
	> _r3tmp/pl.col
bedtools map -a _r3tmp/mn.bed -b ${CPS}.mn.bedgraph -o collapse,collapse -c 2,4 \
	> _r3tmp/mn.col

# Find the position with the maximum read count within blocks
awk '{n = split($4, pos, ","); split($5,readcount,","); \
	maxpos = 1; for(i = 1; i<=n; ++i) if(readcount[ i ]>readcount[ maxpos ]) maxpos=i; \
	print $1"\t"pos[ maxpos ]-300"\t"pos[ maxpos ]"\t"$1":"pos[ maxpos ]"\t"readcount[ maxpos ]"\t+"}' \
	_r3tmp/pl.col > _r3tmp/pl.max

awk '{n = split($4, pos, ","); split($5,readcount,","); \
	maxpos = 1; for(i = 1; i<=n; ++i) if(readcount[ i ]>readcount[ maxpos ]) maxpos=i; \
	print $1"\t"pos[ maxpos ]"\t"pos[ maxpos ]+300"\t"$1":"pos[ maxpos ]"\t"readcount[ maxpos ]"\t-"}' \
	_r3tmp/mn.col > _r3tmp/mn.max

# Merge the plus and minus strand max peak results into bed files
cat _r3tmp/pl.max _r3tmp/mn.max | \
	awk '$2>0&&$3!="."' | sort -k1,1 -k2,2n > _r3tmp/3cps.preref.bed

# If fasta file is designated, cross validate for the lack of genomic poly(A) sequences downstream
if [ "$FA" != false ]; then
	awk '$6=="+"{print $1"\t"$3"\t"$3+30"\t"$4"\t"$5"\t+";next}
		{print $1"\t"$2-30"\t"$2"\t"$4"\t"$5"\t-";}' _r3tmp/3cps.preref.bed \
		> _r3tmp/3cps.fa.bed	
	bedtools getfasta -s -tab -fi ${FA} -bed _r3tmp/3cps.fa.bed > _r3tmp/3cps.seq.fa
	awk 'NR==FNR{seq=toupper($2); if(sub("AAAAAAAAAA", "NNNNNNNNNN", seq)) flag[FNR]=1;next}
		!flag[FNR]{print $0}' _r3tmp/3cps.seq.fa _r3tmp/3cps.preref.bed \
		> $TED/annotation/3cps.ref.bed
fi

# Intersect with the original 3CPS positions in the median PAL table
# extract original 3CPS region 
awk 'NR==FNR{id[$1]=1;next}
	id[$4]{if($6=="+") print $1"\t"$3-100"\t"$3+100"\t"$4"\t"$5"\t+";
		else print $1"\t"$2-100"\t"$2+100"\t"$4"\t"$5"\t-";}' \
	$TED/table/medianpal.txt $TED/annotation/transcripts.bed13 | \
	sort -k1,1 -k2,2n -k3,3n > _r3tmp/3cpsreg.bed
awk '{split($4,a,":");print $1"\t"a[2]"\t"a[2]+1"\t"$4"\t"$5"\t"$6}' $TED/annotation/3cps.ref.bed | \
	sort -k1,1 -k2,2n -k3,3n > _r3tmp/3cpsref.bed 
bedtools intersect -a _r3tmp/3cpsreg.bed -b _r3tmp/3cpsref.bed -s -wao | sort -k4,4 -k11,11nr | \
	sort -k4,4 -u | awk '$7!="."' > _r3tmp/is.bed
cat _r3tmp/is.bed | \
	awk '$12=="+"{print $4"\t"$8-($2+100)"\t"$11;next}
		{print $4"\t"$2+100-$8"\t"$11}' > _r3tmp/redef.bed
cat _r3tmp/redef.bed | \
	awk 'NR==FNR{ref[$1]=$0;next}
		ref[$1]{print ref[$1];next}
		{print $1"\t0\t0"}' - $TED/table/medianpal.txt > \
	$TED/table/3cps_redef.txt
rm -rf _r3tmp
