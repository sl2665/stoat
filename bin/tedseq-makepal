#!/bin/bash
# Parse TED-seq bam files based on soft clipping status.

if [ $# -lt 2 ]; then
	echo -e "Usage:\t tedseq-makepal [options] -a <bam> -b <bed>"
	echo -e "Options:"
	echo -e "\t-bin\tbin size (default=1)"
	echo -e "\t-win\twindow size (default=500)"
	exit
fi

BIN=1
WIN=500

if [ ! -d _tedseq_tmp ]; then
	mkdir _tedseq_tmp
fi

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-a)
		BAM="$2"
		shift
		;;
		-b)
		BED="$2"
		shift
		;;
		-bin)
		BIN="$2"
		shift
		;;
		-win)
		WIN="$2"
		shift
		;;
		--default)
		;;
		*)
		;;
	esac
	shift
done

samtools view -H $BAM > _tedseq_tmp/_sam.header.tmp
samtools view $BAM | awk '!($2 AND 16) {
		if($6~/^[0-9]+S/) {
			split($6,a,"S");
			outfile = "_tedseq_tmp/_sam." a[1] "S.tmp";
			print $0 > outfile;}
		else {print $0 > "_tedseq_tmp/_sam.0S.tmp";}
		next;
	}
	{
		if($6~/[0-9]+S$/) {n=split($6,a,"M");
			outfile = "_tedseq_tmp/_sam." a[n] ".tmp";
			print $0 > outfile;}
		else {print $0 > "_tedseq_tmp/_sam.0S.tmp";}
	}'
for f in _tedseq_tmp/_sam.*S.tmp; do
	outsuf=${f#_tedseq_tmp/_sam.*}
	cat _tedseq_tmp/_sam.header.tmp $f > _tedseq_tmp/_sam2.$outsuf
	samtools view -Sb _tedseq_tmp/_sam2.$outsuf > _tedseq_tmp/_bam.$outsuf
done

# Generate bedgraph files
for f in _tedseq_tmp/_bam.*S.tmp; do
	outsuf=${f#_tedseq_tmp/_bam.*}
	bedtools genomecov -ibam _tedseq_tmp/_bam.$outsuf -bg -5 -strand + > _tedseq_tmp/_bg.pl.$outsuf
	bedtools genomecov -ibam _tedseq_tmp/_bam.$outsuf -bg -5 -strand - > _tedseq_tmp/_bg.mn.$outsuf
done

# Generate last window size bed file with 1 base extension of the blocks
bed12toend3 $BED $WIN | sort -k1,1 -k2,2n | \
awk '{split($11,exonSizes,",");split($12,exonStarts,",");exonCount=$10;A="";B="";
	if($6=="+") {
		for(i=exonCount;i>1;--i) {
			--exonStarts[i];++exonSizes[i];--exonSizes[i-1];}
		for(i=1;i<=exonCount;++i) {
			A=A exonSizes[i]",";B=B exonStarts[i]",";}
		OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,A,B > "_tedseq_tmp/_bed.500.pl.tmp";
	} else {
		for(i=1;i<exonCount;++i) {
			++exonSizes[i];--exonSizes[i+1];++exonStarts[i+1];}
		for(i=1;i<=exonCount;++i) {
			A=A exonSizes[i]",";B=B exonStarts[i]",";}
		OFS="\t";print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,A,B > "_tedseq_tmp/_bed.500.mn.tmp";
	}
}'

# Generate matrices
for f in _tedseq_tmp/_bg.pl.*S.tmp; do
	outsuf=${f#_tedseq_tmp/_bg.pl.*}
	if [ -s _tedseq_tmp/_bg.pl.$outsuf ] && [ -s _tedseq_tmp/_bed.500.pl.tmp ]; then
		bgmatrix -a _tedseq_tmp/_bed.500.pl.tmp -b _tedseq_tmp/_bg.pl.$outsuf -split -bin 1 > _tedseq_tmp/_mat.pl.$outsuf 
	fi
	if [ -s _tedseq_tmp/_bg.mn.$outsuf ] && [ -s _tedseq_tmp/_bed.500.mn.tmp ]; then
		bgmatrix -a _tedseq_tmp/_bed.500.mn.tmp -b _tedseq_tmp/_bg.mn.$outsuf -split -bin 1 > _tedseq_tmp/_mat.mn.$outsuf 
	fi
done

cp _tedseq_tmp/_mat.pl.0S.tmp _tedseq_tmp/_mat.pl.tmp
cp _tedseq_tmp/_mat.mn.0S.tmp _tedseq_tmp/_mat.mn.tmp
mv _tedseq_tmp/_mat.pl.0S.tmp _tedseq_tmp/_mat0.pl.tmp
mv _tedseq_tmp/_mat.mn.0S.tmp _tedseq_tmp/_mat0.mn.tmp

# Collapse matrices
for f in _tedseq_tmp/_mat.pl.*S.tmp; do
	outsuf=${f#_tedseq_tmp/_mat.pl.*}
	if [ -s _tedseq_tmp/_mat.pl.$outsuf ]; then
		paste _tedseq_tmp/_mat.pl.tmp _tedseq_tmp/_mat.pl.$outsuf > _tedseq_tmp/_mat.pl.pasted.tmp
		awk -v os=$outsuf 'BEGIN{split(os,suf,"S"); shiftBase=0+suf[1];}
			{	printf $1+$(1+NF/2+shiftBase);
				for(i=2;i<=NF/2-shiftBase;++i) printf "\t"$i+$(i+NF/2+shiftBase);
				for(i=NF/2-shiftBase+1;i<=NF/2;++i) printf "\t"$i; printf "\n";}' \
			_tedseq_tmp/_mat.pl.pasted.tmp> _tedseq_tmp/_mat.pl.tmp
	fi
	if [ -s _tedseq_tmp/_mat.mn.$outsuf ]; then
		paste _tedseq_tmp/_mat.mn.tmp _tedseq_tmp/_mat.mn.$outsuf > _tedseq_tmp/_mat.mn.pasted.tmp
		awk -v os=$outsuf 'BEGIN{split(os,suf,"S"); shiftBase=0+suf[1];}
			{	printf $1;
				for(i=2;i<=shiftBase;++i) printf "\t"$i;
				for(i=shiftBase+1;i<=NF/2;++i) printf "\t"$i+$(i+NF/2-shiftBase);
				printf "\n";}' \
			_tedseq_tmp/_mat.mn.pasted.tmp> _tedseq_tmp/_mat.mn.tmp
	fi
done

# Finalize output
cut -f4 _tedseq_tmp/_bed.500.pl.tmp | paste - _tedseq_tmp/_mat.pl.tmp > _tedseq_tmp/_mat.pl.id.tmp
cut -f4 _tedseq_tmp/_bed.500.mn.tmp | paste - <(awk '{printf $NF; for(i=NF-1;i>=1;--i) printf "\t"$i; printf "\n";}' _tedseq_tmp/_mat.mn.tmp ) > _tedseq_tmp/_mat.mn.id.tmp
awk 'NR==1{fc=NF}FNR==1{++fid}fid==1{mat[$1]=$0;next}fid==2{mat[$1]=$0;next}
	mat[$4]{print mat[$4];next}{printf $4;for(i=1;i<fc;++i) printf "\t0";printf "\n"}' _tedseq_tmp/_mat.pl.id.tmp _tedseq_tmp/_mat.mn.id.tmp $BED
rm -rf _tedseq_tmp
