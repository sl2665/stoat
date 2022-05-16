#!/bin/bash

if [ $# -lt 2 ]; then
	echo -e "Usage:\t proseq_matrix.sh [options] -a <bed> -p <pl.bedgraph> -m <mn.bedgraph>"
	echo -e "Options:"
	echo -e "\t-bin\tbin size (default=1000)"
	echo -e "\t-5\t5' flanking size (default=2000)"
	echo -e "\t-3\t3' flanking size (default=10000)"
	exit
fi

BIN=1000
FL5=2000
FL3=10000

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-a)
		BED="$2"
		shift
		;;
		-p)
		PL="$2"
		shift
		;;
		-m)
		MN="$2"
		shift
		;;
		-bin)
		BIN="$2"
		shift
		;;
		-5)
		FL5="$2"
		shift
		;;
		-3)
		FL3="$2"
		shift
		;;
		--sdir) SDIR="$2"; shift ;;
		--default)
		;;
		*)
		;;
	esac
	shift
done

awk -v f5=$FL5 -v f3=$FL3 '$6=="+"{print $1"\t"$2-f5"\t"$3+f3"\t"$4"\t"$5"\t+";next}
	{print $1"\t"$2-f3"\t"$3+f5"\t"$4"\t"$5"\t-"}' $BED | \
	awk '$2>0' | sort -k1,1 -k2,2n > .bed.tmp

$SDIR/bgmatrix -bin $BIN -a .bed.tmp -b $PL > .pl.mat.tmp
$SDIR/bgmatrix -bin $BIN -a .bed.tmp -b $MN > .mn.mat.tmp

awk 'FNR==1{++fid}
	fid==1{strand[FNR]=$6;next}
	fid==2{if(strand[FNR]=="+") {ss[FNR]=$0}
		else {as[FNR]=$NF;for(i=NF-1;i>=1;--i) as[FNR]=as[FNR]"\t"$i}
		gsub(/\t/,",",ss[FNR]); gsub(/\t/,",",as[FNR]); next}
	fid==3{if(strand[FNR]=="+") {as[FNR]=$0}
		else {ss[FNR]=$NF;for(i=NF-1;i>=1;--i) ss[FNR]=ss[FNR]"\t"$i}
		gsub(/\t/,",",ss[FNR]); gsub(/\t/,",",as[FNR]); next}
	{print $4"\t"ss[FNR]"\t"as[FNR]}' \
	.bed.tmp .pl.mat.tmp .mn.mat.tmp .bed.tmp

rm .*.tmp
