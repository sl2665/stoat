#!/bin/bash
# Generates a table of PRO-seq expression levels

if [ $# -lt 2 ]; then
	echo -e "Usage:\t proseq-getexpr [options] -p <PRO-seq filename base> -g <gene annotation bed12>"
	echo -e "Options:"
	echo -e "\t-w\tPromoter proximal range (default = 500 bp)"
	echo -e "\t--gbs\tGene body start (default = 1000 bp)"
	echo -e "\t--gbe\tGene body end (default = 5000 bp)"
	exit
fi

if [ ! -d _proseq_tmp ]; then
	mkdir _proseq_tmp
fi

RNG=500
SDIR="~/Work/tools/stoat/bin"
GBS=1000
GBE=5000
while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-p) PRO="$2"; shift ;;
		-g) BED="$2"; shift ;;
		-w) RNG="$2"; shift ;;
		--gbs) GBS="$2"; shift ;;
		--gbe) GBE="$2"; shift ;;
		--sdir) SDIR="$2"; shift ;;
		--default) ;;
		*) ;;
	esac
	shift
done

sort -k1,1 -k2,2n -k3,3n ${BED} | cut -f1-12 > _proseq_tmp/bed

# Gene body counts
bedtools map -a _proseq_tmp/bed -b ${PRO}.pl.bedgraph -o count -c 4 | cut -f4,6,13 > _proseq_tmp/_pl.gb.tmp
bedtools map -a _proseq_tmp/bed -b ${PRO}.mn.bedgraph -o count -c 4 | cut -f4,13 > _proseq_tmp/_mn.gb.tmp

# Effective gene body counts
awk -v gbs=${GBS} -v gbe=${GBE} '$6=="+"{ \
	print $1"\t"$2+gbs"\t"$2+gbe"\t"$4"\t"$5"\t"$6;next}
	{print $1"\t"$3-gbe"\t"$3-gbs"\t"$4"\t"$5"\t"$6}' _proseq_tmp/bed \
	| awk '$2>0' | sort -k1,1 -k2,2n > _proseq_tmp/gbbed

bedtools map -a _proseq_tmp/gbbed -b ${PRO}.pl.bedgraph -o count -c 4 | cut -f4,6,7 > _proseq_tmp/_pl.egb.tmp
bedtools map -a _proseq_tmp/gbbed -b ${PRO}.mn.bedgraph -o count -c 4 | cut -f4,7 > _proseq_tmp/_mn.egb.tmp

# Exon counts
bedtools map -a _proseq_tmp/bed -b ${PRO}.pl.bedgraph -split -o count -c 4 | cut -f4,6,13 > _proseq_tmp/_pl.ex.tmp
bedtools map -a _proseq_tmp/bed -b ${PRO}.mn.bedgraph -split -o count -c 4 | cut -f4,13 > _proseq_tmp/_mn.ex.tmp

# Total read counts & coverage
awk '$4>0{s+=$4;++c;next}$4<0{s-=$4;++c}END{print s"\t"c}' ${PRO}.pl.bedgraph ${PRO}.mn.bedgraph > \
	_proseq_tmp/_total_readcount.tmp

# Gene length and exon length information
awk '{s=0; split($11,a,","); for(i=1;i<=$10;++i) s+=a[i]; print $4"\t"$3-$2"\t"s}' _proseq_tmp/bed > \
	_proseq_tmp/_lengths.tmp

# Promoter proximal counts
awk -v range=${RNG} '$6=="+"{print $1"\t"$2"\t"$2+range"\t"$4"\t"$5"\t"$6;next}
	$3>range{print $1"\t"$3-range"\t"$3"\t"$4"\t"$5"\t"$6}' _proseq_tmp/bed | sort -k1,1 -k2,2n > _proseq_tmp/_pp.bed.tmp
bedtools map -a _proseq_tmp/_pp.bed.tmp -b ${PRO}.pl.bedgraph -o sum -c 4 | cut -f4,6,7 > _proseq_tmp/_pl.pp.tmp
bedtools map -a _proseq_tmp/_pp.bed.tmp -b ${PRO}.mn.bedgraph -o sum -c 4 | cut -f4,7 > _proseq_tmp/_mn.pp.tmp

# Merge two strands
paste _proseq_tmp/_pl.gb.tmp _proseq_tmp/_mn.gb.tmp | awk '$2=="+"{print $1"\t"0+$3;next}{print $1"\t"0+$5}' > \
	_proseq_tmp/_gb.tmp
paste _proseq_tmp/_pl.egb.tmp _proseq_tmp/_mn.egb.tmp | awk '$2=="+"{print $1"\t"0+$3;next}{print $1"\t"0+$5}' > \
	_proseq_tmp/_egb.tmp
paste _proseq_tmp/_pl.pp.tmp _proseq_tmp/_mn.pp.tmp | awk '$2=="+"{print $1"\t"0+$3;next}{print $1"\t"0-$5}' > \
	_proseq_tmp/_pp.tmp
paste _proseq_tmp/_pl.ex.tmp _proseq_tmp/_mn.ex.tmp | awk '$2=="+"{print $1"\t"0+$3;next}{print $1"\t"0+$5}' > \
	_proseq_tmp/_ex.tmp

# Paste all information and process into RPKM data
awk -v range=${RNG} -v gbs=${GBS} -v gbe=${GBE} 'BEGIN{print "id\tpp\tgb\tex\tRPKMpp\teRPKMgb\teRPKMex\tgbs\teRPKMgbs"}
	FNR==1{++fid}
	fid==1{tct=$1;tcv=$2;RPKMfactor=1000000/tct;eRPKMfactor=1000000*tct/tcv/tcv;next}
	fid==2{gl[$1]=$2;exl[$1]=$3;next}
	fid==3{pp[$1]=$2;next}
	fid==4{gb[$1]=$2;next}
	fid==5{ex[$1]=$2;next}
	fid==6{egb[$1]=$2;next}
	{ 	printf $4"\t"0+pp[$4]"\t"0+gb[$4]"\t"0+ex[$4]"\t";
		printf (0+pp[$4])*RPKMfactor*1000/range"\t"(0+gb[$4])*eRPKMfactor*1000/(0.0001+gl[$4])"\t";
		printf (0+ex[$4])*eRPKMfactor*1000/(0.0001+exl[$4])"\t";
		print 0+egb[$4]"\t"(0+egb[$4])*eRPKMfactor*1000/(gbe-gbs)}' \
		_proseq_tmp/_total_readcount.tmp _proseq_tmp/_lengths.tmp _proseq_tmp/_pp.tmp _proseq_tmp/_gb.tmp _proseq_tmp/_ex.tmp _proseq_tmp/_egb.tmp ${BED} 

rm -rf _proseq_tmp
