#/bin/bash
# Align TED2 data to make bam files

if [ $# -lt 4 ]; then
	echo -e ""
	echo -e "tool:    stoat align-ted2"
	echo -e "version: 0.1.190924"
	echo -e ""
	echo -e "usage:   stoat align-ted2 [options] -1 <read 1> -2 <read 2> -r <reference genome>"
	echo -e ""
	echo -e "options:"
	echo -e "         -o   output filename base (default = ted2.out)"
	echo -e ""
	exit
fi

BASE="ted2.out"
BARCODES=(
	"TAC"
	"CGT"
	"CTG"
	"ACT"
	"GAG"
	"GTT"
	"ACA"
	"CCT"
	"CGG"
	"GAT"
)

if [ ! -d _ted2_tmp ]; then
	mkdir _ted2_tmp
fi

while [[ $# -ge 1 ]]; do
	key="$1"
	case $key in
		-1)
		READ1="$2"
		shift
		;;
		-2)
		READ2="$2"
		shift
		;;
		-r)
		REF="$2"
		shift
		;;
		-b)
		BASE="$2"
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

# UMI tag read 2 to read 1, and extract sample barcode from read 1. 
awk -v r2=$READ2 \
'{	id=$1;	getline;
	bc=substr($1,1,3);	seq=substr($1,5);	getline;	getline;
	phred=substr($1,5);
	getline < r2; getline < r2;
	umi=substr($1,1,8);	polya=substr($1,9,24);
	getline < r2; getline < r2;
	if(umi != "NNNNNNNN") printf id":"polya":"bc":"umi"\n"seq"\n+\n"phred"\n"}' $READ1 > _ted2_tmp/umitag.fastq

# Alignemnt to the reference genome
star --genomeDir ${REF} \
		--readFilesIn _ted2_tmp/umitag.fastq \
		--outFilterMultimapNmax 1 --runThreadN 8 \
		--outFileNamePrefix _ted2_tmp/ 

# Collapse same UMIs at the same mapped positions 
samtools view -S _ted2_tmp/Aligned.out.sam | awk '{n=length($1); print substr($1,n-11,12)"\t"$0;}' > _ted2_tmp/sam.tmp
samtools view -SH _ted2_tmp/Aligned.out.sam > _ted2_tmp/header.tmp
sort -k4,4 -k5,5n -k1,1 -u _ted2_tmp/sam.tmp | cut -f2- > _ted2_tmp/umi.unique.sam.tmp

# Split to different files by the barcodes
awk '{n=length($1); print $0 > "_ted2_tmp/"substr($1,n-11,3)".bctmp"}' _ted2_tmp/umi.unique.sam.tmp

for i in {1..10}; do
	BC=${BARCODES[$i-1]}
	# Convert Barcode sam to bam
	 cat _ted2_tmp/header.tmp _ted2_tmp/${BC}.bctmp > _ted2_tmp/bctmp.sam
	 samtools view -Sb _ted2_tmp/bctmp.sam > ${BASE}.${i}.bam
done

# Remove temp directory
rm -rf _ted2_tmp
