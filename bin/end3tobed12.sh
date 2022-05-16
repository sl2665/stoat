# 3CPS positions to bed12 format
CP3=$1
REF=$2
mkdir _t

# Last exons of ref
awk 'substr($1, 1, 1)!="#"{OFS="\t";
	split($11,blockSizes,",");
	split($12,blockStarts,",");
	blockCount=$10;
	if($6=="+")
    print $1, $2 + blockStarts[blockCount], $2 + blockStarts[blockCount] + blockSizes[blockCount], $4, $5, $6;
	else
		print $1, $2 + blockStarts[1], $2 + blockStarts[1] + blockSizes[1], $4, $5, $6;
}' $REF \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -u \
| sort -k1,1 -k2,2n -k3,3n -k6,6 > _t/lastexon.bed

# Find ref last exons for 3CPS positions
bedtools intersect -a $1 -b _t/lastexon.bed \
# Then find ref gene and truncate the last exon
| awk 'NR==FNR{cps[$10]=$3; next}
	cps[$4]{OFS = "\t";
		if($6=="+") {
			truncSize = $3 - cps[$4];
			split($11, blockSizes, ",");
			blockCount = $10;
			blockSizes[blockCount] = blockSizes[blockCount] - truncSize;
			printf $1, $2, cps[$4], $4, $5, $6, $7, $8, $9, blockCount "\t";
			for(i=1;i<=blockCount;i++) printf blockSizes[i] ",";
			print "\t"$12;
		} else {
		  truncSize = cps[$4] - $2;
			split($11, blockSizes, ",");
			split($12, blockStarts, ",");
			blockCount = $10;
      printf $1, cps[$4], $3, $4, $5, $6, $7, $8, $9, blockCount "\t";
      blockSizes[1] = blockSizes[1] - truncSize;
      for(i=2;i<=blockCount;i++) blockStarts[i] = blockStarts[i] - truncSize;
      for(i=1;i<=blockCount;i++) printf blockSizes[i] ",";
      printf "\t";
      for(i=1;i<=blockCount;i++) printf blockStarts[i] ",";
      printf "\n";
      }
    next; }
    { print $0; }
  ' - $2
  
# rm -rf _t
