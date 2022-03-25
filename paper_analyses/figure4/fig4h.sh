# Merge all gene prediction tracks
cd geneprediction
cat augustusGene.bed  geneid.bed  genscan.bed  sgpGene.bed  sibGene.bed > all.bed
cd ..


# Get all predicted 5' splice sites
awk '{OFS="\t";
	split($11,blockSizes,",");
	split($12,blockStarts,",");
	blockCount=$10;
	if($6=="+") for(i=1;i<=blockCount-1;i++)
		print $1, $2 + blockStarts[i] + blockSizes[i], $2 + blockStarts[i] + blockSizes[i] + 1, ".", $5, $6;
	else for(i=2;i<=blockCount;i++)
		print $1, $2 + blockStarts[i], $2 + blockStarts[i] + 1, ".", $5, $6;
}' geneprediction/all.bed | \
sort -k1,1 -k2,2n -k3,3n -u > geneprediction/5ss.bed

# 5' splice sites in the last exons
awk '{OFS="\t";
	split($11,blockSizes,",");
	split($12,blockStarts,",");
	blockCount=$10;
	if($6=="+")
		print $1, $2 + blockStarts[blockCount], $2 + blockStarts[blockCount] + blockSizes[blockCount], $4, $5, $6;
	else
		print $1, $2 + blockStarts[1], $2 + blockStarts[1] + blockSizes[1], $4, $5, $6;
}' bed/gencode.v26.geneName.bed13 | \
sort -k1,1 -k2,2n -k3,3n -k6,6 -u | \
sort -k1,1 -k2,2n -k3,3n -k6,6 > geneprediction/lastexon.bed

bedtools intersect -s -wa -wb -a geneprediction/5ss.bed -b geneprediction/lastexon.bed \
| awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$6}' \
> geneprediction/5ss.lastex.bed

# All introns
awk 'NR>1{OFS="\t";
	split($11,blockSizes,",");
	split($12,blockStarts,",");
	blockCount=$10;
	for(i=1;i<blockCount;i++)
		print $1, $2 + blockStarts[i] + blockSizes[i], $2 + blockStarts[i+1], $4, $5, $6;
}' bed/gencode.v26.geneName.bed13  \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -u \
> geneprediction/allintron.bed

# Find immediate upstream intronic 5' splice site in the same intron for intronic TSSs
# Label which intron 5'ss overlap
bedtools intersect -wa -wb -s -a geneprediction/5ss.bed -b geneprediction/allintron.bed \
| awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$8"-"$9"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 \
> geneprediction/5ss.intron.bed
# Label which intron intronic TSSs overlap
bedtools intersect -wa -wb -s -a bed/dTSSaCPS.bed6 -b geneprediction/allintron.bed \
| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"-"$9"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 \
> geneprediction/pCPS.intron.bed
# Find immediate upstream 5'ss in the same intron
| awk '$5==$10{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -u \
> geneprediction/pCPS.5ss.bed

# Find immediate upstream 3'UTR 5' splice site in the same intron for intronic TSSs
# Label which intron 5'ss overlap
bedtools intersect -wa -wb -s -a geneprediction/5ss.bed -b geneprediction/lastexon.bed \
| awk '{print $1"\t"$2"\t"$3"\t"$8"-"$9"\t0\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 \
> geneprediction/5ss.lastexon.bed
# Label which last exon TSSs overlap
bedtools intersect -wa -wb -s -a bed/dTSSamCPS.bed6 -b geneprediction/lastexon.bed \
| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"-"$9"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 \
> geneprediction/mCPS.lastexon.bed
# Find immediate upstream 5'ss in the same last exon
bedtools closest -id -s -D a -a geneprediction/mCPS.lastexon.bed -b geneprediction/5ss.lastexon.bed \
| awk '$5==$10{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -u \
> geneprediction/mCPS.5ss.bed

# Find immediate upstream 3'UTR 5' splice sites for all CPSs
# Label which last exon TSSs overlap
awk '$6=="+"{print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t+";next}
	{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t-"}' bed/gencode.v26.geneName.bed13 \
| bedtools intersect -wa -wb -s -a - -b geneprediction/lastexon.bed \
| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"-"$9"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 \
> geneprediction/CPS.lastexon.bed
# Find immediate upstream 5'ss in the same last exon
bedtools closest -id -s -D a -a geneprediction/CPS.lastexon.bed -b geneprediction/5ss.lastexon.bed \
| awk '$5==$10{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 -u \
> geneprediction/CPS.5ss.bed
# Find immediate upstream 5'ss in the previous intron
awk '$6=="+"{print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6; next}
  {print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6}' geneprediction/allintron.bed \
| sort -k1,1 -k2,2n -k3,3n \
| bedtools closest -id -s -D a -a geneprediction/CPS.lastexon.bed -b - \
| awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 -u \
> geneprediction/CPS.up5ss.bed
bedtools closest -id -s -D a -a geneprediction/CPS.lastexon.bed -b geneprediction/5ss.lastexon.bed \
| awk '$5!=$10{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$6}' \
| sort -k1,1 -k2,2n -k3,3n -k6,6 -u \
> geneprediction/CPS.no5ss.bed
