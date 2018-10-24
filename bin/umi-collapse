###############################################################################
# Generate UMI collapsed bam file
###############################################################################

samtools view -S $1 | awk '{print substr($1,length($1)-7,8)"\t"$0;}' > sam.tmp
samtools view -SH $1 > umi.unique.sam.tmp
sort --parallel=4 -S 50% -k4,4 -k5,5n -k1,1 -u sam.tmp | cut -f2- >> umi.unique.sam.tmp
samtools view -Sb umi.unique.sam.tmp > "${1%.*}".uniqueUMI.bam
rm *.tmp
