GTF=$1
ODIR=$2

# Convert gtf to bed12
gtfToGenePred $GTF $ODIR/_tmp/pred -ignoreGroupsWithoutExons
genePredToBed $ODIR/_tmp/pred $ODIR/_tmp/bed
# Match transcript ids to gene names
awk -F"\t" '{n=split($9,a," "); gid=""; tid=""; name="";
	for(i=0;i<=n/2;++i) { 
		c = a[i*2+1];
		if(c == "gene_id") gid=a[i*2+2];
		else if(c == "transcript_id") tid=a[i*2+2];
		else if(c == "gene_name") name=a[i*2+2];
	}
	if(tid!="") {gid=substr(gid,2,length(gid)-3);tid=substr(tid,2,length(tid)-3);
		name=substr(name,2,length(name)-3);geneName[tid]=name;geneID[tid]=gid;}
	}END{for(i in geneID) print i"\t"geneID[i]"\t"geneName[i];}' \
	$GTF > $ODIR/_tmp/names
# Output bed file
awk 'NR==FNR{name[$1]=$3;next}{print $0"\t"name[$4]}' $ODIR/_tmp/names \
	$ODIR/_tmp/bed > $ODIR/annotation/transcripts.bed13	 
