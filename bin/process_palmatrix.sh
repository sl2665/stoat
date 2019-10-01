# arguments : directory startpos maxlen
awk -v spos=$2 -v maxl=$3 -v dir=$1 '{end=spos+maxl;if(spos+maxl>NF-1) end=NF-2;
	s=0;for(i=spos;i<=end;++i) s+=$(i+2);
	ss=0; for(i=spos;i<=end;++i) {ss+=$(i+2);if(ss>=s/2) {mlen=i-spos; break;}}
	if(s>5) {print $0 > dir"/table/palmatrix.expr.txt";
		print $1"\t"mlen > dir"/table/medianpal.txt";}}' \
	${1}/table/palmatrix.txt

Rscript --vanilla --quiet -e 'arg=commandArgs(trailingOnly=T); dir=arg[1];spos=as.numeric(arg[2]);maxl=as.numeric(arg[3]);' \
	-e 'mat=read.table(paste0(dir,"/table/palmatrix.expr.txt"),header=F,fill=T,flush=T);' \
	-e 'colnames(mat)=c("id",(-spos):(500-spos));' \
	-e 'library(dplyr);library(tidyr);mat.long=mat%>%gather(pos,count,-id)%>%mutate(pos=as.numeric(pos),count=as.numeric(count));' \
	-e 'saveRDS(list(mat=mat,long=mat.long),paste0(dir,"/Rdata/PALdata.rds"));' \
	${1} ${2} ${3} &> $1/_tmp/Routput
