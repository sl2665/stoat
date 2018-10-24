awk '{id=$1;getline;tag=substr($1,1,8);seq=substr($1,9);getline;getline;phred=substr($1,9);if(length(seq)>=16) printf id":"tag"\n"seq"\n+\n"phred"\n"}' $1 > "${1%.*}".UMItag.fastq
