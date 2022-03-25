# Use gene id list for all genes
stoat pap -d data/TEDseq.HEK.rep1 data/TEDseq.HeLa data/TEDseq.THP1 \
	-g fig3a/genelist/genes.all.txt --de HEK293 HeLa THP1 -o pdf/fig3a.pdf \
	-w 6 -h 3

