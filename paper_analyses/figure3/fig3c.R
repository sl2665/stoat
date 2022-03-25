# Pre-run fig3.R before running this script 
library(pheatmap)

sl.col = function(n = 7, sat = 1, lum = 0.85) {
	col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(n) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

m3c = mat.stab("HEK")
hm.s = cor(m3c[,-1])

library(pheatmap)
pdf("pdf/fig3c.pdf", width = 3.5, height = 4)
pheatmap(hm.s, labels_col = c("AU content",
	"Codon optimality",
	"Gene body PRO-seq",
	"Gene length",
	"mRNA half life",
	expression(3*minute*UTR~length),
	"Poly(A) length",
	"Promoter PRO-seq",
	"RNA-seq expression",
	expression(2*degree~structure),
	"TED-seq expression"),
	labels_row = c("AU", "CO", "GB", "GL","mHL","3UL","PAL","PP","RNA","2S","TED"),
	breaks = -43:43/43*0.7,
	color = sl.col(101, 1, 1)[7:93], 
	treeheight_row = 10,
	treeheight_col = 10,
	angle_col = 315
	)
dev.off()
