# Pre-run fig3.R before running this script 

sl.col = function(n = 7, sat = 1, lum = 0.85) {
	col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(n) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

m3s = mat.sample("HeLa")
hm.hela = cor(m3s[,-1])

pdf("pdf/fig3g.pdf", width = 3, height = 3.5)
pheatmap(hm.hela,
	 main = "HeLa",
	 labels_col = c("AU content",
	"Codon optimality",
	"Gene body PRO-seq",
	"Gene length",
	expression(3*minute*UTR~length),
	"Poly(A) length",
	"Promoter PRO-seq",
	expression(2*degree~structure),
	"TED-seq expression"),
	labels_row = c("AU", "CO", "GB", "GL","3UL","PAL","PP","2S","TED"),
	breaks = -43:43/43*0.7,
	color = sl.col(101, 1, 1)[7:93], 
	treeheight_row = 10,
	treeheight_col = 10,
	angle_col = 315)
dev.off()
