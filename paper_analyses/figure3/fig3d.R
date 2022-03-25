# Pre-run fig3.R before running this script 

source("rscript/scatterPlot.R")

library(RColorBrewer)
sl.col = function(n = 7, sat = 1, lum = 0.85) {
	col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(n) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

aucor = bind_rows(TED.long %>%
		       filter(sample == "HEK" & type == "pal") %>%
		       select(-sample),
	       		AU.long) %>%
	spread(type, val) %>%
	filter_at(vars(-id), all_vars(is.finite(.))) %>%
	mutate(x = au*100, y = pal) %>%
	filter(y > 0 & y < 250)


pdf("pdf/fig3d.pdf", width = 3, height = 3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g = plot.cor(aucor,
	label_x = expression(AU~content~"(%)"),
	label_y = expression(poly(A)~length),
	diag = T,
	xlim = c(20, 90),
	ylim = c(0, 200),
	red = 1000,
	col = sl.col(11),
	lin = TRUE,
	cor = T,
	fit = "l")
print(g)
dev.off()

