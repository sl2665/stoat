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

lencor = bind_rows(TED.long %>%
		       filter(sample == "HEK" & type == "pal") %>%
		       select(-sample),
	       		LEN.long) %>%
	spread(type, val) %>%
	filter_at(vars(-id), all_vars(is.finite(.))) %>%
	mutate(x = log10(len), y = pal) %>%
	filter(y > 0 & y < 250)


pdf("pdf/fig3e.pdf", width = 3, height = 3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g = plot.cor(lencor,
	label_x = expression(log[10]~3*minute~UTR~length),
	label_y = expression(poly(A)~length),
	diag = T,
	xlim = c(1, 4.5),
	ylim = c(0, 200),
	red = 1000,
	col = sl.col(11),
	lin = TRUE,
	cor = T,
	fit="l")
print(g)
dev.off()

