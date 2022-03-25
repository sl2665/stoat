setwd("~/Work/labuser/sl2665/rev/fig/fig4")
rm(list = ls())


# CPS count table
cpsCount = data.frame(
CPScount = factor(c(
1,
2,
3,
4,
5,
">5"), levels = c(1:5, ">5")
),
RNAseq = c(
11583,
842,
246,
82,
33,
45),
TEDseq= c(
10415,
2146,
664,
211,
101,
90))
cpsCount = cpsCount %>%
	gather(type, count, -CPScount)

sl.col = function(n = 7, sat = 1, lum = 0.85) {
	col = colorRampPalette((brewer.pal(11, "Spectral")))(n) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

#sl.fill = sl.col(20,0.9,1.1)[c(18,5,8,16,14,10,12,20)]
#sl.cols = sl.col(20,1,0.7)[c(18,5,8,16,14,10,12,20)]
sl.fill = sl.col(8,0.9,1.1)[2:7] %>% rev
sl.cols = sl.col(8,1,0.7)[2:7] %>% rev

library(ggplot2)
library(RColorBrewer)
pdf("pdf/fig4b.pdf", width=5, height=1.8)
par(mar=c(2,2,1,1),mgp=c(1.5,0.4,0))
g = ggplot(cpsCount, aes(x=type, y=count, group = CPScount, col = CPScount, fill = CPScount)) +
	geom_bar(stat = "identity", position = "fill") +
	coord_flip() +
	theme_bw() +
	scale_fill_manual(values = sl.fill,
			  guide = guide_legend(title=expression(3*minute*CPS~count),
					       title.position = "left",
					       nrow = 1)) +
	scale_color_manual(values = sl.cols,
			  guide = guide_legend(title=expression(3*minute*CPS~count),
					       title.position = "left",
					       nrow = 1)) +
	theme(legend.position="bottom", legend.key = element_rect()) +
	labs(x = "", y = "Fraction")

print(g)
dev.off()

