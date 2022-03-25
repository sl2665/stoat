setwd("~/Work/labuser/sl2665/rev/fig/fig1/fig1c/")
rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

sl.col = function(n = 7, sat = 1, lum = 1) {
	col = colorRampPalette(brewer.pal(11, "Spectral"))(n) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

sl.fill = sl.col(20)[c(18,5,8)]
sl.color = sl.col(20,1,0.7)[c(18,5,8)]

TEDseq.dirname = c("TEDseq.HEK.rep1",
		   "TEDseq.HeLa",
		   "TEDseq.THP1")

TEDseq.samplename = NULL
for(dir in TEDseq.dirname) {
  sampleInfoTable = read.table(paste("data/", dir, "/sample_info.txt", sep=""),
                               colClasses = c("character", "character"), skip = 0)
  sampleName = sampleInfoTable[1,2]
  sampleType = sampleInfoTable[2,2]
  TEDseq.samplename = c(TEDseq.samplename, sampleName)
}

TEDseq.filename = paste("data/", TEDseq.dirname, "/table/expression.txt", sep = "")

TEDseq.table = lapply(TEDseq.filename, read.table, header = T)
names(TEDseq.table) = TEDseq.samplename

TEDseq.readcounts =
  data.frame(id = TEDseq.table[[1]][,1], sapply(TEDseq.table, function(x) x[,2])) %>%
  filter_at(vars(-id), all_vars(. > 0))

TEDseq.readcounts.long = 
  TEDseq.readcounts %>%
  gather(sample, readcount, -id) %>%
  mutate(sample = ifelse(sample == "HEK.rep1",
			"HEK293",
			sample))

pdf("pdf/fig1c.pdf", width = 2, height = 3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.4, 0))
g = ggplot(TEDseq.readcounts.long) +
  geom_histogram(aes(x = log2(readcount), fill = sample), breaks = -10:60/2) +
  stat_bin(aes(x = log2(readcount) - 0.25, col = sample), breaks = -10:60/2- 0.25, geom = "step", size = 0.4) +
  facet_grid(sample~., scales = "free_y") +
  xlim(0,15) +
  ylab("gene count")+
  xlab(expression(log[2]~read~count)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = sl.fill) +
  scale_color_manual(values = sl.color)
print(g)
dev.off()

