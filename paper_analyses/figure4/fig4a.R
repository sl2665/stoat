setwd("~/Work/labuser/sl2665/rev/fig/fig4")
rm(list = ls())

source("rscript/browser.R")
library(RColorBrewer)
library(dplyr)
library(tidyr)

# Track filename lists
bglist = paste0("bedgraph/",
		c("pro.hek.pl",
		  "pro.hek.mn",
		  "ted.hek.pl",
		  "ted.hek.mn",
		  "rna.hek.pl",
		  "rna.hek.mn"),
		".bg")
color.pm = colorRampPalette(brewer.pal(11, "Spectral"))(20)
color = c(color.pm[c(3,19)],
	  color.pm[c(5,18)],
	  color.pm[c(7,16)])
# Plot reference annotation and PRO-seq
browser.setup(genelist="bed/gencode.v26.uniq.bed",
	      bedgraph = bglist[1:2],
	      col = color[1:2], heights = c(1,1), hlines = 2,
	      label.y = 1, description = "PRO-seq", gene.line = 1)
browser.setgene("B3GALT6")
browser.setpos("chr1", 1230000, 1237000)
browser.read()
ymax = c(400,400)
ymax = browser.print(filename="pdf/fig4a1.pdf", nbin=150, width=6, height = 1.8, ymax=ymax)
print(ymax)
# Plot reference annotation and TED-seq
browser.new()
browser.setup(genelist="bed/3cps.ref.bed12",
	      bedgraph = bglist[3:4],
	      col = color[3:4], heights = c(1,1), hlines = 2,
	      label.y = 1, description = "TED-seq", gene.line = 1,
	      gene.col = color[3:4])
browser.setpos("chr1", 1230000, 1237000)
browser.read()
ymax = c(200,200)
ymax = browser.print(filename="pdf/fig4a2.pdf", nbin=150, width=6, height = 1.8, ymax=ymax)
print(ymax)
browser.new()
browser.setup(genelist="bed/RNAseq.cufflinks.bed",
	      bedgraph = bglist[5:6],
	      col = color[5:6], heights = c(2,2), hlines = 4,
	      label.y = 2, description = "RNA-seq", gene.line = 4,
	      gene.col = color[5:6])
browser.setpos("chr1", 1230000, 1237000)
browser.read()
ymax = c(70,70)
ymax = browser.print(filename="pdf/fig4a3.pdf", nbin=150, width=6, height = 2, ymax=ymax)
print(ymax)
