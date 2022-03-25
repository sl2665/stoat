setwd("~/Work/labuser/sl2665/rev/fig/fig1/fig1b/")
source("rscripts/browser.R")
library(RColorBrewer)
library(dplyr)
library(tidyr)

# Set up gene list
browser.setup(genelist="bed/gencode.v26.uniq.bed") # Set up gene list
# Track heights
heights = rep(0.5, 4)
# Horizontal line positions
hlines = c(1, 2)
# Track label descriptions
description = rev(c("PRO-seq", "TED-seq"))
# Track label positions
label.y = c(0.5, 1.5)
# Track filename lists
bglist.all = paste0("bedgraph/",
		c("pro.hek.pl",
		  "pro.hek.mn",
		  "ted.hek.pl",
		  "ted.hek.mn"),
		".bg")
# Select PRO-seq data
bglist = bglist.all
color.pm = colorRampPalette(brewer.pal(11, "Spectral"))(20)
color = c(color.pm[c(3,19)],
	  color.pm[c(5,18)])

browser.setup(bedgraph = bglist,col = color, heights = heights, hlines = hlines,
	      label.y = label.y, description = description)
#browser.setgene("BTF3")

#BTF3 whols span
browser.setpos("chr5", 73496000, 73508000)
browser.read()
ymax = c(300, 300, 12500, 40)
ymax = browser.print(filename="pdf/fig1b1.pdf", nbin=150, width=6, height = 2.5, ymax=ymax)

#BTF uaRNA poly(A)
browser.setup(genelist="bed/gencode.v26.uniq.bed") # Set up gene list
# Track heights
heights = c(0.25, 0.75)
# Horizontal line positions
hlines = c(1)
# Track label descriptions
description = rev(c("TED-seq"))
# Track label positions
label.y = c(0.5)
# Select data
bglist = bglist.all[3:4]
color.pm = colorRampPalette(brewer.pal(11, "Spectral"))(20)
color = c(color.pm[c(5,18)])
browser.setup(bedgraph = bglist,col = color, heights = heights, hlines = hlines,
	      label.y = label.y, description = description)
browser.setpos("chr5", 73497500, 73498000)
browser.read()
ymax = c(3000,10)
ymax = browser.print(filename="pdf/fig1b2.pdf", nbin=100, width=3, height = 1.5, ymax=ymax)

# BTF transciript poly(A)
heights = c(0.75, 0.25)
browser.setup(bedgraph = bglist,col = color, heights = heights, hlines = hlines,
	      label.y = label.y, description = description)
browser.setpos("chr5", 73505150, 73505650)
#browser.setpos("chr5", 73498400, 73498900)
browser.read()
ymax = c(3000, 10)
ymax = browser.print(filename="pdf/fig1b3.pdf", nbin=100, width=3, height = 1.5, ymax=ymax)
