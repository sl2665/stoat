setwd("~/Work/labuser/sl2665/rev/fig/fig1/fig1f/")
rm(list = ls())
source("rscript/heatmap.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

readmatrix = read.table("mat/ted.HEK1.matrix", stringsAsFactors=F, header=F, fill=T, flush=T)
readmatrix[is.na(readmatrix)]=-1
readmatrix = as.matrix(readmatrix)
lrm2 = log10(readmatrix+1)
lrm2[!is.finite(lrm2)] = -10000
lrm3 = lrm2[,471:495]

sl.col = function(n = 7, sat = 1, lum = 1) {
	col = rev(colorRampPalette(brewer.pal(11, "Spectral"))(n)) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

red_hmap = function(mat, nrow = 400, ncol = 400) {
	nc = ncol(mat)
	nr = nrow(mat)
	
	red_col_mat = sapply(1:ncol, function(i) {
		start = floor(nc / ncol * (i-1)) + 1
		end = floor(nc / ncol * i)
		count = end - start + 1
		v = rep(0, nc)
		v[start:end] = 1/count
		return(v)})
	red_row_mat = t(sapply(1:nrow, function(i) {
		start = floor(nr / nrow * (i-1)) + 1
		end = floor(nr / nrow * i)
		count = end - start + 1
		v = rep(0, nr)
		v[start:end] = 1/count
		return(v)}))
	return(red_row_mat %*% mat %*% red_col_mat)
}

lrm2.red = red_hmap(lrm2, 600, 400)
lrm3.red = red_hmap(lrm3, 600, 25)
# TED-seq heatmap

plotfig2 = function(filename = "pdf/fig1f.pdf", width = 2.5, height = 3, black.background = T) {
  pdf(filename, width = width, height = height)
  par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.5,0.4,0))
  nR = nrow(readmatrix)
  plot(0,0,type='n',xlim=c(0,8),ylim=c(0, nR),
       axes=FALSE,xaxs='i',yaxs='i',xlab=expression("Position on transcript 3"*minute %->% 5*minute),ylab='Transcripts sorted by length')
#  rect(-1, 0, 155, nR, col = "grey82", border = NA) 
  axis(1, at=c(0:5, 5.5, 6.5, 7.5), 
       labels = c("", "-4 kb", "", "-2 kb", "", expression(3*minute*CPS), "0", "100", "200"),
       cex.axis=0.7, tck=-0.02, lwd=0, lwd.tick=1)
  axis(2, at=nR - seq(0, nR, by = 2500), labels = seq(0, nR, by = 2500),
       cex.axis=0.7, tck=-0.02, lwd=0, lwd.tick=1, par(las=0))
  y = nR - seq(0, nR, by = 2500)
  abline(v = 0:5, col = "grey80")
  segments(x0 = rep(0, length(y)), y0 = y, x1 = rep(5, length(y)), col = "grey80")
  heatmap.1col(lrm2.red, col = 'p', xlim = c(0, 5), ylim=c(0, nR), zlim=c(0, 1.5), black = black.background, lwd = 0.5)
  heatmap.1col(lrm3.red, col = 'p', xlim = c(5.5, 8), ylim=c(0, nR), zlim=c(0, 1.5), black = black.background, lwd = 0.5)
  dev.off()
}

plotfig2("pdf/fig1f.pdf", bl = F)
