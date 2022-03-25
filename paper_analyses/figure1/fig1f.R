setwd("~/Work/labuser/sl2665/rev/fig/fig1/fig1e/")
rm(list = ls())
source("rscript/heatmap.R")
library(dplyr)
library(tidyr)
library(ggplot2)

readmatrix = read.table("mat/pro.HEK1.matrix", stringsAsFactors=F, header=F, fill=T, flush=T)
readmatrix[is.na(readmatrix)]=-1
readmatrix = as.matrix(readmatrix)
lrm = log10(readmatrix+1)
lrm[!is.finite(lrm)] = -1000000

# Function to reduce heatmap
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

lrm.red = red_hmap(lrm, 600, 1000)
scbar = matrix( rep(rev(0:100)/100, 10), ncol = 10)

plotfig = function(filename = "pdf/fig1e.pdf", width = 4, height = 3, black.background = T) {
  pdf(filename, width = width, height = height)
  par(mar=c(2.5,2.5,0.5,0.5),mgp=c(1.5,0.4,0))
  nR = nrow(readmatrix)
  plot(0,0,type='n',xlim=c(-1,155),ylim=c(0, nR),
       axes=FALSE,xaxs='i',yaxs='i',xlab=expression("Position on gene 5"*minute %->% 3*minute),ylab='Genes sorted by length')
#  rect(-1, 0, 155, nR, col = "grey82", border = NA) 
  axis(1, at=c(0, 25, 50, 75, 100, 125, 150),
       labels = c("TSS", "25 kb", "50 kb", "75 kb", "100 kb", "125 kb", "150 kb"),
       cex.axis=0.7, tck=-0.02, lwd=0, lwd.tick=1)
  axis(2, at=nR - seq(0, nR, by = 5000), labels = seq(0, nR, by = 5000),
       cex.axis=0.7, tck=-0.02, lwd=0, lwd.tick=1, par(las=0))
  abline(v = c(0, 25, 50, 75, 100, 125, 150), col = "grey80")
  abline(h = nR - seq(0, nR, by = 5000), col = "grey80")
  heatmap.1col(lrm.red, col = "p", xlim = c(-1, 155), ylim=c(0, nR), zlim=c(0, 0.8), bl = black.background, lwd = 0.5)
  heatmap.1col(scbar, col = "p", xlim = c(140, 147.5), ylim=c(nR*0.4, nR*0.95),
	       zlim=c(0, 1), bl = black.background, lwd = 0.5)
  dev.off()
}

plotfig("pdf/fig1e.pdf", bl = F)


