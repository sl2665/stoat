setwd("~/Work/labuser/sl2665/rev/fig/fig2/")
rm(list = ls())

library(dplyr)
library(tidyr)

library(ggplot2)
library(RColorBrewer)
source("rscript/scatterPlot.R")

sl.col = function(n = 7, sat = 1, lum = 0.85) {
	col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(n) %>%
	col2rgb %>%
	rgb2hsv
	col[2,] = col[2,] * sat
	col[3,] = col[3,] * lum
	col[3,col[3,]>1] = 1
	return(hsv(col[1,], col[2,], col[3,])) }

#sl.fill = sl.col(20,0.9,1.1)[c(18,5,8,16,14,10,12,20)]
#sl.cols = sl.col(20,1,0.7)[c(18,5,8,16,14,10,12,20)]


# Read and process stability data
st.4s.1 = read.csv("Stability/GSE99517_HEK.4SU.HL.DRUID.one.csv")
st.4s.2 = read.csv("Stability/GSE99517_HEK.4SU.HL.DRUID.two.csv")
st.aa.1 = read.csv("Stability/GSE99517_Aman.HEK.HL.one.csv")
st.aa.2 = read.csv("Stability/GSE99517_Aman.HEK.HL.two.csv")
st.ad.1 = read.csv("Stability/GSE99517_ActD.HEK.HL.one.csv")
st.ad.2 = read.csv("Stability/GSE99517_ActD.HEK.HL.two.csv")
st.4s = st.4s.1 %>%
  select(geneName = X, hl1 = std) %>%
  inner_join(st.4s.2 %>% select(geneName = X, hl2 = std), by = "geneName") %>%
  mutate(hl.4s = (hl1 + hl2) / 2)
st.aa = st.aa.1 %>%
  select(geneName = X, hl1 = std) %>%
  inner_join(st.aa.2 %>% select(geneName = X, hl2 = std), by = "geneName") %>%
  mutate(hl.aa = (hl1 + hl2) / 2)
st.ad = st.ad.1 %>%
  select(geneName = X, hl1 = std) %>%
  inner_join(st.ad.2 %>% select(geneName = X, hl2 = std), by = "geneName") %>%
  mutate(hl.ad = (hl1 + hl2) / 2)
st.all = st.4s %>% select(geneName, hl.4s) %>%
  inner_join(st.ad %>% select(geneName, hl.ad), by = "geneName") %>%
  inner_join(st.aa %>% select(geneName, hl.aa), by = "geneName")

# Read and process stoat data
pro = read.table("Unique/PROseq.HEK.eRPKMgb.txt", header = 1)
rna = read.table("Unique/RNAseq.HEK.isoforms.txt", header = 1)
ted = read.table("Unique/TEDseq.HEK.rcRPKM.txt", header = 1)
uniq.genes = rna %>%
  group_by(geneName) %>%
  summarise(count = n()) %>%
  filter(count == 1) %>%
  ungroup() %>%
  select(geneName)

pro = pro %>%
  mutate(eRPKM = (eRPKMgb.1 + eRPKMgb.2)/2)
rna = rna %>%
  mutate(FPKM = (FPKM.DMSO1 + FPKM.DMSO2)/2)
ted = ted %>%
  mutate(RPKM = (RPKM.DMSO1 + RPKM.DMSO2)/2)

stoat = rna %>% select(id, geneName, FPKM) %>%
  inner_join(pro %>% select(id, eRPKM), by = "id") %>%
  inner_join(ted %>% select(id, RPKM), by = "id") %>%
  mutate(si.rna = FPKM/eRPKM, si.ted = RPKM/eRPKM) %>%
  filter_at(vars(FPKM, RPKM), all_vars(. > 0.001 & is.finite(.)))

stoat.m = stoat %>%
  group_by(geneName) %>%
  summarise(FPKM = mean(FPKM),
            RPKM = mean(RPKM),
            eRPKM = mean(eRPKM),
            si.rna = mean(si.rna),
            si.ted = mean(si.ted))

stoat.u = stoat %>%
  inner_join(uniq.genes)

# Correlation between RNA-seq and TED-seq in uniq genes
print("figS2b")
pdf("pdf/figS2b.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = stoat.u$FPKM, y = stoat.u$RPKM),
         expression(log[10]*RNA-seq),
         expression(log[10]*TED-seq),
         diag = T,
         xlim = c(-2, 3),
         ylim = c(-1, 4),
         red = 0,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# Correlation between RNA-seq and TED-seq in all genes
print("figS2a")
pdf("pdf/figS2a.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = stoat$FPKM, y = stoat$RPKM),
         expression(log[10]*RNA-seq),
         expression(log[10]*TED-seq),
         diag = T,
         xlim = c(-3, 4),
         ylim = c(-1, 5),
         red = 400,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# Correlation between RNA-seq and TED-seq in gene averages
print("fig2a")
pdf("pdf/fig2a.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = stoat.m$FPKM, y = stoat.m$RPKM),
         expression(log[10]*RNA-seq),
         expression(log[10]*TED-seq),
         diag = T,
         xlim = c(-2, 3),
         ylim = c(-1, 4),
         red = 400,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# Correlation between PRO-seq and RNA-seq in gene averages
print("fig2b")
pdf("pdf/fig2b.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = stoat.m$eRPKM, y = stoat.m$FPKM),
         expression(log[10]*PRO-seq),
         expression(log[10]*RNA-seq),
         diag = T,
         xlim = c(-3, 2),
         ylim = c(-2, 3),
         red = 400,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# Correlation between PRO-seq and TED-seq in gene averages
print("fig2c")
pdf("pdf/fig2c.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = stoat.m$eRPKM, y = stoat.m$RPKM),
         expression(log[10]*PRO-seq),
         expression(log[10]*TED-seq),
         diag = T,
         xlim = c(-3, 2),
         ylim = c(-1, 4),
         red = 400,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# compare between stability data and stoat
comp = stoat.m %>%
  inner_join(st.all, by = "geneName") %>%
  filter_at(vars(-geneName), all_vars(. > 0 & is.finite(.)))

# Correlation between stability index and half life data (modify these plots using ggplot scripts)
# Using TED-seq/PRO-seq
print("fig2d")
pdf("pdf/fig2d.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = comp$si.ted, y = comp$hl.4s),
         expression(log[10]*Stability~Index[TED-seq]),
         expression(log[10]*Half~Life[4*sU]),
         diag = T,
         #xlim = c(-1.5, 2.5),
         #ylim = c(-0.5, 1.5),
         xlim = c(-2, 3),
         ylim = c(-1, 1.5),
         red = 400,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# Correlation between stability index and half life data (modify these plots using ggplot scripts)
# Using RNA-seq/PRO-seq
print("fig2e")
pdf("pdf/fig2e.pdf", width = 3, height = 3)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = comp$si.rna, y = comp$hl.4s),
         expression(log[10]*Stability~Index[RNA-seq]),
         expression(log[10]*Half~Life[4*sU]),
         diag = T,
         xlim = c(-1.5, 2),
         ylim = c(-1, 1.5),
         red = 1000,
	 col = sl.col(11),
	 cor = T)
print(g)
dev.off()

# Scaled normalized data for limnear modeling
data = data.frame(
  hl = scale(log10(comp$hl.4s)),
  hl.aa = scale(log10(comp$hl.aa)),
  hl.ad = scale(log10(comp$hl.ad)),
  si.ted = scale(log10(comp$si.ted)),
  si.rna = scale(log10(comp$si.rna)),
  ted = scale(log10(comp$RPKM)),
  rna = scale(log10(comp$FPKM)),
  pro = scale(log10(comp$eRPKM)))

# Linear modeling
lm.rna = lm(hl ~ rna, data = data)
lm.ted = lm(hl ~ ted, data = data)
lm.pro = lm(hl ~ pro, data = data)
lm.si.rna = lm(hl ~ si.rna, data = data)
lm.si.ted = lm(hl ~ si.ted, data = data)

# merge correlation coefficients and standard errors
coef = bind_rows(data.frame(Model = "RNA", t(summary(lm.rna)$coefficients[2,])),
	data.frame(Model = "TED", t(summary(lm.ted)$coefficients[2,])),
	data.frame(Model = "PRO", t(abs(summary(lm.pro)$coefficients[2,]))),
	data.frame(Model = "RNA/PRO", t(summary(lm.si.rna)$coefficients[2,])),
	data.frame(Model = "TED/PRO", t(summary(lm.si.ted)$coefficients[2,]))) %>%
	mutate(Model = factor(Model, levels = c("RNA", "TED", "PRO", "RNA/PRO", "TED/PRO")))

# Plot barplot with error bars
print("fig2f")
pdf("fig2f.pdf", width = 3, height =3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g=ggplot(data = coef, aes(Model, Estimate)) +
	geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error), width = 0.4) +
	geom_col(aes(col = Model, fill = Model)) +
	scale_color_manual(values = sl.col(7,1,0.7)[2:6]) +
	scale_fill_manual(values = sl.col(7,0.9,1.1)[2:6]) +
	ylab("|Correlation coefficient|") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
	theme(legend.key = element.rect())
print(g)
dev.off()

