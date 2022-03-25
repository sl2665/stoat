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

hek.list = Reduce(bind_rows, lapply(read.table("miRNA/HEK/samples.txt", header = F)[,1],
	function(a) data.frame(sample = a,
							read.table(paste0("miRNA/HEK/",a,".txt"),
										header = F,
										col.names=c("id","val")))))

hela.list = Reduce(bind_rows, lapply(read.table("miRNA/HeLa/samples.txt", header = F)[,1],
	function(a) data.frame(sample = a,
							read.table(paste0("miRNA/HeLa/",a,".txt"),
										header = F,
										col.names=c("id","val")))))

hek.rep = hek.list %>%
	spread(sample, val) %>%
	separate(id, c("spe","mir","name","strand"), sep = "-") %>%
	unite(id, mir, name, sep = "-") %>%
	select(-spe) %>%
	mutate(strand = ifelse(is.na(strand), "5p", strand)) %>%
	mutate_if(is.numeric, log10) %>%
	filter_if(is.numeric, is.finite)

hela.rep = hela.list %>%
	spread(sample, val) %>%
	separate(id, c("spe","mir","name","strand"), sep = "-") %>%
	unite(id, mir, name, sep = "-") %>%
	select(-spe) %>%
	mutate(strand = ifelse(is.na(strand), "5p", strand)) %>%
	mutate_if(is.numeric, log10) %>%
	filter_if(is.numeric, is.finite)

hek.mean = hek.list %>%
	group_by(sample) %>%
	summarise(mean = 2^mean(log2(val))) %>%
	ungroup

hela.mean = hela.list %>%
	group_by(sample) %>%
	summarise(mean = 2^mean(log2(val))) %>%
	ungroup

hek.norm = hek.list %>%
	inner_join(hek.mean, by = "sample") %>%
	mutate(val = val / mean) %>%
	select(-mean)

hela.norm = hela.list %>%
	inner_join(hela.mean, by = "sample") %>%
	mutate(val = val / mean) %>%
	select(-mean)

hek.sum = hek.norm %>%
	group_by(id) %>%
	summarise(mean = mean(val),
		sd = sd(val)) %>%
	ungroup()

hela.sum = hela.norm %>%
	group_by(id) %>%
	summarise(mean = mean(val),
		sd = sd(val)) %>%
	ungroup()

hek.tot = hek.list %>%
	group_by(id) %>%
	summarise(mean = mean(val),
		sd = sd(val)) %>%
	ungroup()

hela.tot = hela.list %>%
	group_by(id) %>%
	summarise(mean = mean(val),
		sd = sd(val)) %>%
	ungroup()

miRNA.all = inner_join(
	hek.sum %>% mutate(HEK293 = log10(mean)) %>% select(-sd, -mean),
	hela.sum %>% mutate(HeLa = log10(mean)) %>% select(-sd, -mean),
	by = "id") %>%
	mutate(ratio = HEK293 - HeLa) %>%
	arrange(ratio)

miRNA.tot = inner_join(
	hek.tot %>% mutate(HEK293 = log10(mean)) %>% select(-sd, -mean),
	hela.tot %>% mutate(HeLa = log10(mean)) %>% select(-sd, -mean),
	by = "id") %>%
	mutate(ratio = HEK293 - HeLa) %>%
	arrange(ratio)

# Read PRO-seq derived data
hek.pro = read.table("miRNA/readcount/HEK.count.txt", header = F, col.names = c("id","name","r1","r2","r3"))
hela.pro = read.table("miRNA/readcount/HeLa.count.txt", header = F, col.names = c("id","name","r1","r2","r3"))
thp1.pro = read.table("miRNA/readcount/THP1.count.txt", header = F, col.names = c("id","name","r1","r2","r3"))

# Convert miRNA names
miRNA.id = miRNA.tot %>%
	separate(id, c("spe","mir","name","strand"), sep = "-") %>%
	unite(id, mir, name, sep = "-") %>%
	select(-spe)

pro.id = function(a) {
	a %>%
	mutate(name = tolower(substring(name, 4))) %>%
	mutate(name = ifelse(substring(name, 1, 3)=="let",
		paste0("let-", substring(name, 4)),
		paste0("miR-", name))) %>%
	separate(name, c("mir", "rest"), sep = "-", remove = F, extra = "merge") %>%
	mutate(rest = sub("\\-[0-9]+", "", rest)) %>%
	mutate(rest = sub("([a-z]+)","\\1_", rest)) %>%
	mutate(rest = sub("_[0-9]*","",rest)) %>%
	unite(id, mir, rest, sep = "-") %>%
	group_by(id) %>%
	summarise(r1 = sum(r1), r2 = sum(r2), r3 = sum(r3)) %>%
	return
}

hek.pro.id = hek.pro %>% pro.id
hela.pro.id = hela.pro %>% pro.id
thp1.pro.id = thp1.pro %>% pro.id
trc = c(hek = 58.564786,
	hela = 6.662422,
	thp1 = 3.006169)

# Merge all data after RPM normalization
hek.all = hek.pro.id %>%
	inner_join(miRNA.id, by = "id") %>%
	mutate(hek.pro = log10(r1/trc['hek']))
hela.all = hela.pro.id %>%
	inner_join(miRNA.id, by = "id") %>%
	mutate(hela.pro = log10(r1/trc['hela']))
thp1.all = thp1.pro.id %>%
	mutate(thp.pro = log10(r1/trc['thp1']))
all = inner_join(
	hek.all %>% select(id, strand, HEK293, HeLa, ratio, hek.pro),
	hela.all %>% select(id, strand, hela.pro),
	by = c("id", "strand")) %>%
	mutate(pro.ratio = hek.pro - hela.pro,
		hek.psi = HEK293 - hek.pro,
		hela.psi = HeLa - hela.pro) %>%
	filter_at(vars(-id, -strand), all_vars(is.finite(.)))

# estimate HEK smRNA-seq using HeLa PSI
all.est = all %>%
	mutate(hek.est = hek.pro + hela.psi,
		hela.est = hela.pro + hek.psi)

# Make plots
source("rscript/scatterPlot.R")
# PRO-seq vs smRNA-seq
pdf("pdf/fig2g.pdf", width = 2.5, height = 2.5)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = hek.all$hek.pro, y = hek.all$HEK293),
    expression(log[10]~"PRO-seq (RPKM)"),
    expression(log[10]~"smRNA-seq (RPM)"),
    diag = T,
    xlim = c(-2.5, 1.5),
    ylim = c(0.5, 5),
    red = 100,
	lin = T,
	col = sl.col(11)[1],
	cor = T)
print(g)
dev.off()

# Processing-stability index HEK293 vs HeLa
pdf("pdf/fig2h.pdf", width = 2.5, height = 2.5)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = all$hek.psi, y = all$hela.psi),
    expression(log[10]~"PSI (HEK293)"),
    expression(log[10]~"PSI (HeLa)"),
    diag = T,
    xlim = c(0, 5),
    ylim = c(0, 5),
    red = 100,
	lin = T,
	col = sl.col(11)[3],
	cor = T)
print(g)
dev.off()

# Estimated smRNA-seq 
pdf("pdf/fig2i.pdf", width = 2.5, height = 2.5)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = all.est$hek.est, y = all.est$HEK293),
    expression(log[10]~"PRO-seq estimated miRNA"),
    expression(log[10]~"smRNA-seq (RPM)"),
    diag = T,
    xlim = c(0.5, 5),
    ylim = c(0.5, 5),
    red = 100,
	lin = T,
	col = sl.col(11)[8],
	cor = T)
print(g)
dev.off()

# Fold change correlation
pdf("pdf/fig2j.pdf", width = 2.5, height = 2.5)
par(mar = c(2,2,1,1), mgp = c(1.5, 0.4, 0))
g=plot.cor(data.frame(x = all$pro.ratio, y = all.est$ratio),
    expression(italic(Delta)~log[10]~"PRO-seq (HEK293/HeLa)"),
    expression(italic(Delta)~log[10]~"smRNA-seq (HEK293/HeLa)"),
    diag = T,
    xlim = c(-1.5, 1.5),
    ylim = c(-2, 2),
    red = 100,
	lin = T,
	col = sl.col(11)[11],
	cor = T)
print(g)
dev.off()


