setwd("~/Work/labuser/sl2665/rev/fig/fig2/")
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Load miRNA expression data in HEK and HeLa
hek.tot = Reduce(bind_rows, lapply(read.table("miRNA/HEK/samples.txt", header = F)[,1],
	function(a) data.frame(sample = a,
							read.table(paste0("miRNA/HEK/",a,".txt"),
										header = F,
										col.names=c("id","val"))))) %>%
	group_by(id) %>%
	summarise(mean = mean(val)) %>%
	ungroup %>%
	separate(id, c("spe","mir","name","strand"), sep = "-") %>%
	unite(id, mir, name, sep = "-") %>%
	select(-spe) %>%
	mutate(strand = ifelse(is.na(strand), "5p", strand))

hela.tot = Reduce(bind_rows, lapply(read.table("miRNA/HeLa/samples.txt", header = F)[,1],
	function(a) data.frame(sample = a,
							read.table(paste0("miRNA/HeLa/",a,".txt"),
										header = F,
										col.names=c("id","val"))))) %>%
	group_by(id) %>%
	summarise(mean = mean(val)) %>%
	ungroup %>%
	separate(id, c("spe","mir","name","strand"), sep = "-") %>%
	unite(id, mir, name, sep = "-") %>%
	select(-spe) %>%
	mutate(strand = ifelse(is.na(strand), "5p", strand))

miRNA.tot = inner_join(
	hek.tot %>% mutate(HEK293 = mean) %>% select(-mean),
	hela.tot %>% mutate(HeLa = mean) %>% select(-mean),
	by = c("id", "strand"))

# Load PRO-seq estimated miRNA expression data
# Read PRO-seq derived data
hek.pro = read.table("miRNA/readcount/HEK.count.txt", header = F, col.names = c("id","name","r1","r2","r3"))
hela.pro = read.table("miRNA/readcount/HeLa.count.txt", header = F, col.names = c("id","name","r1","r2","r3"))
thp1.pro = read.table("miRNA/readcount/THP1.count.txt", header = F, col.names = c("id","name","r1","r2","r3"))

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
  inner_join(miRNA.tot, by = "id") %>%
  mutate(hek.pro = r1/trc['hek'])
hela.all = hela.pro.id %>%
  inner_join(miRNA.tot, by = "id") %>%
  mutate(hela.pro = r1/trc['hela'])
thp1.all = thp1.pro.id %>%
  mutate(thp.pro = r1/trc['thp1'])
all = inner_join(
  hek.all %>% select(id, strand, HEK293, HeLa, hek.pro),
  hela.all %>% select(id, strand, hela.pro),
  by = c("id", "strand")) %>%
  mutate(pro.ratio = hek.pro/hela.pro,
         hek.psi = log10(HEK293/hek.pro),
         hela.psi = log10(HeLa/hela.pro)) %>%
  filter_at(vars(-id, -strand), all_vars(is.finite(.)))

# estimate HEK smRNA-seq using HeLa PSI
all.est = all %>%
  mutate(hek.est = hek.pro * 10^hela.psi,
         hela.est = hela.pro * 10^hek.psi) %>%
  select(id, strand, hek.est, hela.est)

# Load targetscan data
ts = read.table("miRNA/targetscan.txt", header = T, sep = "\t") %>%
	mutate(p = as.numeric(as.character(Aggregate.PCT))) %>%
	mutate(id = Transcript.ID, name = Gene.Symbol, fam = miRNA.family,
		mir = Representative.miRNA) %>%
	select(id, name, fam, mir, p) %>%
	filter(is.finite(p))

# miRNA family info
mirfam = read.table("miRNA/mirfam.txt", header = T) %>%
	separate(name, c("spe","mir","n","strand"), sep = "-") %>%
	unite(id, mir, n, sep = "-") %>%
	select(-spe) %>%
	mutate(strand = ifelse(is.na(strand), "5p", strand))

# sum miRNA expression level by mir family
hek.fam = hek.tot %>%
	inner_join(mirfam, by = c("id", "strand")) %>%
	group_by(fam) %>%
	summarise(expr = sum(mean)) %>%
	ungroup

hela.fam = hela.tot %>%
	inner_join(mirfam, by = c("id", "strand")) %>%
	group_by(fam) %>%
	summarise(expr = sum(mean)) %>%
	ungroup

# miRNA family expression level by PRO-seq
all.pfam = all.est %>%
  inner_join(mirfam, by = c("id", "strand")) %>%
  group_by(fam) %>%
  summarise(hek.pfam = sum(hek.est),
            hela.pfam = sum(hela.est)) %>%
  ungroup

# tag targetscan result with miRNA expression level for the target mRNAs
hek.ts.expr = ts %>%
	inner_join(hek.fam, by = "fam") %>%
	select(id, name, fam, expr, p)

hela.ts.expr = ts %>%
	inner_join(hela.fam, by = "fam") %>%
	select(id, name, fam, expr, p)

all.ts.expr = ts %>%
  inner_join(all.pfam, by = "fam") %>%
  select(id, name, fam, hek.pfam, hela.pfam, p)

# Calculate miRNA repression index 
thresfx = function(x) {
	log10(x+1) %>%
		return
}
hek.ts.ri = hek.ts.expr %>%
	group_by(id, name) %>%
	summarise(RI = sum(thresfx(expr) * p)) %>%
	ungroup %>%
	arrange(-RI)

hela.ts.ri = hela.ts.expr %>%
	group_by(id, name) %>%
	summarise(RI = sum(thresfx(expr) * p)) %>%
	ungroup %>%
	arrange(-RI)

ri.pfam.all = all.ts.expr %>%
  group_by(id, name) %>%
  summarise(hek.RI = sum(thresfx(hek.pfam) * p),
            hela.RI = sum(thresfx(hela.pfam) * p)) %>%
  ungroup %>%
  mutate(dRI = log2(hek.RI/hela.RI)) %>%
  filter(is.finite(dRI)) %>%
  arrange(-abs(dRI))

ri.all = inner_join(hek.ts.ri %>% mutate(hek = RI) %>% select(id, name, hek),
	hela.ts.ri %>% mutate(hela = RI) %>% select(id, name, hela),
	by = c("id", "name")) %>%
	mutate(dRI = log2(hek/hela)) %>%
	filter(is.finite(dRI)) %>%
	arrange(-abs(dRI))

# Stability data
st.4s.1 = read.csv("Stability/GSE99517_HEK.4SU.HL.DRUID.one.csv")
st.4s.2 = read.csv("Stability/GSE99517_HEK.4SU.HL.DRUID.two.csv")
st.4s = st.4s.1 %>%
  select(geneName = X, hl1 = std) %>%
  inner_join(st.4s.2 %>% select(geneName = X, hl2 = std), by = "geneName") %>%
  mutate(hl.4s = (hl1 + hl2) / 2)

# RI-stability
ris = inner_join(ri.all,
		st.4s %>% mutate(name = geneName) %>% select(name, hl.4s),
		by = "name") %>%
	filter_if(vars(is.numeric(.)), all_vars(is.finite(.)))

# Read and process stoat data
pro = read.table("miRNA/readcount/PROseq.HEK.txt", header = 1)
ted = read.table("Unique/TEDseq.HEK.rcRPKM.txt", header = 1)

pro = pro %>%
  mutate(eRPKM = eRPKMgbs)
ted = ted %>%
  mutate(RPKM = (RPKM.DMSO1 + RPKM.DMSO2)/2)

stoat = pro %>% select(id, eRPKM) %>%
  inner_join(ted %>% select(id, RPKM), by = "id") %>%
  mutate(si = RPKM/eRPKM) %>%
  filter(eRPKM > 0.1)

geneName = read.table("bed/gencode.v26.geneName.bed13", header = F)
geneName$len = geneName$V3-geneName$V2
geneName$tlen = sapply(geneName$V11, function(a) {sum(as.numeric(strsplit(as.character(a), ",")[[1]]))})
geneName = geneName[c(4,13,14,15)]
colnames(geneName) = c("id", "geneName", "len", "tlen")
geneName = geneName %>%
  inner_join( read.table("miRNA/gencode.v26.3utrlen.txt", header = F, col.names = c("id", "u3len")), by = "id") 

# Final stoat data with RI calculated from smRNA-seq
stoat.m = geneName %>%
  inner_join(stoat, by = "id") %>%
  group_by(geneName) %>%
  summarise(RPKM = mean(RPKM),
            eRPKM = mean(eRPKM),
            si = mean(si),
            len = mean(len/1000),
            tlen = mean(tlen/1000),
            u3len = mean(u3len/1000)) %>%
  ungroup %>%
  inner_join(ris %>% mutate(geneName = name) %>%
               select(geneName, hek, hl.4s),
             by = "geneName")

stoat.ms = stoat.m %>%
  inner_join(ri.pfam.all %>% mutate(geneName = name) %>% select(geneName, hek.RI),
             by = "geneName")

source("rscript/scatterPlot.R")
sl.col = function(n = 7, sat = 1, lum = 0.85) {
  col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(n) %>%
    col2rgb %>%
    rgb2hsv
  col[2,] = col[2,] * sat
  col[3,] = col[3,] * lum
  col[3,col[3,]>1] = 1
  return(hsv(col[1,], col[2,], col[3,])) }

# Correlation between miRNA repression index and 4sU stability
pdf("pdf/fig2k.pdf", width = 3, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g=plot.cor(data.frame(x = stoat.m$hek , y = log10(stoat.m$hl.4s)),
           "miRNA repression index (smRNA-seq)",
           expression(log[10]~"4sU half-life"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim=c(0,30), ylim=c(-0.5, 1.5))
print(g)
dev.off()

# Control correlation between gene length and 4sU stability
pdf("pdf/figS2k.pdf", width = 3, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g=plot.cor(data.frame(x = (stoat.m$u3len) , y = log10(stoat.m$hl.4s)),
           expression(3*minute~"UTR length (kb)"),
           expression(log[10]~"4sU half-life"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim = c(0,5), ylim=c(-0.5, 1.5))
print(g)
dev.off()

# Correlation between miRNA repression index and stability index
pdf("pdf/fig2l.pdf", width = 3, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g=plot.cor(data.frame(x = stoat.m$hek , y = log10(stoat.m$si)),
           "miRNA repression index (smRNA-seq)",
           expression(log[10]~"Stability Index (STOAT)"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim=c(0,30), ylim=c(-2, 3))
print(g)
dev.off()

# Control correlation between gene length and stability index
pdf("pdf/figS2l.pdf", width = 3, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g=plot.cor(data.frame(x = stoat.m$u3len , y = log10(stoat.m$si)),
           expression(log[10]~"Transcript length (kb)"),
           expression(log[10]~"Stability Index (STOAT)"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim=c(0,10), ylim=c(-2, 3))
print(g)
dev.off()

# Correlation between STOAT miRNA repression index and 4sU stability
pdf("pdf/fig2m.pdf", width = 3, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g=plot.cor(data.frame(x = stoat.ms$hek.RI , y = log10(stoat.ms$hl.4s)),
           "miRNA repression index (STOAT)",
           expression(log[10]~"4sU half-life"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim=c(0,30), ylim=c(-0.5, 1.5))
print(g)
dev.off()

# Correlation between STOAT miRNA repression index and STOAT SI
pdf("pdf/fig2n.pdf", width = 3, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g=plot.cor(data.frame(x = stoat.ms$hek.RI , y = log10(stoat.ms$si)),
           "miRNA repression index (STOAT)",
           expression(log[10]~"Stability Index (STOAT)"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim=c(0,30), ylim=c(-2, 3))
print(g)
dev.off()

