setwd("~/Work/labuser/sl2665/rev/fig/fig2/")
rm(ls = list())
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# Load miRNA and repression index data
source("rscript/fig2.repression.R")
pro.id = Reduce(inner_join,
  list(
    all %>%
      mutate(psi = (hek.psi + hela.psi)/2) %>%
      select(id, strand, psi, hek.psi, hela.psi),
    hek.pro.id %>% mutate(HEK = r1/trc['hek']) %>% select(id, HEK),
    hela.pro.id %>% mutate(HeLa = r1/trc['hela']) %>% select(id, HeLa),
    thp1.pro.id %>% mutate(THP1 = r1/trc['thp1']) %>% select(id, THP1)))

stoat.mir = pro.id %>%
  mutate(HEK = HEK * 10^psi,
         HeLa = HeLa * 10^psi,
         THP1 = THP1 * 10^psi) %>%
  select(-psi, -hela.psi, -hek.psi)

stoat.mirfam = stoat.mir %>%
  inner_join(mirfam, by = c("id", "strand")) %>%
  group_by(fam) %>%
  summarise(HEK = sum(HEK),
            HeLa = sum(HeLa),
            THP1 = sum(THP1)) %>%
  ungroup

stoat.RI = ts %>%
  inner_join(stoat.mirfam, by = "fam") %>%
  select(id, name, fam, HEK, HeLa, THP1, p) %>%
  group_by(id, name) %>%
  summarise(HEK = sum(log10(HEK+1) * p),
            HeLa = sum(log10(HeLa+1) * p),
            THP1 = sum(log10(THP1+1) * p)) %>%
  ungroup %>%
  inner_join(geneName %>% select(id, len, tlen), by = "id")

# Read PRO-seq data
pro.long = Reduce(bind_rows,
             lapply(c("HEK", "HeLa", "THP1"),
                    function(i) {
                      data.frame(sample = i,
                                 read.table(
                                   paste0("miRNA/readcount/PROseq.",
                                          i, ".txt"),
                                   header = T,
                                   stringsAsFactors = F)) %>%
                        return })) %>%
  select(sample, id, eRPKMgbs)
                        
# Read TED-seq data
ted.long = Reduce(bind_rows,
                  lapply(c("HEK.rep1", "HeLa", "THP1"),
                         function(i) {
                           data.frame(sample = i,
                                      read.table(
                                        paste0("data/TEDseq.",
                                               i,
                                               "/table/expression.txt"),
                                        header = T,
                                        stringsAsFactors = F)) %>%
                             return })) %>%
  select(sample, id, RPKMted) %>%
  mutate(sample = ifelse(sample == "HEK.rep1", "HEK", sample))

# Stability index
stoat.si = inner_join(pro.long,
                      ted.long,
                      by = c("sample", "id")) %>%
  mutate(si = log10(RPKMted/eRPKMgbs)) %>%
  filter(is.finite(si)) %>%
  inner_join(geneName %>% select(id, geneName, tlen), by = "id") %>%
#  group_by(sample, geneName) %>%
#  summarise(si = mean(si)) %>%
#  ungroup %>%
  select(sample, id, geneName, si) %>%
  spread(sample, si) %>%
  drop_na %>%
  mutate(name = geneName, si.HEK = HEK, si.HeLa = HeLa, si.THP1 = THP1) %>%
  select(-geneName, -HEK, -HeLa, -THP1)

# Merge stability index with Repression index
stoat.sri = stoat.RI %>%
  mutate(ri.HEK = HEK, ri.HeLa = HeLa, ri.THP1 = THP1) %>%
  select(-HEK, -HeLa, -THP1) %>%
  inner_join(stoat.si, by = c("id", "name")) %>%
  filter(ri.HEK > 0 & ri.HeLa > 0 & ri.THP1 > 0)

# Zscore of RI and SI
z.sri = stoat.sri %>%
  mutate_at(vars(contains("si")), funs(scale))

# Correlation between STOAT miRNA repression index and STOAT SI in HEK, THP1, 
library(gridExtra)
pdf("pdf/figS2n.pdf", width = 8, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g1=plot.cor(data.frame(x = stoat.sri$ri.HEK , y = stoat.sri$si.HEK),
           "miRNA RI (STOAT)",
           expression(log[10]~"SI (STOAT)"),
           lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
           xlim=c(0,20), ylim=c(-2, 3)) + ggtitle("HEK293")
g2=plot.cor(data.frame(x = stoat.sri$ri.HeLa , y = stoat.sri$si.HeLa),
            "miRNA RI (STOAT)",
            expression(log[10]~"SI (STOAT)"),
            lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
            xlim=c(0,20), ylim=c(-2, 3)) + ggtitle("HeLa")
g3=plot.cor(data.frame(x = stoat.sri$ri.THP1 , y = stoat.sri$si.THP1),
            "miRNA RI (STOAT)",
            expression(log[10]~"SI (STOAT)"),
            lin = T, diag = T, cor = T, fit = "l", col = sl.col(11),
            xlim=c(0,20), ylim=c(-2, 3)) + ggtitle("THP1")
grid.arrange(g1, g2, g3, ncol = 3)
dev.off()

# Violin plots of STOAT miRNA repression index
sl.col = function(n = 7, sat = 1, lum = 0.85) {
  col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(n) %>%
    col2rgb %>%
    rgb2hsv
  col[2,] = col[2,] * sat
  col[3,] = col[3,] * lum
  col[3,col[3,]>1] = 1
  return(hsv(col[1,], col[2,], col[3,])) }

sl.fill = rev(sl.col(20,0.9,1.1))[c(18,5,8,16,14,10,12,20)]
sl.cols = rev(sl.col(20,1,0.7))[c(18,5,8,16,14,10,12,20)]

st.v = stoat.sri %>%
  select(id, ri.HEK, ri.HeLa, ri.THP1) %>%
  gather(sample, mRI, -id) %>%
  separate(sample, c("pre", "cell"), sep = "\\.") %>%
  select(-pre) %>%
  inner_join(stoat.sri %>%
               select(id, si.HEK, si.HeLa, si.THP1) %>%
               gather(sample, SI, -id) %>%
               separate(sample, c("pre", "cell"), sep = "\\.") %>%
               select(-pre),
             by = c("id", "cell")) %>%
  mutate(cell = ifelse(cell == "HEK", "HEK293",cell)) %>%
  mutate(Type = ifelse(mRI > 10, "mRI > 10", "mRI < 10")) %>%
  mutate(Type = factor(Type, levels = c("mRI < 10", "mRI > 10")))

# Calculate p-values
fig2o.pval = c(HEK293 = t.test(st.v %>% filter(cell == "HEK293" & mRI > 10) %>% select(SI) %>% unlist,
                               st.v %>% filter(cell == "HEK293" & mRI < 10) %>% select(SI) %>% unlist)$p.value,
               HeLa = t.test(st.v %>% filter(cell == "HeLa" & mRI > 10) %>% select(SI) %>% unlist,
                             st.v %>% filter(cell == "HeLa" & mRI < 10) %>% select(SI) %>% unlist)$p.value,
               THP1 = t.test(st.v %>% filter(cell == "THP1" & mRI > 10) %>% select(SI) %>% unlist,
                             st.v %>% filter(cell == "THP1" & mRI < 10) %>% select(SI) %>% unlist)$p.value)
print(fig2o.pval)

pdf("pdf/fig2o.pdf", width = 3.5, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g = ggplot(st.v, aes(x = Type, y = SI)) +
  geom_violin(aes(col = Type, fill = Type), width = 0.9) +
  geom_boxplot(aes(col = Type), width = 0.4, outlier.shape = NA) +
  theme_bw() +
  scale_color_manual(values = sl.cols) +
  scale_fill_manual(values = sl.fill) +
  ylab(expression(log[10]~"Stability Index")) +
  xlab("miRNA Repression Index (mRI)") +
  facet_grid(. ~ cell, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.key = element_rect()) +
  labs(col = "mRI group", fill = "mRI group")
print(g)
dev.off()


# Differentially miRNA control
thp1.diff = z.sri %>% 
  mutate(dRI = ri.THP1 - ri.HEK, dSI = (si.THP1 - si.HEK) * log(10) / log(2)) %>%
  select(id, name, len, tlen, ri.HEK, ri.THP1, dRI, dSI) %>%
  mutate(Type = ifelse((ri.HEK > 10 & ri.THP1 < 10) & dRI < -7, "Lower mRI\nin THP1","No\ndifference")) %>%
  mutate(Type = factor(Type, levels = c("No\ndifference", "Lower mRI\nin THP1"))) %>%
  arrange(dRI) 

# Calculate p-values
fig2p.pval  = t.test(thp1.diff %>% filter(Type == "No\ndifference") %>% select(dSI) %>% unlist,
                     thp1.diff %>% filter(Type == "Lower mRI\nin THP1") %>% select(dSI) %>% unlist)$p.value
                     
fig2p.pval

pdf("pdf/fig2p.pdf", width = 2.5, height = 3)
par(mar=c(2,2,1,1), mgp=c(1.5,0.4,0))
g = ggplot(thp1.diff, aes(x = Type, y = dSI)) +
  geom_violin(aes(col = Type, fill = Type), width = 0.9) +
  geom_boxplot(aes(col = Type), width = 0.4, outlier.shape = NA) +
  theme_bw() +
  scale_color_manual(values = sl.cols) +
  scale_fill_manual(values = sl.fill) +
  ylab(expression(log[2]~"SI (THP1 / HEK293)")) +
  xlab("miRNA Repression Index (mRI)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5)) +
  theme(legend.key = element_rect()) +
  labs(col = expression(italic(Delta)*"mRI group"),
       fill = expression(italic(Delta)*"mRI group"))
print(g)
dev.off()

# Track an example RNA
transcript = "FEM1C"

mirfam.net = stoat.mirfam %>%
  filter(fam %in% (ts %>%
           filter(name == transcript) %>%
           filter(p > 0.5) %>%
          select(fam) %>% unlist)) %>%
  mutate_if(is.numeric, funs(log10(. + 1)))

mir.net = mirfam %>%
  inner_join(mirfam.net %>% select(fam),
             by= "fam") %>%
  inner_join(stoat.mir,
             by = c("id", "strand"))

stoat.net = stoat.si %>%
  filter(name == transcript)

stoat.raw.net = inner_join(
  pro.long %>% filter(id %in% (stoat.net %>% select(id) %>% unlist)),
  ted.long %>% filter(id %in% (stoat.net %>% select(id) %>% unlist)),
  by = c("id", "sample"))
