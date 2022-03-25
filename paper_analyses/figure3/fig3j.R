# AU rich RNA binding protein expression levels

ARElist = data.frame(
  ABP = c("AUF1", "TTP", "ZFP36L1", "ZFP36L2",
          "KSRP", "HuR", "HuD", "GAPDH", "LDHM"),
  ENST = paste0("ENST00000", c(313899, 594442, 282388, 336440, 398148,
           407627, 371827, 229239, 422447))
)

ABP.exp = TEDseq.exp.long %>%
  separate(id, c("ENST", "isoform"), sep = "\\.") %>%
  inner_join(ARElist, by = "ENST") %>%
  select(ABP, sample, val) %>%
  filter(sample != "HEK2" & ABP %in% c("HuR", "TTP", "ZFP36L1", "ZFP36L2"))

library(RColorBrewer)
library(ggplot2)

sl.col = function(n = 7, sat = 1, lum = 0.85) {
  col = colorRampPalette((brewer.pal(11, "Spectral")))(n) %>%
    col2rgb %>%
    rgb2hsv
  col[2,] = col[2,] * sat
  col[3,] = col[3,] * lum
  col[3,col[3,]>1] = 1
  return(hsv(col[1,], col[2,], col[3,])) }

sl.fill = sl.col(20,0.9,1.1)[c(18,5,8,16,14,10,12,20)]
sl.cols = sl.col(20,1,0.7)[c(18,5,8,16,14,10,12,20)]

# Plot barplot
pdf("pdf/fig3j.pdf", width = 3, height = 3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g = ggplot(data = ABP.exp %>% filter(ABP %in% c("HuR", "TTP")), aes(sample, val, fill = sample)) +
  geom_col(aes(col = sample), position = "dodge", width = 0.8) +
  scale_color_manual(values = sl.cols) +
  scale_fill_manual(values = sl.fill) +
  ylab("Normalized RPM") +
  xlab("Cell type") +
  scale_x_discrete(labels = c("HEK293", "HeLa", "THP1"))+
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~ABP, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
print(g)
dev.off()

