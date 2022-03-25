# Pre-run fig3.R before running this script 

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

mat.common = mat.all %>% filter(id %in% mat.all.id$id)
hek = mat.common %>% filter(sample == "HEK")
hela = mat.common %>% filter(sample == "HeLa")
thp = mat.common %>% filter(sample == "THP1")

lm.hek.au = lm(pal ~ au, data = hek)
lm.hek.len = lm(pal ~ len, data = hek)
lm.hela.au = lm(pal ~ au, data = hela)
lm.hela.len = lm(pal ~ len, data = hela)
lm.thp.au = lm(pal ~ au, data = thp)
lm.thp.len = lm(pal ~ len, data = thp)

# merge correlation coefficients and standard errors
coef = bind_rows(data.frame(Model = "AU", Sample = "HEK", t(summary(lm.hek.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "HEK", t(summary(lm.hek.len)$coefficients[2,])),
	data.frame(Model = "AU", Sample = "HeLa", t(summary(lm.hela.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "HeLa", t(summary(lm.hela.len)$coefficients[2,])),
	data.frame(Model = "AU", Sample = "THP1", t(summary(lm.thp.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "THP1", t(summary(lm.thp.len)$coefficients[2,])))

# Plot barplot with error bars
pdf("pdf/fig3i.pdf", width = 3, height =3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g = ggplot(data = coef, aes(Sample, Estimate, fill = Model)) +
	geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error),
		width = 0.2, position = position_dodge(width=0.9)) +
	geom_col(aes(col=Model), position = "dodge", width = 0.8) +
	scale_color_manual(values = sl.cols) +
	scale_fill_manual(values = sl.fill) +
	ylab("Correlation coefficient") +
	scale_x_discrete(labels = c("HEK293", "HeLa", "THP1"))+
	theme_bw() +
	theme(legend.key = element_rect())
print(g)
dev.off()

# Plot barplot for HEK replicates
mat.comrep = mat.allreps %>% filter(id %in% mat.allreps.id$id)
hek.a = mat.comrep %>% filter(sample == "HEK")
hek3.a = mat.comrep %>% filter(sample == "HEK2")
hela.a = mat.comrep %>% filter(sample == "HeLa")
thp.a = mat.comrep %>% filter(sample == "THP1")

lm2.hek.au = lm(pal ~ au, data = hek.a)
lm2.hek.len = lm(pal ~ len, data = hek.a)
lm2.hek3.au = lm(pal ~ au, data = hek3.a)
lm2.hek3.len = lm(pal ~ len, data = hek3.a)
lm2.hela.au = lm(pal ~ au, data = hela.a)
lm2.hela.len = lm(pal ~ len, data = hela.a)
lm2.thp.au = lm(pal ~ au, data = thp.a)
lm2.thp.len = lm(pal ~ len, data = thp.a)

# merge correlation coefficients and standard errors
coef2 = bind_rows(data.frame(Model = "AU", Sample = "HEK", t(summary(lm2.hek.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "HEK", t(summary(lm2.hek.len)$coefficients[2,])),
	data.frame(Model = "AU", Sample = "HEK2", t(summary(lm2.hek3.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "HEK2", t(summary(lm2.hek3.len)$coefficients[2,])),
	data.frame(Model = "AU", Sample = "HeLa", t(summary(lm2.hela.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "HeLa", t(summary(lm2.hela.len)$coefficients[2,])),
	data.frame(Model = "AU", Sample = "THP1", t(summary(lm2.thp.au)$coefficients[2,])),
	data.frame(Model = "Len", Sample = "THP1", t(summary(lm2.thp.len)$coefficients[2,])))

# Plot barplot with error bars
pdf("pdf/figS3i.pdf", width = 3, height =3)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
g = ggplot(data = coef2, aes(Sample, Estimate, fill = Model)) +
	geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error),
		width = 0.2, position = position_dodge(width=0.9)) +
	geom_col(aes(col=Model), position = "dodge", width = 0.8) +
	scale_color_manual(values = sl.cols) +
	scale_fill_manual(values = sl.fill) +
	ylab("Correlation coefficient") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(legend.key = element_rect())
print(g)
dev.off()
