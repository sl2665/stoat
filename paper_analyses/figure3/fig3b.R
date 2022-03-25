# Pre-run fig3.R before running this script 
library(ggplot2)
library(dplyr)
library(tidyr)

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

pdf("pdf/fig3b.pdf", width = 5, height = 2)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))

gg = ggplot(TED.long %>% filter(sample !="HEK2" & type == "pal") %>%
	    mutate(sample = ifelse(sample == "HEK", "HEK293", sample))) +
    xlim(0,200) +
    xlab(expression(Median~poly*(A)~length)) +
    ylab("Transcript count") +
    theme_bw() +
    scale_color_manual(values = sl.cols) +
    scale_fill_manual(values = sl.fill) +
    theme(legend.position = "none")
  
print(gg + 
  geom_histogram(aes(x = val, fill = sample), breaks = seq(-10, 300, by = 5)) +
  stat_bin(aes(x = val - 2.5, col = sample), breaks = seq(-10, 300, by = 5) - 2.5, geom = "step", size = 0.4) +
  facet_wrap(~ sample, ncol = 3, scales = "free_y")
  )

dev.off()

pdf("pdf/figS3b.pdf", width = 8, height = 2)
par(mar = c(2, 2, 1, 1), mgp = c(1.5, 0.5, 0))
gg = ggplot(TED.long %>% filter(type == "pal") %>%
	    mutate(sample = ifelse(sample == "HEK", "HEK1", sample))) +
    xlim(0,250) +
    xlab(expression(Median~poly*(A)~length)) +
    ylab("Transcript count") +
    theme_bw() +
    scale_color_manual(values = sl.cols[c(1,4,2,3)]) +
    scale_fill_manual(values = sl.fill[c(1,4,2,3)]) +
    theme(legend.position = "none")
  
print(gg + 
  geom_histogram(aes(x = val, fill = sample), breaks = seq(-10, 300, by = 5)) +
  stat_bin(aes(x = val - 2.5, col = sample), breaks = seq(-10, 300, by = 5) - 2.5, geom = "step", size = 0.4) +
  facet_wrap(~ sample, ncol = 4, scales = "free_y")
  )
dev.off()
