setwd("~/Work/labuser/sl2665/rev/fig/fig1/fig1i/")
rm(list = ls())
library(dplyr)
library(tidyr)

sample = c("HEK1",
	   "HeLa")
# Read PAL and expression tables
pal.list = lapply(paste0("table/", sample, ".txt"), read.table,
		  header = F, stringsAsFactors = F, col.names = c("id", "pal"))
exp.list = lapply(paste0("table/", sample, ".exp.txt"), read.table, header = T, stringsAsFactors = F)

# conjoin PAL and expression tables
pal = inner_join(pal.list[[1]] %>% mutate(HEK1 = pal) %>% select(-pal),
		 pal.list[[2]] %>% mutate(HeLa = pal) %>% select(-pal))
exp = inner_join(exp.list[[1]] %>% mutate(HEK1.exp = ted) %>% select(id, HEK1.exp),
		 exp.list[[2]] %>% mutate(HeLa.exp = ted) %>% select(id, HeLa.exp))
tab = inner_join(pal, exp)

# Read gene id and nama table
id = read.table("table/transcripts.txt", header = T, stringsAsFactors = F)
tab = inner_join(id, tab)

# PAL difference table sorted by PAL difference and sort/filter for the most different, highly expressed
tab.diff = tab %>%
	mutate(pal.diff = HeLa - HEK1) %>%
	arrange(pal.diff) %>%
	filter(abs(pal.diff) > 40) %>%
	filter(HEK1.exp > 200 & HeLa.exp > 200)
write.table(tab.diff, "PAL_difference.txt", quote = F, row.names = F, sep = "\t")
