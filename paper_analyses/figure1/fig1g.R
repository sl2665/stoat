setwd("~/Work/labuser/sl2665/rev/fig/fig1/fig1g/")
rm(list = ls())

library(dplyr)
library(tidyr)
library(DESeq2)

# Sample list
samples = c("HEK1",
	    "HEK2",
	    "HeLa",
	    "THP1")

# Load PRO-seq expression tables and merge to one table. Use promoter reads
pro.list = lapply(paste0("table/", samples, ".txt"), read.table, header = T, stringsAsFactors = F)
pro = lapply(1:length(samples), function(i) pro.list[[i]] %>%
	     select(id, pp) %>% mutate(!!samples[i] := pp) %>%
	     select(-pp))
pro.mg = Reduce(inner_join, pro)

# Convert merged table to a matrix comaptible with DESeq2
pro.mat = pro.mg[,-1]
rownames(pro.mat) = pro.mg$id
sample.info = data.frame(sample = c("HEK", "HEK", "HeLa", "THP1"))
rownames(sample.info) = colnames(pro.mat)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = pro.mat,
			      colData = sample.info,
			      design = ~ sample)
dds <- DESeq(dds)

# Find differentially expressed genes between HEK and HeLa
res <- results(dds, contrast = c("sample", "HEK", "HeLa"))

# Filter genes by fdr < 0.05
res.table = data.frame(res) %>%
	mutate(id = row.names(res)) %>%
	arrange(padj) %>%
	filter(padj < 0.05)

# Assign gene name to ID
tr = read.table("table/transcripts.txt", header = T, stringsAsFactors = F)
res.table = tr %>%
	inner_join(res.table)

# Further filter for highly expressed ones and sort by fold change
res.filter = res.table %>%
	filter(baseMean > 100) %>%
	arrange(log2FoldChange)

write.table(res.filter, "PROseq_difference.txt", quote = F, row.names = F, sep = "\t")
