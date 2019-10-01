library(ggplot2)
library(dplyr)
library(tidyr)

# Read figure arguments
arg = commandArgs(trailingOnly=T)
filename = arg[1]
width = as.numeric(arg[2])
height = as.numeric(arg[3])

# Read palmat file
spos = 200
mat = read.table("_tmp/palmat",header=F,fill=T,flush=T)
colnames(mat) = c("id",(-spos):(500-spos))

# Convert to long format
mat.long = mat %>%
  gather(pos,count,-id) %>%
  mutate(pos=as.numeric(pos),count=as.numeric(count)) %>%
  separate(id, c("gene", "group"), sep = ";")

mat.longres = mat.long %>% 
  mutate(pos = round(pos/5)*5) %>%

# Open pdf file
pdf(filename, width = width, height = height)
par(mar = c(2,2,1,2), mgp = c(1.5, 0.5, 0))

# If no group has been assigned, plot individual genes
if(length(unique(mat.long$group)) == 1) {
  pal = mat.longres %>%
    group_by(pos, gene) %>%
    summarise(readCount = sum(count))
  med = mat.long %>%
  	filter(pos >= 0 & pos <=250) %>%
    group_by(gene) %>%
	summarise(mLen = median(rep(pos, count)))
  ggplot(pal, aes(x = pos, y = readCount)) +
	geom_vline(data = med, aes(xintercept = mLen), col = "grey") +
    geom_step(aes(col = gene), stat = "identity") +
    facet_grid(gene~., scales = "free_y") +
    xlim(0,250) +
    xlab(expression(Poly*(A)*length[TED-seq])) +
    ylab("TED-seq UMI counts") +
    theme_bw() +
	scale_color_brewer(palette="Dark2")
} else {
# If more than one group has been assigned, plot sum of all reads in the group
  pal = mat.longres %>%
    group_by(pos, group) %>%
    summarise(readCount = sum(count))
  med = mat.long %>%
  	filter(pos >= 0 & pos <=250) %>%
    group_by(group) %>%
	summarise(mLen = median(rep(pos, count)))
  ggplot(pal, aes(x = pos, y = readCount)) +
	geom_vline(data = med, aes(xintercept = mLen), col = "grey") +
    geom_step(aes(col = group), stat = "identity") +
    facet_grid(group~., scales = "free_y") +
    xlim(0,250) +
    xlab(expression(Poly*(A)*length[TED-seq])) +
    ylab("TED-seq UMI counts") +
    theme_bw() +
	scale_color_brewer(palette="Dark2")
}

dev.off()
