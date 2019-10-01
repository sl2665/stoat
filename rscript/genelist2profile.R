library(ggplot2)
library(dplyr)
library(tidyr)

# Read figure arguments
arg = commandArgs(trailingOnly=T)

sampleFilename = arg[1]
geneFilename = arg[2]
width = as.numeric(arg[3])
height = as.numeric(arg[4])
if(arg[5] == "true") {
  overlaySample = TRUE
} else {
  overlaySample = FALSE
}
if(arg[6] == "true") {
  overlayGene = TRUE
} else {
  overlayGene = FALSE
}
filename = arg[7]

# Read TEDseq dir names
TEDseq.info = read.csv(sampleFilename,
                            header = F,
                            stringsAsFactors = F)
TEDseq.dirname = TEDseq.info[,1]
TEDseq.description = TEDseq.info[,2]
TEDseq.count = length(TEDseq.dirname)

readPALmatrix = function(i) {
  matobj = readRDS(paste0(TEDseq.dirname[i], "/Rdata/PALdata.rds"))
  return(data.frame(matobj$long, sample = TEDseq.description[i]))
}

PAL.list = lapply(1:TEDseq.count, readPALmatrix)
PAL.long = bind_rows(PAL.list)
rm(PAL.list)

# Read gene list
geneAnnotation = read.table("_tmp/annot",
                            header = F, stringsAsFactors = F,
                            col.names = c("id", "name"))

genes.info = read.table(geneFilename, header = F, stringsAsFactors = F)

# find gene id first
if(all(genes.info[,1] %in% geneAnnotation$id)) {
  # All gene id
  colnames(genes.info)[1] = "id"
} else {
  # Need to find gene name
  colnames(genes.info)[1] = "name"
  genes.info = inner_join(genes.info, geneAnnotation, by = "name")
}

# Filter for the PAL matrix with the gene id
pal = PAL.long %>%
  inner_join(genes.info, by = "id")

pdf(filename, width = width, height = height)
par(mar = c(2,2,1,2), mgp = c(1.5, 0.5, 0))

# if genes names are provided, swap id to name
pal = pal %>% {
  if("name" %in% names(.)) mutate(., id = name) %>%
    select(-name)
  }
genes.info = genes.info %>% {
  if("name" %in% names(.)) select(., -name)
  }

if(ncol(genes.info)==1) {
  # Display individual genes
  med = pal %>%
    filter(pos >= 0 & pos <=250) %>%
    group_by(id, sample) %>%
    summarise(mLen = median(rep(pos, count))) %>%
    ungroup()
  pal2 = pal %>%
    mutate(pos = round(pos/5)*5) %>%
    group_by(pos, id, sample) %>%
    summarise(readCount = sum(count)) %>%
    ungroup()
  gg = ggplot(pal2, aes(x = pos, y = readCount)) +
    xlim(0,250) +
    xlab(expression(Poly*(A)*length[TED-seq])) +
    ylab("TED-seq UMI counts") +
    theme_bw() +
    scale_color_brewer(palette="Dark2")
  if(overlaySample) {
    print(gg + 
            geom_vline(data = med, aes(xintercept = mLen, col = sample), linetype = "dashed") +
            geom_step(aes(col = sample), stat = "identity") +
            facet_grid(id ~., scales = "free_y"))    
  } else if(overlayGene) {
    print(gg + 
            geom_vline(data = med, aes(xintercept = mLen, col = gene), linetype = "dashed") +
            geom_step(aes(col = id), stat = "identity") +
            facet_grid(sample~., scales = "free_y"))    
  } else {
    print(gg + 
            geom_vline(data = med, aes(xintercept = mLen), col = "grey") +
            geom_step(aes(), stat = "identity") + 
            facet_wrap(sample~id, scales = "free_y", nrow = TEDseq.count))    
  }    
} else {
  # Display by gene groups
  pal2 = pal %>%
    mutate(group = V2) %>%
    select(-V2)
  med = pal2 %>%
    filter(pos >= 0 & pos <=250) %>%
    group_by(group, sample) %>%
    summarise(mLen = median(rep(pos, count))) %>%
    ungroup()
  pal2 = pal2 %>%
    mutate(pos = round(pos/5)*5) %>%
    group_by(pos, group, sample) %>%
    summarise(readCount = sum(count)) %>%
    ungroup()
  gg = ggplot(pal2, aes(x = pos, y = readCount)) +
    xlim(0,250) +
    xlab(expression(Poly*(A)*length[TED-seq])) +
    ylab("TED-seq UMI counts") +
    theme_bw() +
    scale_color_brewer(palette="Dark2")
  
  if(overlaySample) {
    print(gg + 
            geom_vline(data = med, aes(xintercept = mLen, col = sample), linetype = "dashed") +
            geom_step(aes(col = sample), stat = "identity") +
            facet_grid(group ~., scales = "free_y"))
  } else if(overlayGene) {
    print(gg + 
            geom_vline(data = med, aes(xintercept = mLen, col = group), linetype = "dashed") +
            geom_step(aes(col = group), stat = "identity") +
            facet_grid(sample~., scales = "free_y"))
  } else {
    print(gg + 
            geom_vline(data = med, aes(xintercept = mLen), col = "grey") +
            geom_step(aes(), stat = "identity") + 
            facet_wrap(group~sample, scales = "free_y", ncol = TEDseq.count))
  }
}

dev.off()
