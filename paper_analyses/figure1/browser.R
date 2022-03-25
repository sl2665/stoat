read.genelist=function(filename)
{
    table=read.delim(filename,header=F)
    if(ncol(table)==12)
		colnames(table)=c('chr','start','end','name','score','strand','cdsStart','cdsEnd','rgb','exonCount','exonStart','exonEnd')
	else
	{
		table=table[,1:6]
		colnames(table)=c('chr','start','end','name','score','strand')
	}
    return(table)
}

find.gene=function(genelist,name,isoform=1)
{
    gene=list()
    id=which(genelist[,4]==name)[isoform]
    if(is.na(id)) return(NULL)
	gene$name=genelist[id,4]
    gene$chr=genelist[id,1]
    gene$start=genelist[id,2]
    gene$end=genelist[id,3]
    gene$strand=genelist[id,6]
    if(ncol(genelist)==12)
	{
		gene$cdsStart=genelist[id,7]
		gene$cdsEnd=genelist[id,8]
		gene$exonCount=genelist[id,10]

		gene$exonStart=as.numeric(strsplit(toString(genelist[id,12]),",")[[1]])+gene$start
		gene$exonEnd=gene$exonStart+as.numeric(strsplit(toString(genelist[id,11]),",")[[1]])
	}
	else
	{
		gene$cdsStart=gene$end+1
		gene$cdsEnd=gene$start-1
		gene$exonCount=1
		gene$exonStart=c(gene$start)
		gene$exonEnd=c(gene$end)
	}
    return(gene)
}

genes.inrange=function(genelist,chr,start,end)
{
    inRange=(genelist[,2]>start & genelist[,2]<end)|(genelist[,3]>start & genelist[,3]<end)|(genelist[,2]<start & genelist[,3]>end)
    return(unique(genelist[genelist[,1]==chr & inRange,4]))
}

get.reads=function(bedgraph,chr,start,end,bin=1)
{
    system(paste('tabix ',bedgraph,".gz ",chr,":",formatC(start,format='d'),"-",formatC(end,format='d'),' > temp.tmp',sep=""))
#    system(paste('~/Sandbox/procap/code/bsearchbg',bedgraph,chr,formatC(start,format='d'),formatC(end,format='d'),formatC(bin,format='d'),'> temp.tmp'))
    n=ceiling((end-start)/bin)+1
    reads=numeric(n)
    
	data=try(read.delim('temp.tmp',header=F),silent=T)
	#system('rm temp.tmp')
    nr=nrow(data)
	if(is.null(nr)) return(reads[1:n])
	if(bin==1) for(i in 1:nr) {range=(data[i,2]+1):data[i,3]-start;if(all(range>0)) reads[range]=data[i,4] }
    else for(i in 1:nr) { range=floor((data[i,2]+1):data[i,3]-start)/bin+1; if(all(range>0)) reads[range]=reads[range]+data[i,4] }
    return(reads[1:n])
}

plot.reads=function(reads,start,y=0,height=1,bin=1,ymax=NULL,col='#000000')
{
    data=reads
    n=length(data)
    if(is.null(ymax)) ymax=1.05*max(data)
    if(ymax==0) ymax=-1.05*min(data)
    data=data/ymax
    data[data > 1]=1
    data[data < -1]=-1
    rect((1:n-1)*bin+start,y,(1:n)*bin+start,y+data*height,col=col,border=NA)
    return(ymax)
}

plot.gene=function(gene,y=0,height=1,col='#000000')
{
    lines(c(gene$start,gene$end),rep(y,2),col=col,lend=1)
    rect(gene$exonStart,rep(y-height*0.3,gene$exonCount),gene$exonEnd,rep(y+height*0.3,gene$exonCount),col=col,border=NA)
    if(gene$cdsStart<gene$cdsEnd)
    {
        cdsExonStart=c()
        cdsExonEnd=c()
        firstCodingExon=which(gene$exonStart>gene$cdsStart)[1]-1
        lastCodingExon=which(gene$exonEnd>=gene$cdsEnd)[1]
        if(is.finite(firstCodingExon)&is.finite(lastCodingExon))
        {
            if(firstCodingExon==lastCodingExon)
            {
                cdsExonStart=gene$cdsStart
                cdsExonEnd=gene$cdsEnd
            }
            else
            {
                cdsExonStart=c(gene$cdsStart,gene$exonStart[(firstCodingExon+1):lastCodingExon])
                cdsExonEnd=c(gene$exonEnd[firstCodingExon:(lastCodingExon-1)],gene$cdsEnd)
            }
            rect(cdsExonStart,rep(y-height*0.5,gene$exonCount),cdsExonEnd,rep(y+height*0.5,gene$exonCount),col=col,border=NA)
        }
    }
}

plot.browser=function(reads,genes,start,end,bin=1,col=NULL,gene.col=NULL,y=NULL,height=NULL,ymax=NULL,vline=NULL,line.gene=1)
{
    # Setup parameters
    par(ps=10,cex=1,cex.main=1,mar=c(1,4,2,1),mgp=c(1.5,0.2,0))
    n=ncol(reads)
    if(is.null(height)) height=rep(1,n)
    if(is.null(y))
    {
        y=c()
        for(i in 1:(n-1)) y[i]=sum(height[(i-1):n])
        y[n]=0
    }
    if(is.null(col)) col=rep('#000000',n)
    if(is.null(gene.col)) gene.col=rep('#000000',n)

    # Plot area
    plot(0,0,type='n',xlim=c(start,end),ylim=c(0,y[1]+height[1]+0.1+0.4*line.gene),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i',bty='n')
	if(!is.null(vline)) for(i in 1:length(vline)) lines(rep(vline[i],2),c(0,y[1]+height[1]),col='#c0c0c0')
    max=c()
    if(is.null(ymax)) for(i in 1:n) max[i]=plot.reads(reads[,i],start,y=y[i],bin=bin,col=col[i],height=height[i])
	else for(i in 1:n) max[i]=plot.reads(reads[,i],start,y=y[i],bin=bin,col=col[i],ymax=ymax[i],height=height[i])
	# Plot gene
	gene.Disp.mat = matrix(T, ncol=100, nrow=line.gene+1)

	if(length(genes)>0) for(i in 1:length(genes))
    {
		# Find empty space for displaying gene
		gene.rx.start = max(floor((genes[[i]]$start-start)/(end-start)*100),0) + 1
		gene.rx.end = min(floor((genes[[i]]$end-start)/(end-start)*100),99) + 1
		gene.Disp.line = 0
		for(j in 1:line.gene)
			if(all(gene.Disp.mat[j,gene.rx.start:gene.rx.end])) {
				gene.Disp.line=j
				break
			}
		if(gene.Disp.line == 0) next	# No empty space -> skip displaying the gene
		gene.Disp.mat[gene.Disp.line,gene.rx.start:gene.rx.end] = F
		geneypos=y[1]+height[1]+0.4*(line.gene-gene.Disp.line)+0.25
        
		plot.gene(genes[[i]],y=geneypos,height=0.3,col=gene.col[i])
        name=paste(genes[[i]]$name," ")
        if(length(genes)<10)
        {
            if(genes[[i]]$strand=='+')
            {
                text(genes[[i]]$start,geneypos,labels=name,cex=0.8,adj=1)
            }
            else
            {
                text(genes[[i]]$end,geneypos,labels=name,cex=0.8,adj=0)
            }
        }
    }
    return(max)
}

browser.param=list()
browser.setup=function(genelist=NULL,bedgraph=NULL,description=NULL,heights=NULL,col=NULL,hlines=NULL,label.y=NULL,readcount=NULL,gene.line=1)
{
    if(!is.null(genelist)) browser.param$genelist<<-read.genelist(genelist)
    if(!is.null(bedgraph))
    {
		browser.param$bedgraph<<-bedgraph
        if(is.null(heights)) browser.param$heights<<-rep(1,length(bedgraph))
        if(is.null(description)) browser.param$description<<-rep('',length(bedgraph))
    }
    if(!is.null(description)) browser.param$description<<-description
    if(!is.null(label.y)) browser.param$label.y<<-label.y
    if(!is.null(heights)) browser.param$heights<<-heights
    if(!is.null(col)) browser.param$color<<-col
    if(!is.null(hlines)) browser.param$hlines<<-hlines
	if(!is.null(readcount)) browser.param$readcount<<-readcount
    
    if(is.null(browser.param$plus.strand.gene.color)) browser.param$plus.strand.gene.color<<-'#b2182b'
    if(is.null(browser.param$minus.strand.gene.color)) browser.param$minus.strand.gene.color<<-'#2166ac'
   	if(is.null(browser.param$gene.line)) browser.param$gene.line<<-gene.line 
}

browser.setpos=function(chr,start,end)
{
    browser.param$chr<<-chr
    browser.param$start<<-start
    browser.param$end<<-end
    genes=genes.inrange(browser.param$genelist,chr,start,end)
    browser.param$dispgene<<-list()
    if(length(genes)>0)
    {
        if(!is.null(genes)) for(i in 1:length(genes))
			for(j in 1:browser.param$gene.line) {
				found.gene = find.gene(browser.param$genelist,genes[i],isoform=j)
				if(!is.null(found.gene)) browser.param$dispgene[[length(browser.param$dispgene)+1]]<<-found.gene
		}
        r=browser.param$plus.strand.gene.color
        b=browser.param$minus.strand.gene.color
        browser.param$genecolor<<-c()
        for(i in 1:length(browser.param$dispgene))
            if(browser.param$dispgene[[i]]$strand=='+') browser.param$genecolor<<-c(browser.param$genecolor,r)
            else browser.param$genecolor<<-c(browser.param$genecolor,b)
    }
}

browser.setgene=function(genename)
{
    g=find.gene(browser.param$genelist,genename)
    if(!is.na(g$name))
    {
        length=g$end-g$start
        browser.setpos(g$chr,floor(g$start-length/10),ceiling(g$end+length/10))
    }
}

browser.read=function(nbin=400,bin=NULL)
{
    browser.param$reads<<-c()
    if(is.null(bin)) browser.param$bin<<-ceiling((browser.param$end-browser.param$start)/nbin)
    else browser.param$bin<<-bin
    for(i in 1:length(browser.param$bedgraph)) browser.param$reads<<-cbind(browser.param$reads, get.reads(browser.param$bedgraph[i],browser.param$chr,browser.param$start,browser.param$end,browser.param$bin))
}

browser.draw=function(ymax=NULL,nbin=400,bin=NULL,read=T,vline=NULL)
{
	if(read)
	{
		if(!is.null(bin)) browser.read(bin=bin)
		else browser.read(nbin=nbin)
	}
    ncol=ncol(browser.param$reads)
    max=apply(browser.param$reads,2,max)
    min=apply(browser.param$reads,2,min)
    max[max==min]=1
    y=rep(0,ncol)-browser.param$heights*min/(max-min)
    if(ncol > 1 ) for(i in 1:(ncol-1)) y[i]=y[i]+sum(browser.param$heights[(i+1):ncol])
    if(is.null(browser.param$color)) col=rep('#000000',ncol)
    else col=browser.param$color
    max = plot.browser(browser.param$reads, browser.param$dispgene, browser.param$start, 
		       browser.param$end, bin=browser.param$bin, col=col,
		       gene.col = browser.param$genecolor, y=y, ymax=ymax,
		       height=browser.param$heights,vline=vline,line.gene=browser.param$gene.line)
    axis(3,lwd=0,lwd.tick=0.5,tck=-0.015,cex=0.4,at=axTicks(1),labels=formatC(axTicks(1),format='d',big.mark=','))
    text(rep(1.015*browser.param$start-0.015*browser.param$end,length(browser.param$description)), browser.param$label.y, labels=browser.param$description,xpd=T,adj=1)
    mtext(browser.param$chr,line=-1.8,adj=0,outer=T)
    abline(h=browser.param$hlines,lwd=0.5)
    box(lwd=0.5)
    return(max)
}

browser.setstranded=function(bedgraph,genelist=NULL,description=NULL,readcount=NULL)
{
    if(!is.null(genelist)) browser.setup(genelist=genelist)
    bgnames=c()
    col=c()
    n=length(bedgraph)    
    hlines=(1:n)*2
    label.y=(n:1)*2-1
	if(is.null(description)) description=bedgraph
	if(!is.null(readcount)) browser.setup(readcount=rbind(readcount,readcount))
    for(i in 1:n)
    {
        bgnames[2*i-1]=paste(bedgraph[i],'.pl.bedgraph',sep='')
        bgnames[2*i]=paste(bedgraph[i],'.mn.bedgraph',sep='')
        col=c(col,'#b2182b','#2166ac')
    }
    browser.setup(bedgraph=bgnames,hlines=hlines)
    browser.setup(description=description,label.y=label.y,col=col)
}

browser.print=function(filename='test.pdf',width=8,height=4,nbin=400,bin=NULL,ymax=NULL,normCount=NULL,read=T)
{
    pdf(filename,width=width,height=height)
	if(!is.null(normCount)) ymax=browser.param$readcount/1000000*normCount
    max=browser.draw(nbin=nbin,bin=bin,ymax=ymax,read=read)
    dev.off()
	return(max)
}

