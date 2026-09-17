#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else 'step0'
out <- file.path(root,'plots'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
s <- read.delim(file.path(root,'region_stats.tsv'),check.names=FALSE)
m <- read.delim(file.path(root,'manifest.tsv'),check.names=FALSE)
cstat <- read.delim(file.path(root,'chromosome_status.tsv'),check.names=FALSE)
chroms <- paste0('chr',c(1:22,'X','Y'))
families <- c('HSat2','HSat3','background')
colors <- c(HSat2='#267EAE',HSat3='#E59432',background='#718096')
ref <- m$hap[m$sample=='CHM13']; eligible <- m$hap[m$status=='PASS']
haps <- c(sort(setdiff(eligible,ref)),intersect(ref,eligible))
write.table(data.frame(row_from_top=seq_along(haps),hap=haps),file.path(out,'heatmap_row_order.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
start <- function(name,w=2200,h=1650){png(file.path(out,name),width=w,height=h,res=180);par(family='sans',col.axis='#354052',col.lab='#354052',col.main='#172B4D',fg='#354052')}
finish <- function(){dev.off()}
agg <- aggregate(cbind(bp,ACGT_bp,non_ACGT_bp,valid_16mer_starts)~region,s,sum)
a <- agg[match(families,agg$region),]
start('01_step0_overview.png')
par(mfrow=c(2,2),mar=c(6,5,4,2),oma=c(4,0,4,0))
counts <- c(sum(m$status=='PASS'),sum(grepl('incomplete',m$status)),sum(grepl('invalid_annotation',m$status)))
x <- barplot(counts,names.arg=c('Eligible','Incomplete / failed','Invalid coordinates'),col=c('#268878','#C59130','#C75857'),ylim=c(0,max(counts)*1.15),ylab='Input sequence sets',main='Input selection',cex.names=.85)
text(x,counts,labels=counts,pos=3,font=2)
vals <- a$bp/1e9
x <- barplot(vals,names.arg=families,col=colors[families],log='y',ylim=c(1,max(vals)*2),ylab='Total sequence length (Gb; log scale)',main='Annotated target and background')
text(x,vals,labels=sprintf('%.2f Gb',vals),pos=3,cex=.9)
vals <- a$valid_16mer_starts/1e9
x <- barplot(vals,names.arg=families,col=colors[families],log='y',ylim=c(1,max(vals)*2),ylab='Valid 16-mer starts (billions; log scale)',main='Usable sequence windows')
text(x,vals,labels=sprintf('%.3f B',vals),pos=3,cex=.9)
vals <- 100*a$non_ACGT_bp/a$bp
x <- barplot(vals,names.arg=families,col=colors[families],ylim=c(0,max(vals)*1.3+1e-6),ylab='Non-ACGT bases (%)',main='Sequence ambiguity')
text(x,vals,labels=sprintf('%.5f%%',vals),pos=3,cex=.9)
mtext('STEP 0 | HSat2/3 input preparation',outer=TRUE,side=3,line=1.5,cex=1.5,font=2)
mtext(sprintf('%d eligible sets (including CHM13) | %s set-chromosome pairs',length(eligible),format(nrow(s)/5,big.mark=',')),outer=TRUE,side=3,line=.1,cex=.95)
mtext('Background = T2T-passed sequence outside the existing HSat2/3 BED; it may contain other satellites.',outer=TRUE,side=1,line=1,cex=.85)
mtext('Valid starts are sequence opportunities, not observed k-mer counts. Region boundaries and non-ACGT bases break windows.',outer=TRUE,side=1,line=2.3,cex=.8)
finish()

start('02_hsat23_by_chromosome.png',2500,1900)
par(mfrow=c(3,1),mar=c(3.5,5,3,2),oma=c(4,0,4,0))
for(fam in c('HSat2','HSat3')){
 d <- s[s$region==fam & !(s$hap %in% ref),]
 vals <- lapply(chroms,function(c) d$bp[d$chromosome==c]/1e6)
 boxplot(vals,names=sub('chr','',chroms),col=adjustcolor(colors[fam],alpha.f=.6),outline=TRUE,pch=16,cex=.35,ylab=paste(fam,'length (Mb)'),main=paste(fam,'per eligible haplotype and chromosome'),las=1)
 r <- s[s$region==fam & s$hap %in% ref,]
 points(match(r$chromosome,chroms),r$bp/1e6,pch=18,col='#922D50',cex=1)
 legend('topright',legend=c('Haplotype distribution (zeros included)','CHM13 reference'),fill=c(adjustcolor(colors[fam],alpha.f=.6),NA),pch=c(NA,18),col=c(colors[fam],'#922D50'),bty='n',cex=.85)
}
n <- s[s$region=='HSat23' & !(s$hap %in% ref),]
counts <- vapply(chroms,function(c) sum(n$chromosome==c),integer(1))
x <- barplot(counts,names.arg=sub('chr','',chroms),col='#78909C',ylim=c(0,max(counts)*1.2),ylab='Eligible haplotypes',xlab='',main='Denominator: eligible T2T-passed haplotypes')
text(x,counts,labels=counts,pos=3,cex=.75)
mtext('HSat2/3 annotation length by chromosome',outer=TRUE,side=3,line=1.5,cex=1.5,font=2)
mtext('Unavailable chromosomes are excluded, not counted as zero. Boxes: median and IQR; whiskers: 1.5 x IQR.',outer=TRUE,side=1,line=1,cex=.85)
mtext('CHM13 is shown separately and excluded from the haplotype distributions and denominators.',outer=TRUE,side=1,line=2.3,cex=.85)
finish()

start('03_haplotype_chromosome_heatmaps.png',2700,2800)
layout(matrix(1:6,nrow=2,byrow=TRUE),heights=c(12,1.2))
par(oma=c(3,2,5,1))
matrices <- list()
for(fam in c('HSat2','HSat3')){
 z <- matrix(NA_real_,length(haps),length(chroms),dimnames=list(haps,chroms))
 d <- s[s$region==fam,]
 z[cbind(match(d$hap,haps),match(d$chromosome,chroms))] <- d$bp/1000
 matrices[[fam]] <- z
}
coverage <- ifelse(is.na(matrices$HSat2),0,1)
heat <- function(z,palette,zlim,title){
 par(mar=c(4,4.5,3,1))
 plot.new();plot.window(xlim=c(.5,24.5),ylim=c(.5,length(haps)+.5),xaxs='i',yaxs='i')
 rect(.5,.5,24.5,length(haps)+.5,col='#DBDFE5',border=NA)
 image(seq_along(chroms),seq_along(haps),t(z[nrow(z):1,,drop=FALSE]),col=palette,zlim=zlim,add=TRUE)
 axis(1,at=1:24,labels=sub('chr','',chroms),las=2,cex.axis=.85)
 idx <- unique(c(1,seq(100,length(haps),100),length(haps)))
 axis(2,at=length(haps)+1-idx,labels=ifelse(idx==length(haps),'CHM13',idx),las=1,cex.axis=.8)
 title(main=title,xlab='Chromosome',ylab='Haplotype row (alphabetical; CHM13 last)')
 box()
}
heat(coverage,c('#DBDFE5','#268878'),c(0,1),'T2T-passed / eligible chromosome')
limit <- ceiling(log10(1+max(unlist(matrices),na.rm=TRUE)))
pal2 <- colorRampPalette(c('#FFFFFF','#CFEAF4','#267EAE','#103957'))(256)
pal3 <- colorRampPalette(c('#FFFFFF','#FCE5B8','#E59432','#71380E'))(256)
heat(log10(1+matrices$HSat2),pal2,c(0,limit),'HSat2 annotated length')
heat(log10(1+matrices$HSat3),pal3,c(0,limit),'HSat3 annotated length')
par(mar=c(2,3,1,1));plot.new();legend('center',legend=c('Eligible chromosome','Unavailable (not zero)'),fill=c('#268878','#DBDFE5'),bty='n',cex=.95)
for(pal in list(pal2,pal3)){
 par(mar=c(3,4.5,1,1))
 image(seq(0,limit,length.out=256),1,matrix(seq(0,limit,length.out=256),ncol=1),col=pal,axes=FALSE,xlab='',ylab='')
 axis(1,at=0:limit,labels=format(10^(0:limit)-1,scientific=FALSE,trim=TRUE),cex.axis=.85)
 mtext('Annotated length (kb; log10(1 + kb) colors)',side=1,line=2,cex=.8)
}
mtext('Haplotype x chromosome | annotation coverage and HSat2/3 length',outer=TRUE,side=3,line=2.5,cex=1.5,font=2)
mtext('574 eligible sequence sets. Rows are alphabetical, not clustered; row IDs are provided in heatmap_row_order.tsv.',outer=TRUE,side=3,line=1,cex=.95)
mtext('Gray = unavailable chromosome; white = eligible chromosome with zero annotated family length.',outer=TRUE,side=1,line=1,cex=.9)
finish()
writeLines(c('STEP 0 plots (PNG, 180 dpi)','',
'01_step0_overview.png: input selection, target/background length, valid 16-mer starts and non-ACGT fraction. Includes CHM13.',
'02_hsat23_by_chromosome.png: per-chromosome HSat2/3 length distributions and eligible-haplotype denominators. CHM13 is a separate diamond.',
'03_haplotype_chromosome_heatmaps.png: eligibility and family lengths, with the same row order and color scale for both families.',
'heatmap_row_order.tsv: exact row identities from top to bottom.',
'',
'No clustering or k-mer abundance analysis is performed here. Missing chromosomes are never treated as biological zeros.',
'Background means outside existing HSat2/3 annotations on eligible chromosomes, not confirmed non-CenSat sequence.',
'Reproduce from the project directory: Rscript scripts/plot_step0.R'),file.path(out,'README.txt'))
