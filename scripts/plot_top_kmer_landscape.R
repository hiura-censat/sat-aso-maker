root<-'step4/top_kmer_landscape'; out<-file.path(root,'plots');dir.create(out,showWarnings=FALSE)
ids<-c('k16_3e476af8','k16_a0d0562c','k16_4118ba0e');chr<-c(paste0('chr',1:22),'chrX','chrY');grp<-c(paste0('G',1:7),'NA')
b<-read.delim(file.path(root,'chm13_1Mb_bins.tsv'),na.strings='');h<-read.delim(file.path(root,'hap_total_exact_counts.tsv'),na.strings='');g<-read.delim(file.path(root,'group_chromosome_summary.tsv'),na.strings='');s<-read.delim(file.path(root,'chm13_chromosome_summary.tsv'),na.strings='')
rank<-read.delim('step4/ranked/mismatch_evaluated_ranked.tsv');names(rank)<-make.names(names(rank))
regioncolors<-c(HSat3='#D47A33',HSat2='#247BA0',HSat23_ambiguous='#9759A6',background_or_boundary='#8B9098')
for(kid in ids){
 seq<-rank$canonical_target_5to3[match(kid,rank$kmer_id)];d<-b[b$kmer_id==kid,];hh<-h[h$kmer_id==kid & h$hap!='chm13v2.0',];gg<-g[g$kmer_id==kid,];ss<-s[s$kmer_id==kid,];maxMb<-max(ss$chromosome_length_bp)/1e6
 png(file.path(out,paste0(kid,'.png')),width=3400,height=2650,res=180)
 layout(matrix(1:3,ncol=1),heights=c(1.4,.95,1.5));par(oma=c(3,2,4,2))
 par(mar=c(5,7,3,3));plot(NA,xlim=c(0,maxMb),ylim=c(.5,24.5),yaxt='n',xlab='CHM13 position (Mb; exact 16-mer hit)',ylab='',main=sprintf('CHM13 | %s total hits, %s in HSat3',format(sum(ss$total_hits),big.mark=','),format(sum(ss$HSat3_hits),big.mark=',')))
 axis(2,at=24:1,labels=chr,las=2,cex.axis=.7);abline(h=seq(.5,24.5,1),col='gray92');grid(nx=10,ny=NA,col='gray92')
 if(nrow(d))for(r in unique(d$region)){
  x<-d[d$region==r,];points(x$bin_1Mb+.5,25-match(x$chromosome,chr),pch=15,cex=pmin(1.8,.35+.24*log1p(x$exact_hits)),col=regioncolors[r])
 }
 legend('topright',names(regioncolors),pch=15,col=regioncolors,pt.cex=1.2,bty='n',cex=.75)
 par(mar=c(5,6,3,2));hh$group<-factor(ifelse(hh$HSat3_group=='NA','NA: core incomplete',hh$HSat3_group),levels=c(paste0('G',1:7),'NA: core incomplete'))
 boxes<-split(log1p(hh$HSat3_exact_hits),hh$group);boxes<-boxes[vapply(boxes,length,integer(1))>0]
 n<-vapply(boxes,length,integer(1));boxplot(boxes,names=paste0(sub('NA: core incomplete','NA',names(boxes)),' (n=',n,')'),col=c(rep('#E6A85F',7),'gray80')[match(names(boxes),levels(hh$group))],outline=FALSE,las=1,ylab='log1p(HSat3 exact hits per hap)',main='HSat3 common-core regional groups | per-hap HSat3 exact counts')
 stripchart(boxes,vertical=TRUE,method='jitter',pch=16,cex=.2,col=adjustcolor('#333333',alpha.f=.25),add=TRUE)
 par(mar=c(6,7,3,3));z<-matrix(NA_real_,nrow=length(chr),ncol=length(grp),dimnames=list(chr,grp))
 for(i in seq_len(nrow(gg))){x<-gg[i,];if(!is.na(x$median_HSat3_hits))z[x$chromosome,x$HSat3_group]<-log10(1+x$median_HSat3_hits)}
 zz<-z[24:1,,drop=FALSE];image(seq_along(grp),seq_along(chr),t(zz),col=colorRampPalette(c('#f6f4ee','#E8BA73','#C45F31','#783827'))(100),zlim=c(0,max(z,na.rm=TRUE)),axes=FALSE,xlab='HSat3 group (regional)',ylab='',main='Group × chromosome | median HSat3 exact-hit count per measurable hap')
 axis(1,at=seq_along(grp),labels=grp);axis(2,at=seq_along(chr),labels=rev(chr),las=2,cex.axis=.7);box()
 for(x in seq_along(grp))for(y in seq_along(chr))if(is.finite(zz[y,x]))text(x,y,sprintf('%.0f',10^zz[y,x]-1),cex=.48,col=if(zz[y,x]>.65*max(z,na.rm=TRUE))'white' else '#333333')
 mtext('Cells with no measurable HSat3 region are blank. Numbers are median exact-hit counts; colors use log10(1 + count).',side=1,line=4.5,cex=.75)
 mtext(paste('STEP 4 follow-up |',kid,'| canonical target',seq),outer=TRUE,side=3,line=1,font=2,cex=1.3)
 mtext('CHM13 plot includes both DNA orientations; group and chromosome counts use HSat3 annotated intervals and exact matches.',outer=TRUE,side=1,line=1,cex=.78)
 dev.off()
}

hits<-read.delim(gzfile(file.path(root,'chm13_all_exact_hits.tsv.gz')),na.strings='');hits<-hits[hits$chromosome=='chr9',];hits$bin100kb<-floor(hits$start_0based/1e5)
bed<-read.delim('step0/regions/chm13v2.0/HSat3.bed',header=FALSE);bed<-bed[bed$V1=='chr9',]
png(file.path(out,'chr9_100kb_zoom_all_three.png'),width=2800,height=1900,res=180)
par(mfrow=c(3,1),mar=c(4,5,3,2),oma=c(4,0,3,0));cols<-c('#C3532C','#247BA0','#7060A5')
for(i in seq_along(ids)){
 d<-hits[hits$kmer_id==ids[i],];tab<-table(d$bin100kb);x<-seq(480,770);y<-numeric(length(x));y[match(as.integer(names(tab)),x,nomatch=0)[match(as.integer(names(tab)),x,nomatch=0)>0]]<-as.numeric(tab)[match(as.integer(names(tab)),x,nomatch=0)>0]
 plot(x/10,y,type='n',xlim=c(48,77),ylim=c(0,max(y)*1.05+1),xlab=if(i==3)'CHM13 chr9 position (Mb)' else '',ylab='Exact hits / 100 kb',main=ids[i])
 for(j in seq_len(nrow(bed)))rect(bed$V2[j]/1e6,0,bed$V3[j]/1e6,par('usr')[4],col=adjustcolor('#E6A85F',alpha.f=.17),border=NA)
 lines(x/10,y,col=cols[i],lwd=2);grid(col='gray90');legend('topright',c('Exact hits (both orientations)','HSat3 annotation'),col=c(cols[i],'#E6A85F'),lty=c(1,NA),lwd=c(2,NA),pch=c(NA,15),bty='n',cex=.75)
}
mtext('CHM13 chr9 | 100 kb bins for the three top pan-HSat3 16-mers',outer=TRUE,side=3,line=1,font=2,cex=1.2)
mtext('Shaded areas are CHM13 HSat3 BED intervals; curves count exact matches on both DNA orientations.',outer=TRUE,side=1,line=1,cex=.8);dev.off()
