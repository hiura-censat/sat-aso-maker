root<-'step4/all72_landscape';out<-file.path(root,'plots');dir.create(out,showWarnings=FALSE)
q<-read.delim('step4/mismatch/queries.tsv',na.strings='');ct<-read.delim(file.path(root,'chm13_censat_by_chromosome.tsv'),na.strings='');bins<-read.delim(file.path(root,'chm13_100kb_bins.tsv'),na.strings='')
fine<-read.delim('../VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed',header=FALSE,na.strings='');broad<-read.delim('../VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed',header=FALSE,na.strings='')
chrs<-c(paste0('chr',1:22),'chrX','chrY');colors<-c(CenSat='#CACFD4',alphaSat='#D46B74',HSat2='#247BA0',HSat3='#D99036',inside='#232C33',outside='#C53030')
for(kid in q$kmer_id){
 chrdata<-ct[ct$kmer_id==kid,];fam<-q$target_family[q$kmer_id==kid];png(file.path(out,paste0(kid,'_dominant_zoom.png')),width=2600,height=1500,res=180)
 if(sum(chrdata$total_exact_hits)==0){par(mar=c(5,5,5,2));plot.new();text(.5,.58,sprintf('%s | no exact hit on CHM13',kid),cex=1.8,font=2);text(.5,.42,'Haplotype group and chromosome count remain available in the landscape PNG.',cex=1.1);dev.off();next}
 dominant<-chrdata$chromosome[which.max(chrdata$total_exact_hits)];bc<-broad[broad$V1==dominant,];b<-bins[bins$kmer_id==kid & bins$chromosome==dominant,];inside<-b[b$outside_broad_CenSat==0,];outside<-b[b$outside_broad_CenSat==1,];pos<-b$bin_start_0based/1e6
 left<-max(0,min(c(bc$V2/1e6,pos))-2);right<-max(c(bc$V3/1e6,pos))+2;left<-floor(left);right<-ceiling(right)
 layout(matrix(1:2,ncol=1),heights=c(1.4,.65));par(oma=c(3,0,3,0));par(mar=c(4,5,4,2))
 yy<-numeric(max(1,ceiling((right-left)*10)+1));x<-left+(seq_along(yy)-1)/10
 if(nrow(inside)){ix<-inside$bin_100kb-as.integer(left*10)+1;yy[ix]<-inside$exact_hits}
 plot(x,yy,type='l',lwd=2,col=colors['inside'],xlim=c(left,right),ylim=c(0,max(yy,b$exact_hits,na.rm=TRUE)*1.05+1),xlab='',ylab='Exact hits / 100 kb',main=sprintf('%s | CHM13 %s: %s exact hits, %s outside broad CenSat',kid,dominant,format(sum(chrdata$total_exact_hits[chrdata$chromosome==dominant]),big.mark=','),format(sum(chrdata$outside_broad_CenSat[chrdata$chromosome==dominant]),big.mark=',')));grid(col='gray92')
 if(nrow(outside))points((outside$bin_100kb+.5)/10,outside$exact_hits,pch=16,col=colors['outside'],cex=.9)
 legend('topright',c('Inside broad CenSat','Outside broad CenSat'),col=c(colors['inside'],colors['outside']),lty=c(1,NA),pch=c(NA,16),bty='n',cex=.8)
 par(mar=c(5,5,3,2));plot(NA,xlim=c(left,right),ylim=c(0,3),yaxt='n',xlab=paste('CHM13',dominant,'position (Mb)'),ylab='',main='Annotation tracks')
 if(nrow(bc))for(j in seq_len(nrow(bc)))rect(bc$V2[j]/1e6,.1,bc$V3[j]/1e6,.5,col=colors['CenSat'],border=NA)
 for(family in c('alphaSat','HSat2','HSat3')){ff<-fine[fine$V1==dominant & fine$V4==family,];lev<-switch(family,alphaSat=.8,HSat2=1.5,HSat3=2.2);if(nrow(ff))for(j in seq_len(nrow(ff)))rect(ff$V2[j]/1e6,lev,ff$V3[j]/1e6,lev+.4,col=colors[family],border=NA)}
 axis(2,at=c(.3,1,1.7,2.4),labels=c('Broad CenSat','alphaSat','HSat2','HSat3'),las=2,cex.axis=.8);box()
 mtext(sprintf('%s | %s target | dominant CHM13 chromosome zoom',kid,fam),outer=TRUE,side=3,line=1,font=2,cex=1.25)
 mtext('The full-chromosome landscape PNG shows other chromosomes and all outside positions. Exact 16-mer hits only.',outer=TRUE,side=1,line=1,cex=.8);dev.off();cat('zoom',kid,dominant,'\n')
}
