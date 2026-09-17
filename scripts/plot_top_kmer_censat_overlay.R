root<-'step4/top_kmer_landscape';out<-file.path(root,'plots');dir.create(out,showWarnings=FALSE)
ids<-read.delim(file.path(root,'selected_kmers.tsv'))$kmer_id;chr<-c(paste0('chr',1:22),'chrX','chrY')
fine<-read.delim('../VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed',header=FALSE,na.strings='');broad<-read.delim('../VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed',header=FALSE,na.strings='')
hits<-read.delim(gzfile(file.path(root,'chm13_all_exact_hits_with_censat.tsv.gz')),na.strings='');tab<-read.delim(file.path(root,'chm13_censat_inside_outside_by_chromosome.tsv'),na.strings='');tot<-read.delim(file.path(root,'chm13_censat_inside_outside_totals.tsv'),na.strings='');cs<-read.delim(file.path(root,'chm13_chromosome_summary.tsv'),na.strings='')
col<-c(CenSat='#CACFD4',alphaSat='#D46B74',HSat2='#247BA0',HSat3='#D99036',inside='#232C33',outside='#C53030')
for(kid in ids){
 d<-hits[hits$kmer_id==kid,];d$bin100kb<-floor(d$start_0based/1e5);aggregate<-if(nrow(d))stats::aggregate(list(exact_hits=d$start_0based),list(chromosome=d$chromosome,bin100kb=d$bin100kb,outside=d$outside_broad_censat),length) else data.frame(chromosome=character(),bin100kb=integer(),outside=integer(),exact_hits=integer())
 t<-tab[tab$kmer_id==kid,];tt<-tot[tot$kmer_id==kid,];ss<-cs[cs$kmer_id==kid,]
 png(file.path(out,paste0(kid,'_censat_overlay.png')),width=3400,height=2700,res=180)
 par(mar=c(5,8,5,9),oma=c(3,0,3,0));plot(NA,xlim=c(0,285),ylim=c(.5,24.5),yaxt='n',xaxt='n',xlab='CHM13 position (Mb; 100 kb hit bins)',ylab='',main=sprintf('%s | CenSat inside %s; outside %s | fine HSat3 %s',kid,format(tt$inside_broad_CenSat,big.mark=','),format(tt$outside_broad_CenSat,big.mark=','),format(tt$fine_HSat3_hits,big.mark=',')))
 axis(1,at=seq(0,250,25));axis(2,at=24:1,labels=chr,las=2,cex.axis=.78);abline(v=seq(0,250,25),col='gray92');abline(h=seq(.5,24.5,1),col='gray92')
 for(ci in seq_along(chr)){
  ch<-chr[ci];y<-25-ci;len<-ss$chromosome_length_bp[ss$chromosome==ch]/1e6;segments(0,y-.37,len,y-.37,col='gray85',lwd=2)
  bc<-broad[broad$V1==ch,];if(nrow(bc))rect(bc$V2/1e6,y-.33,bc$V3/1e6,y-.21,col=col['CenSat'],border=NA)
  for(family in c('alphaSat','HSat2','HSat3')){
   ff<-fine[fine$V1==ch & fine$V4==family,];level<-switch(family,alphaSat=-.15,HSat2=.01,HSat3=.17)
   if(nrow(ff))rect(ff$V2/1e6,y+level,ff$V3/1e6,y+level+.12,col=col[family],border=NA)
  }
  hh<-aggregate[aggregate$chromosome==ch,];if(nrow(hh))points((hh$bin100kb+.5)/10,rep(y+.38,nrow(hh)),pch=16,cex=pmin(1.6,.34+.22*log1p(hh$exact_hits)),col=ifelse(hh$outside==1,col['outside'],col['inside']))
  r<-t[t$chromosome==ch,];text(255,y,sprintf('in %s | out %s',format(r$inside_broad_CenSat,big.mark=','),format(r$outside_broad_CenSat,big.mark=',')),adj=0,cex=.72,col=if(r$outside_broad_CenSat>0)col['outside'] else '#384148')
 }
 text(255,24.45,'CenSat hits by chromosome',adj=0,cex=.72,font=2)
 legend('topleft',legend=c('Broad CenSat boundary','Fine alphaSat','Fine HSat2','Fine HSat3','Exact hits inside CenSat','Exact hits outside CenSat'),fill=c(col['CenSat'],col['alphaSat'],col['HSat2'],col['HSat3'],NA,NA),pch=c(NA,NA,NA,NA,16,16),col=c(rep(NA,4),col['inside'],col['outside']),border=NA,bty='n',cex=.72)
 mtext('Each chromosome has separate tracks: broad CenSat, alphaSat, HSat2, HSat3, and exact-hit bins. Right column shows all exact hits inside/outside broad CenSat.',outer=TRUE,side=1,line=1,cex=.76)
 dev.off()
}

# Independent view of the chromosome with the most selected exact hits.
counts<-sort(table(hits$chromosome),decreasing=TRUE);focus<-if(length(counts))names(counts)[1] else 'chr9'
chosen<-hits[hits$chromosome==focus,];bc<-broad[broad$V1==focus,];ff<-fine[fine$V1==focus & fine$V4=='HSat3',]
x<-if(nrow(chosen))seq(max(0,floor(min(chosen$start_0based)/1e5)-10),floor(max(chosen$start_0based)/1e5)+10) else seq(0,50)
png(file.path(out,'dominant_censat_and_hits_separate_tracks.png'),width=2800,height=2100,res=180)
layout(matrix(1:4,ncol=1),heights=c(1,1,1,.5));par(oma=c(3,0,3,0));pal<-c('#C3532C','#247BA0','#7060A5')
for(i in seq_along(ids)){
 d<-chosen[chosen$kmer_id==ids[i],];bn<-floor(d$start_0based/1e5);tab100<-table(bn);y<-numeric(length(x));ix<-match(as.integer(names(tab100)),x);y[ix[!is.na(ix)]]<-as.numeric(tab100)[!is.na(ix)]
 par(mar=c(3,5,3,2));plot(x/10,y,type='l',ylim=c(0,max(y)+1),lwd=2,col=pal[i],xlab='',ylab='Hits / 100 kb',main=sprintf('%s | %s exact hits',ids[i],format(nrow(d),big.mark=',')));grid(col='gray92')
}
par(mar=c(5,5,2,2));plot(NA,xlim=range(x/10),ylim=c(0,2),yaxt='n',xlab=paste('CHM13',focus,'position (Mb)'),ylab='',main='Annotation intervals')
for(j in seq_len(nrow(bc)))rect(bc$V2[j]/1e6,.2,bc$V3[j]/1e6,.6,col=col['CenSat'],border=NA)
for(j in seq_len(nrow(ff)))rect(ff$V2[j]/1e6,1.2,ff$V3[j]/1e6,1.6,col=col['HSat3'],border=NA)
axis(2,at=c(.4,1.4),labels=c('Broad CenSat','Fine HSat3'),las=2,cex.axis=.8);box()
mtext(paste('CHM13',focus,'| exact hits and CenSat annotation on separate tracks'),outer=TRUE,side=3,line=1,font=2,cex=1.2)
dev.off()
