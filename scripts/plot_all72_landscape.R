root<-'step4/all72_landscape';out<-file.path(root,'plots');dir.create(out,showWarnings=FALSE)
args<-commandArgs(trailingOnly=TRUE);pilot<-if(length(args))as.integer(args[1]) else 72
q<-read.delim('step4/mismatch/queries.tsv',na.strings='');rank<-read.delim('step4/ranked/mismatch_evaluated_ranked.tsv',na.strings='');total<-read.delim(file.path(root,'chm13_censat_totals.tsv'),na.strings='');chrtable<-read.delim(file.path(root,'chm13_censat_by_chromosome.tsv'),na.strings='');bins<-read.delim(file.path(root,'chm13_100kb_bins.tsv'),na.strings='');hap<-read.delim(gzfile(file.path(root,'hap_total_exact_counts.tsv.gz')),na.strings='');gchr<-read.delim(file.path(root,'group_chromosome_summary.tsv'),na.strings='')
fine<-read.delim('../VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed',header=FALSE,na.strings='');broad<-read.delim('../VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed',header=FALSE,na.strings='')
chrs<-c(paste0('chr',1:22),'chrX','chrY');colors<-c(CenSat='#CACFD4',alphaSat='#D46B74',HSat2='#247BA0',HSat3='#D99036',inside='#232C33',outside='#C53030')
for(kid in head(q$kmer_id,pilot)){
 t<-total[total$kmer_id==kid,];ct<-chrtable[chrtable$kmer_id==kid,];b<-bins[bins$kmer_id==kid,];h<-hap[hap$kmer_id==kid & hap$hap!='chm13v2.0',];gc<-gchr[gchr$kmer_id==kid,];r<-rank[rank$kmer_id==kid,]
 fam<-as.character(t$target_family);grp<-paste0('G',1:8)
 png(file.path(out,paste0(kid,'_landscape.png')),width=3400,height=3350,res=180)
 layout(matrix(1:3,ncol=1),heights=c(1.7,.78,1.28));par(oma=c(3,0,4,0))
 par(mar=c(5,7,3,9));plot(NA,xlim=c(0,285),ylim=c(.5,24.5),yaxt='n',xaxt='n',xlab='CHM13 position (Mb; exact hits in 100 kb bins)',ylab='',main=sprintf('CHM13 | %s inside broad CenSat | %s outside | %s fine %s',format(t$inside_broad_CenSat,big.mark=','),format(t$outside_broad_CenSat,big.mark=','),format(t[[paste0('fine_',fam,'_hits')]],big.mark=','),fam))
 axis(1,at=seq(0,250,25));axis(2,at=24:1,labels=chrs,las=2,cex.axis=.7);abline(v=seq(0,250,25),col='gray93');abline(h=seq(.5,24.5,1),col='gray93')
 for(ci in seq_along(chrs)){
  ch<-chrs[ci];y<-25-ci;bc<-broad[broad$V1==ch,];if(nrow(bc))rect(bc$V2/1e6,y-.33,bc$V3/1e6,y-.21,col=colors['CenSat'],border=NA)
  for(family in c('alphaSat','HSat2','HSat3')){ff<-fine[fine$V1==ch & fine$V4==family,];lev<-switch(family,alphaSat=-.15,HSat2=.01,HSat3=.17);if(nrow(ff))rect(ff$V2/1e6,y+lev,ff$V3/1e6,y+lev+.12,col=colors[family],border=NA)}
  bb<-b[b$chromosome==ch,];if(nrow(bb))points((bb$bin_100kb+.5)/10,rep(y+.38,nrow(bb)),pch=16,cex=pmin(1.6,.3+.21*log1p(bb$exact_hits)),col=ifelse(bb$outside_broad_CenSat==1,colors['outside'],colors['inside']))
  z<-ct[ct$chromosome==ch,];text(255,y,sprintf('in %s | out %s',format(z$inside_broad_CenSat,big.mark=','),format(z$outside_broad_CenSat,big.mark=',')),adj=0,cex=.64,col=if(z$outside_broad_CenSat>0)colors['outside'] else '#3A444C')
 }
 legend('topleft',legend=c('Broad CenSat','Fine alphaSat','Fine HSat2','Fine HSat3','Exact inside','Exact outside'),fill=c(colors['CenSat'],colors['alphaSat'],colors['HSat2'],colors['HSat3'],NA,NA),pch=c(NA,NA,NA,NA,16,16),col=c(rep(NA,4),colors['inside'],colors['outside']),border=NA,bty='n',cex=.7)
 par(mar=c(5,6,3,2))
 sizes<-tapply(gc$measurable_haps,gc$chromosome,sum)
 if(length(sizes))barplot(sizes,las=2,col=if(fam=='HSat2')'#94C4D6' else '#E6A85F',ylab='Assigned measurable haps',main=sprintf('%s groups fitted separately by chromosome',fam)) else {plot.new();title('No supported chromosome groups')}
 par(mar=c(5,7,3,2));zz<-matrix(NA_real_,nrow=24,ncol=length(grp),dimnames=list(chrs,grp))
 for(i in seq_len(nrow(gc))){x<-gc[i,];if(x$group%in%grp&&!is.na(x$median_target_hits))zz[x$chromosome,x$group]<-log10(1+x$median_target_hits)}
 mat<-zz[24:1,,drop=FALSE];mx<-max(1,mat,na.rm=TRUE);image(seq_along(grp),seq_along(chrs),t(mat),col=colorRampPalette(c('#f6f4ee','#E8BA73','#C45F31','#783827'))(100),zlim=c(0,mx),axes=FALSE,xlab='Group (local to each chromosome)',ylab='',main=sprintf('Group × chromosome | median %s exact-hit count per measurable hap',fam))
 axis(1,at=seq_along(grp),labels=grp);axis(2,at=seq_along(chrs),labels=rev(chrs),las=2,cex.axis=.7);box()
 for(x in seq_along(grp))for(y in seq_along(chrs))if(is.finite(mat[y,x]))text(x,y,sprintf('%.0f',10^mat[y,x]-1),cex=.44,col=if(mat[y,x]>.65*mx)'white' else '#333333')
 mtext(sprintf('%s | %s | canonical target %s | %s | %s',kid,fam,r$canonical_target_5to3,r$candidate_type,r$category_robustness_tier),outer=TRUE,side=3,line=1,font=2,cex=1.25)
 mtext('Exact 16-mer hits only. Group counts are inside the target-family BED; missing chromosomes are excluded from medians.',outer=TRUE,side=1,line=1,cex=.76)
 dev.off();cat('plotted',kid,'\n')
}

if(pilot>=72){
 candidate<-read.delim(file.path(root,'chm13_candidate_censat_summary.tsv'),na.strings='');hit<-candidate[candidate$CHM13_outside_broad_CenSat>0,];hit<-hit[order(-hit$CHM13_outside_broad_CenSat),]
 png(file.path(out,'all72_censat_outside_overview.png'),width=2700,height=1700,res=180)
 par(mfrow=c(1,2),mar=c(5,5,4,1),oma=c(3,0,3,0))
 plot(log10(1+candidate$CHM13_total_exact_hits),candidate$CHM13_outside_broad_CenSat,pch=21,bg=ifelse(candidate$target_family=='HSat2',colors['HSat2'],colors['HSat3']),col='white',cex=1.1,xlab='log10(1 + CHM13 total exact hits)',ylab='Outside broad CenSat exact hits',main='72 k-mers | CenSat outside hits');grid(col='gray90');legend('topright',c('HSat2','HSat3'),pch=21,pt.bg=c(colors['HSat2'],colors['HSat3']),bty='n')
 if(nrow(hit))barplot(rev(hit$CHM13_outside_broad_CenSat),names.arg=rev(sub('k16_','',hit$kmer_id)),horiz=TRUE,las=1,col=ifelse(rev(hit$target_family)=='HSat2',colors['HSat2'],colors['HSat3']),xlab='Exact hits outside broad CenSat',main=sprintf('%d k-mers with outside hits',nrow(hit)),cex.names=.65) else plot.new()
 mtext('CHM13 exact 16-mer counts. Candidate list sorted by outside hit count is in chm13_candidate_censat_summary.tsv.',outer=TRUE,side=1,line=1,cex=.8);dev.off()
}
