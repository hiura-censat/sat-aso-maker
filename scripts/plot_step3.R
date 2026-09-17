root<-'step3';out<-file.path(root,'plots');dir.create(out,showWarnings=FALSE)
d<-read.delim(file.path(root,'coverage.tsv'))
png(file.path(out,'01_chromosome_coverage.png'),width=2600,height=1300,res=180)
par(mar=c(6,5,4,1));barplot(t(as.matrix(d[,2:4])),beside=TRUE,names.arg=d$chromosome,las=2,col=c('#8794a5','#247BA0','#E49A32'),ylab='Non-reference haplotypes',ylim=c(0,650),main='STEP 3 | Chromosome availability and annotated target windows')
abline(h=573*.85,lty=2,col='gray40');legend('topright',c('T2T available','HSat2 valid windows >0','HSat3 valid windows >0','85% eligibility threshold'),fill=c('#8794a5','#247BA0','#E49A32',NA),border=NA,bty='n',cex=.8)
mtext('Core HSat2: chr7 + chr10 (489 complete haps). Core HSat3: chr5 + chr7 + chr10 + chr20 (442 complete haps).',side=1,line=4.5,cex=.8)
dev.off()
pal<-c('#247BA0','#E49A32','#268878','#AB5BB4','#C54C4C','#6B63A6','#64823F','#333333')
for(family in c('HSat2','HSat3')){
 folder<-file.path(root,'clustering',paste0(family,'_core'));p<-read.delim(file.path(folder,'pca.tsv'));a<-read.delim(file.path(folder,'assignments.tsv'));s<-read.delim(file.path(folder,'cluster_selection.tsv'));ex<-read.delim(file.path(folder,'pca_explained.tsv'));hc<-readRDS(file.path(folder,'hierarchy.rds'));f<-read.delim(file.path(folder,'features.tsv'),check.names=FALSE)
 groups<-a$group[match(p$hap,a$hap)];lev<-sort(unique(groups));col<-ifelse(groups=='unresolved','gray60',pal[match(groups,lev)]);k<-length(lev)
 png(file.path(out,paste0(ifelse(family=='HSat2','02_','03_'),family,'_core_clustering.png')),width=2800,height=2100,res=180)
 par(mfrow=c(2,2),mar=c(4.5,5,3,1),oma=c(3,0,3,0))
 plot(p$PC1,p$PC2,pch=16,cex=.6,col=col,xlab=sprintf('PC1 (%.1f%%)',ex$variance_fraction[1]*100),ylab=sprintf('PC2 (%.1f%%)',ex$variance_fraction[2]*100),main=paste(family,'core-region PCA'))
 legend('topright',lev,col=if(lev[1]=='unresolved')'gray60' else pal[seq_along(lev)],pch=16,bty='n',cex=.8)
 plot(s$k,s$silhouette,type='b',pch=16,col='#247BA0',ylim=c(-.05,1.05),xlab='Number of groups (k)',ylab='Score',main='Partition selection and feature-subsampling stability')
 lines(s$k,s$median_subsample_ARI,type='b',pch=17,col='#268878');lines(s$k,s$median_worst_cluster_Jaccard,type='b',pch=15,col='#E49A32');if(k>1)abline(v=k,lty=2,col='gray40')
 legend('bottomright',c('Silhouette','Median ARI (20 subsamples)','Median worst-group Jaccard'),col=c('#247BA0','#268878','#E49A32'),pch=c(16,17,15),lty=1,bty='n',cex=.7)
 plot(hc,labels=FALSE,hang=-1,main='Ward.D2 hierarchy',xlab=paste(nrow(p),'complete non-reference haplotypes'),sub='',ylab='Height');if(k>1)rect.hclust(hc,k=k,border='gray60')
 xx<-as.matrix(f[,-1]);use<-seq_len(min(20,ncol(xx)));zz<-scale(xx[,use,drop=FALSE]);zz[zz>3]<-3;zz[zz< -3]<- -3;oo<-hc$order
 par(mar=c(4.5,8,3,1));image(seq_len(nrow(xx)),seq_along(use),zz[oo,,drop=FALSE],col=colorRampPalette(c('#247BA0','white','#C54C4C'))(101),zlim=c(-3,3),axes=FALSE,xlab='Haplotypes ordered by hierarchy',ylab='',main='Representative 16-mer profiles')
 axis(2,at=seq_along(use),labels=colnames(xx)[use],las=2,cex.axis=.55);box()
 mtext(paste('STEP 3 |',family,'common-core sequence-pattern groups'),outer=TRUE,side=3,line=1,font=2,cex=1.4)
 mtext('Clustering: log1p(counts per million valid windows), correlated features reduced. Heatmap: per-feature z-scores, blue -3 / white 0 / red +3.',outer=TRUE,side=1,line=.5,cex=.8)
 mtext('Groups are regional patterns, not established whole-genome haplotype groups. Missing core regions are not assigned.',outer=TRUE,side=1,line=1.7,cex=.8)
 dev.off()
}
d<-read.delim(file.path(root,'cluster_concordance.tsv'));v<-d$adjusted_rand_index;cols<-ifelse(grepl('^HSat2',d$comparison),'#247BA0','#E49A32');cols[is.na(v)]<-'gray80'
png(file.path(out,'04_regional_group_concordance.png'),width=2300,height=1600,res=180)
par(mar=c(6.5,10,4,2));bp<-barplot(replace(v,is.na(v),0),names.arg=d$comparison,horiz=TRUE,las=1,col=cols,xlim=c(-.1,1),xlab='Adjusted Rand index on overlapping haplotypes',main='Core-region groups versus other regional analyses',cex.names=.7)
abline(v=0,col='gray50');for(i in which(is.na(v)))text(.02,bp[i],'NA: no comparison partition',adj=0,cex=.6,col='gray40')
mtext('All-available analyses can be affected by missing chromosomes; chromosome-only comparisons are exploratory.',side=1,line=5,cex=.75);dev.off()
if(file.exists(file.path(root,'mismatch','pooled_mismatch_scores.tsv'))){
 d<-read.delim(file.path(root,'mismatch','pooled_mismatch_scores.tsv'));q<-read.delim(file.path(root,'candidates','mismatch_queries.tsv'));n<-nrow(q);labs<-paste0('Q',sprintf('%02d',q$query_index+1),' ',q$family,' ',q$selection_category)
 png(file.path(out,'05_mismatch_specificity.png'),width=3300,height=2200,res=180)
 par(mfrow=c(1,3),mar=c(5,11,4,1),oma=c(3,0,3,0))
 for(field in c('target_vs_background_E','target_vs_other_family_E','strict_passing_nonreference_haps')){
  mat<-matrix(d[[field]],nrow=n,byrow=TRUE)[n:1,,drop=FALSE];hapfield<-field=='strict_passing_nonreference_haps'
  image(0:2,seq_len(n),t(mat),col=if(hapfield)c('#f9ead0','#76b7a8') else c('#e49a95','#f9ead0','#83bbce'),breaks=if(hapfield)c(-.5,515.5,573.5) else c(-100,0,10,100),axes=FALSE,xlab='Maximum mismatches',ylab='',main=switch(field,target_vs_background_E='Target versus background E',target_vs_other_family_E='Target versus other-family E',strict_passing_nonreference_haps='Haps passing strict pan criteria'))
  axis(1,at=0:2);axis(2,at=seq_len(n),labels=rev(labs),las=2,cex.axis=.47);box()
  for(i in seq_len(n))for(j in 1:3)text(j-1,i,if(hapfield)sprintf('%d',mat[i,j]) else sprintf('%.1f',mat[i,j]),cex=.62)
 }
 mtext('STEP 3 | Exhaustive 0-2 mismatch evaluation of 24 representative candidates',outer=TRUE,side=3,line=1,font=2,cex=1.3)
 mtext('E panels: blue >10, beige 0-10, red <0. Hap panel: green >=516/573; strict conditions include family fraction >=99.9%.',outer=TRUE,side=1,line=.5,cex=.8)
 mtext('Both DNA orientations, substitutions only, eligible annotated/background regions. These are not RNA knockdown efficacy scores.',outer=TRUE,side=1,line=1.7,cex=.8);dev.off()
}
