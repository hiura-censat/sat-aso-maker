root<-'step4'; out<-file.path(root,'plots'); dir.create(out,showWarnings=FALSE)
pal<-c(HSat2='#247BA0',HSat3='#E49A32')

d<-read.delim(file.path(root,'ranked','mismatch_evaluated_ranked.tsv'),check.names=FALSE)
png(file.path(out,'01_candidate_score_landscape.png'),width=2600,height=1600,res=180)
par(mfrow=c(1,2),mar=c(6,5,4,1),oma=c(4,0,3,0))
plot(d$exact_composite_score,d$mismatch_score,pch=21,bg=pal[d$target_family],col='white',cex=1.2,xlab='Exact-match composite score',ylab='0-2 mismatch score',main='72 candidates with mismatch evaluation')
legend('bottomright',names(pal),pt.bg=pal,pch=21,bty='n');grid(col='gray90')
ord<-c('HSat2.chromosome','HSat3.chromosome','HSat2.pan','HSat3.pan','HSat2.regional_group','HSat3.regional_group'); grp<-interaction(d$target_family,d$candidate_type)
boxplot(split(d$final_score,grp)[ord],names=c('H2 chr','H3 chr','H2 pan','H3 pan','H2 group','H3 group'),las=2,col=c(pal['HSat2'],pal['HSat3'],pal['HSat2'],pal['HSat3'],pal['HSat2'],pal['HSat3']),ylab='Final prioritization score',xlab='',main='Final score by family and intended use')
mtext('STEP 4 | Integrated exact-match, sequence-quality proxy, and mismatch evidence',outer=TRUE,side=3,line=1,font=2,cex=1.25)
mtext('Scores prioritize experimental follow-up; they are not measured ASO efficacy.',outer=TRUE,side=1,line=.5,cex=.8);dev.off()

p<-read.delim(file.path(root,'mismatch','pooled_mismatch_scores.tsv'))
top<-d[order(-d$final_score),][1:min(30,nrow(d)),]; ids<-top$kmer_id
z<-p[p$kmer_id%in%ids,]; z<-z[order(match(z$kmer_id,ids),z$max_mismatches),]
labs<-paste0(substr(top$kmer_id,5,12),' ',top$target_family,' ',top$candidate_type)
png(file.path(out,'02_top_candidate_mismatch_specificity.png'),width=3200,height=2200,res=180)
par(mfrow=c(1,3),mar=c(5,12,4,1),oma=c(3,0,3,0))
for(field in c('target_vs_background_E','target_vs_other_family_E','strict_passing_nonreference_haps')){
 mat<-matrix(z[[field]],nrow=length(ids),byrow=TRUE)[length(ids):1,,drop=FALSE]; hp<-field=='strict_passing_nonreference_haps'
 image(0:2,seq_len(length(ids)),t(mat),col=if(hp)c('#f4d8d5','#eaf3dd','#68a88d') else c('#d7786f','#f5e8c8','#6aaac3'),breaks=if(hp)c(-.5,99.5,515.5,573.5) else c(-100,0,10,100),axes=FALSE,xlab='Maximum mismatches',ylab='',main=switch(field,target_vs_background_E='Target vs background E',target_vs_other_family_E='Target vs other-family E',strict_passing_nonreference_haps='Strict passing haps'))
 axis(1,at=0:2);axis(2,at=seq_len(length(ids)),labels=rev(labs),las=2,cex.axis=.55);box()
 for(i in seq_len(nrow(mat)))for(j in 1:3)text(j-1,i,if(hp)sprintf('%d',mat[i,j]) else sprintf('%.1f',mat[i,j]),cex=.62)
}
mtext('STEP 4 | Top 30 candidates: exhaustive substitution-neighborhood specificity',outer=TRUE,side=3,line=1,font=2,cex=1.25)
mtext('Both orientations; radius 0, 1, and 2 substitutions. E >10 and >=516/573 strict-passing haps are stringent pan-family criteria.',outer=TRUE,side=1,line=.5,cex=.75);dev.off()

tiers<-c('A_le2mm','B_le1mm','C_exact','D_fails_exact');tab<-table(factor(d$category_robustness_tier,levels=tiers),d$target_family)
png(file.path(out,'03_robustness_tiers.png'),width=2200,height=1500,res=180)
par(mar=c(7,5,4,1));barplot(t(tab),beside=TRUE,col=pal,names.arg=c('A: <=2 mm','B: <=1 mm','C: exact','D: fails exact'),las=2,ylab='Mismatch-evaluated candidates',main='Intended-use robustness tier',ylim=c(0,max(tab)+5));legend('topright',names(pal),fill=pal,bty='n')
mtext('Category criteria: pan-family specificity, regional core-group contrast, or chromosome concentration.',side=1,line=5.3,cex=.8);dev.off()
