root<-'step3';out<-file.path(root,'plots');dir.create(out,showWarnings=FALSE,recursive=TRUE)
d<-read.delim(file.path(root,'coverage.tsv'))
png(file.path(out,'01_chromosome_coverage.png'),width=2600,height=1300,res=180)
par(mar=c(6,5,4,1))
barplot(t(as.matrix(d[,2:4])),beside=TRUE,names.arg=d$chromosome,las=2,col=c('#8794a5','#247BA0','#E49A32'),ylab='Non-reference haplotypes',ylim=c(0,650),main='STEP 3 | Measurable haplotypes by chromosome')
legend('topright',c('T2T available','HSat2 measurable','HSat3 measurable'),fill=c('#8794a5','#247BA0','#E49A32'),bty='n')
dev.off()
groups<-read.delim(file.path(root,'hap_chromosome_groups.tsv'))
cases<-unique(groups[,c('family','chromosome')]);cases<-cases[order(cases$family,match(cases$chromosome,d$chromosome)),]
stats<-data.frame(family=character(),chromosome=character(),measurable=integer(),assigned=integer(),groups=integer())
for(i in seq_len(nrow(cases))){
 x<-subset(groups,family==cases$family[i] & chromosome==cases$chromosome[i] & is_reference==0)
 stats<-rbind(stats,data.frame(family=cases$family[i],chromosome=cases$chromosome[i],measurable=sum(x$measurable==1),assigned=sum(x$group!='NA'),groups=length(unique(x$group[x$group!='NA']))))
 folder<-file.path(root,'clustering',paste0(cases$family[i],'_',cases$chromosome[i]))
 if(!file.exists(file.path(folder,'pca.tsv')))next
 p<-read.delim(file.path(folder,'pca.tsv'));a<-read.delim(file.path(folder,'assignments.tsv'))
 if(!'PC2'%in%names(p))next
 label<-a$group[match(p$hap,a$hap)]
 lev<-sort(unique(label));palette<-setNames(rainbow(length(lev)),lev);if('unresolved'%in%lev)palette['unresolved']<-'gray65'
 png(file.path(out,paste0(cases$family[i],'_',cases$chromosome[i],'_clustering.png')),width=1900,height=1400,res=180)
 plot(p$PC1,p$PC2,pch=16,cex=.65,col=palette[label],xlab='PC1',ylab='PC2',main=paste('Independent chromosome cluster:',cases$family[i],cases$chromosome[i]))
 legend('topright',lev,pch=16,col=palette[lev],bty='n')
 dev.off()
}
write.table(stats,file.path(root,'plots','chromosome_group_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
png(file.path(out,'02_chromosome_group_summary.png'),width=2600,height=1500,res=180)
par(mfrow=c(2,1),mar=c(5,5,4,1))
for(fam in c('HSat2','HSat3')){
 x<-stats[stats$family==fam,]
 barplot(x$groups,names.arg=x$chromosome,las=2,col=if(fam=='HSat2')'#247BA0' else '#E49A32',ylab='Supported groups',main=paste(fam,'independent chromosome partitions'))
}
dev.off()
