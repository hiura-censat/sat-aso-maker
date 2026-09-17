args<-commandArgs(trailingOnly=TRUE);folder<-args[1];reps<-as.integer(args[2]);set.seed(42)
d<-read.delim(file.path(folder,'features.tsv'),check.names=FALSE);haps<-d[[1]];x<-as.matrix(d[,-1]);rownames(x)<-haps
dd<-dist(x);hc<-hclust(dd,method='ward.D2');ks<-2:min(8,nrow(x)-1);labs<-sapply(ks,function(k)cutree(hc,k));colnames(labs)<-paste0('k',ks)
ari<-function(a,b){t<-table(a,b);choose2<-function(v)v*(v-1)/2;n<-sum(t);a2<-sum(choose2(rowSums(t)));b2<-sum(choose2(colSums(t)));exp<-a2*b2/choose2(n);den<-(a2+b2)/2-exp;if(den==0)return(as.numeric(all(outer(a,a,'==')==outer(b,b,'=='))));(sum(choose2(t))-exp)/den}
boot<-matrix(NA_real_,nrow=reps,ncol=length(ks));jac<-boot
if(reps>0)for(b in seq_len(reps)){
 take<-sample(ncol(x),max(1,floor(.8*ncol(x))),replace=FALSE);bh<-hclust(dist(x[,take,drop=FALSE]),method='ward.D2')
 for(j in seq_along(ks)){
  bl<-cutree(bh,ks[j]);base<-labs[,j];boot[b,j]<-ari(base,bl)
  score<-sapply(unique(base),function(g)max(sapply(unique(bl),function(h)sum(base==g&bl==h)/sum(base==g|bl==h))))
  jac[b,j]<-min(score)
 }
}
stats<-data.frame(k=ks,silhouette=sapply(seq_along(ks),function(j)mean(cluster::silhouette(labs[,j],dd)[,3])),min_cluster_size=sapply(seq_along(ks),function(j)min(table(labs[,j]))),median_subsample_ARI=if(reps>0)apply(boot,2,median) else NA_real_,median_worst_cluster_Jaccard=if(reps>0)apply(jac,2,median) else NA_real_)
write.table(data.frame(hap=haps,labs,check.names=FALSE),file.path(folder,'partitions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(stats,file.path(folder,'cluster_selection.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
if(reps>0)write.table(data.frame(replicate=seq_len(reps),boot,check.names=FALSE),file.path(folder,'feature_subsample_ARI.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(hc,file.path(folder,'hierarchy.rds'))
write.table(data.frame(order=seq_along(hc$order),hap=haps[hc$order]),file.path(folder,'dendrogram_order.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
