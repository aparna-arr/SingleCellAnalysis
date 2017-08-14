library(monocle)

clus2 <- cutree(hc,2)
set.seed(0)
pam_tsne_clusters<-pamk(tSNE.pagoda$Y)$pamobject$clustering
set.seed(0)
kmeans_tsne<-kmeans(tSNE.pagoda$Y, 5)$cluster

sep<- do.call("rbind", lapply(tamr2$cnam$"#PC1# GO:0050909", "[[", 1)) 
sep_clean<-sub(".+\\s(.+)", "\\1", sep, perl=TRUE) 
gene_clusters<-sep_clean[grep("geneCluster", sep_clean)]
genes_from_clusters<-clpca$clusters[gene_clusters]  
genes_from_clusters<-matrix(unlist(genes_from_clusters), ncol=1, byrow=TRUE)[,1]

go<-sep_clean[grep("GO", sep_clean)]  
genes_from_go<-unique(get(go, go.env))

all_genes_of_top_aspect<-unique(c(genes_from_go, genes_from_clusters))

subset_taste_genes<-unique(get("GO:0050909", go.env))

pd <- data.frame(batch=batch, pagoda=factor(clus2[colnames(cd)]), pam_tsne=factor(pam_tsne_clusters), kmeans_tsne=factor(kmeans_tsne)) 

pd <- new('AnnotatedDataFrame', data = pd) 
#fd <- data.frame(taste_genes=all_genes_of_top_aspect)

fd<-data.frame(allgenes=rownames(cd))
rownames(fd)<-rownames(cd)
fd<-new('AnnotatedDataFrame', data=fd)

cd.mat<-as.matrix(cd)
#dataset<-newCellDataSet(cd.mat,
#                        phenoData = pd,
#                        featureData = fd,
#                        lowerDetectionLimit = 0.5,
#                        expressionFamily = negbinomial.size()
#                        )

dataset<-new('CellDataSet', exprs=cd.mat, phenoData=pd, featureData=fd, expressionFamily=negbinomial.size(), lowerDetectionLimit=0.5)

dataset <- estimateSizeFactors(dataset)
dataset <- estimateDispersions(dataset)

print("setOrderingFilter")
all_dataset<-setOrderingFilter(dataset, all_genes_of_top_aspect)

print("reduceDimension")
all_dataset<-reduceDimension(all_dataset, max_components=2, method="DDRTree")
set.seed(0)
print("orderCells")
all_dataset<-orderCells(all_dataset)

taste_dataset<-setOrderingFilter(dataset, subset_taste_genes)
print("reduceDimension")
taste_dataset<-reduceDimension(taste_dataset, max_components=2, method="DDRTree")
set.seed(0)
print("orderCells")
taste_dataset<-orderCells(taste_dataset)

#print("reduceDimension")
#dataset<-reduceDimension(dataset, max_components=2, method="DDRTree")
#set.seed(0)
#print("orderCells")
#dataset<-orderCells(dataset)

pdf("spanning_tree.pdf")
plot(tSNE.pagoda$Y,cex=1,pch=19,xlab="",ylab="")
plot(tSNE.pagoda$Y,col=factor(pam_tsne_clusters),cex=1,pch=19,xlab="",ylab="", main="pamk")
plot(tSNE.pagoda$Y,col=factor(kmeans_tsne),cex=1,pch=19,xlab="",ylab="", main="kmeans")
plot(tSNE.pagoda$Y,col=factor(clus2[colnames(cd)]),cex=1,pch=19,xlab="",ylab="", main="pagoda")

plot_cell_trajectory(all_dataset)
plot_cell_trajectory(all_dataset, color_by="pagoda")
plot_cell_trajectory(all_dataset, color_by="batch")
plot_cell_trajectory(all_dataset, color_by="pam_tsne")
plot_cell_trajectory(all_dataset, color_by="kmeans_tsne")

plot_cell_trajectory(taste_dataset)
plot_cell_trajectory(taste_dataset, color_by="pagoda")
plot_cell_trajectory(taste_dataset, color_by="batch")
plot_cell_trajectory(taste_dataset, color_by="pam_tsne")
plot_cell_trajectory(taste_dataset, color_by="kmeans_tsne")
dev.off()
