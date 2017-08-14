set.seed(0)


# de novo gene clusters
#set.seed(0)
print("pagoda.gene.clusters()")
# REPRODUCIBLE
clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = MAX_CORES, plot = TRUE)

# find top aspects, including both all GO terms and de novo gene clustering
#df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

# get top aspects
print("pagoda.top.aspects()")
#REPRODUCIBLE
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))

# cluster the cells by top aspects
print("pagoda.cluster.cells()")
# REPRODUCIBLE
hc <- pagoda.cluster.cells(tam, varinfo)

# remove redundant pathways / terms
print("pagoda.reduce.loading.redundancy()")
#REPRODUCIBLE
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)
print("pagoda.reduce.loading.redundancy()")
#REPRODUCIBLE
tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

# get names of all pathways collapsed under a single pathway

# > names(tamr2$cnam$"#PC1# GO:0050909")
# get names of all de novo clusters

# tamr2 is in a weird format, the below doesn't work!
# > genes<-names(tamr2$cnam$"#PC1# GO:0050909"[grep("geneCluster", names(tamr2$cnam$"#PC1# GO:0050909"))])
# get names of all GO terms
# > go<-names(tamr2$cnam$"#PC1# GO:0050909"[grep("GO", names(tamr2$cnam$"#PC1# GO:0050909"))])
# remove #PC# 
# > go_new<-sub(".+\\s(.+)", "\\1", go, perl=TRUE)
# > genes_new<-sub(".+\\s(.+)", "\\1", genes, perl=TRUE)

# this works:
# > sep<- do.call("rbind", lapply(tamr2$cnam$"#PC1# GO:0050909", "[[", 1))

# > sep_clean<-sub(".+\\s(.+)", "\\1", sep, perl=TRUE)
# > gene_clusters<-sep_clean[grep("geneCluster", sep_clean)]
# > go<-sep_clean[grep("GO", sep_clean)]

# get genes in de novo clusters:

#> genes_from_clusters<-clpca$clusters[gene_clusters]
#> names(genes_from_clusters)
#[1] "geneCluster.31" "geneCluster.8"  "geneCluster.15" "geneCluster.10"
#[5] "geneCluster.4"  "geneCluster.25"

# no convenient way to do this automatically:

#> write.table(genes_from_clusters$geneCluster.25, "gene_cluster_25.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
#> write.table(genes_from_clusters$geneCluster.4, "gene_cluster_4.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
#> write.table(genes_from_clusters$geneCluster.10, "gene_cluster_10.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
#> write.table(genes_from_clusters$geneCluster.15, "gene_cluster_15.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
#> write.table(genes_from_clusters$geneCluster.8, "gene_cluster_8.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
#> write.table(genes_from_clusters$geneCluster.31, "gene_cluster_31.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# cluster cells for tSNE and final plots
print("pagoda.cluster.cells()")
#REPRODUCIBLE
cell.clustering <- pagoda.cluster.cells(tam,varinfo,include.aspects=TRUE,verbose=TRUE,return.details=T)

# get list of cells in each cluster: 2 clusters
#> groups<-cutree(hc,2)
#> write.table(names(groups[groups==1]), "cluster1_cells.txt")
#> write.table(names(groups[groups==2]), "cluster2_cells.txt", row.names=FALSE, col.names=FALSE,quote=FALSE)



# create tSNE object
print("Rtsne")
set.seed(0)
# REPRODUCIBLE WITH SEED SET BEFORE EACH CALL
tSNE.pagoda <- Rtsne(cell.clustering$distance,is_distance=T,initial_dims=100,perplexity=10)

# set colors and cut tree of cell clustering into clusters
col.cols2 <- rbind(groups = cutree(hc, 2))
col.cols3 <- rbind(groups = cutree(hc, 3))

# silhouette analysis: 2 clusters
print("pam() silhouette() 2")
pam_result_2<-pam(t(cd), 2)
sil_2<-silhouette(pam_result_2)

# silhouette analysis: 3 clusters
print("pam() silhouette() 3")
pam_result_3<-pam(t(cd), 3)
sil_3<-silhouette(pam_result_3)

## Plots ##
print("making plots")
pdf("TEST_REPRODUCIBLE_silhouette_2_clusters.pdf")
plot(sil_2)
dev.off()

pdf("TEST_REPRODUCIBLE_silhouette_3_clusters.pdf")
plot(sil_3)
dev.off()

pdf("TEST_REPRODUCIBLE_pagoda_all_top_aspects_2_clusters.pdf")
pagoda.view.aspects(tam, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(col.cols2))
pagoda.view.aspects(tamr, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(col.cols2))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(col.cols2))

dev.off()

pdf("TEST_REPRODUCIBLE_pagoda_aspects_2_clusters.pdf")
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(col.cols2))
dev.off()

pdf("TEST_REPRODUCIBLE_pagoda_tsne_2_clusters.pdf")
par(mfrow=c(1,1), mar = c(2.5,2.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
plot(tSNE.pagoda$Y,col=adjustcolor(col.cols2,alpha=0.5),cex=1,pch=19,xlab="",ylab="")
dev.off()

pdf("TEST_REPRODUCIBLE_pagoda_aspects_3_clusters.pdf")
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = col.cols3)
dev.off()

pdf("TEST_REPRODUCIBLE_pagoda_tsne_3_clusters.pdf")
par(mfrow=c(1,1), mar = c(2.5,2.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
plot(tSNE.pagoda$Y,col=adjustcolor(col.cols3,alpha=0.5),cex=1,pch=19,xlab="",ylab="")
dev.off()

#### END
