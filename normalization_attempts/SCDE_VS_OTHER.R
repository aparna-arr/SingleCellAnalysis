library("scde")

mtec<-read.delim("countfile.tsv", header=T)
#mtec<-read.delim("countfile.tsv", header=T)
rownames(mtec) <- mtec[,1]
mtec<-mtec[,2:ncol(mtec)]

cd_mtec <- clean.counts(mtec, min.lib.size=1000, min.reads = 1, min.detected = 1)

write.table(cd_mtec, file="cd_mtec_after_clean_counts.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

o.ifm <- scde.error.models(counts = cd_mtec, n.cores = 20, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

valid.cells <- o.ifm$corr.a > 0

o.ifm <- o.ifm[valid.cells, ]

o.prior <- scde.expression.prior(models = o.ifm, counts = cd_mtec, length.out = 400, show.plot = FALSE)

o.fpm <- scde.expression.magnitude(o.ifm, counts = cd_mtec)

write.table(o.fpm, file="ofpm_mtec_after_expression_correct.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

# continue with DE functions here

# highly variable genes
nolog_ofpm<-exp(o.fpm)
library(matrixStats)
vars_row <- rowVars(nolog_ofpm)
i1 <-order(vars_row, decreasing=TRUE)
genes <- rownames(nolog_ofpm)
order_genes<-genes[i1]
top_1<-order_genes[1:(length(genes)/100)]
write.table(top_1, file="top_1perc_most_variable_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Skipping to pagoda/gene sets

#knn <- knn.error.models(cd_mtec, k = ncol(cd_mtec)/4, n.cores = 20, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
