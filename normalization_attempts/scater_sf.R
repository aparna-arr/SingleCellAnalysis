library(scater)
library(knitr)
library(RUVSeq)
library(scran)
library(scde)
library(org.Mm.eg.db)
library(gplots)

MAX_CORES = 35

args <- c("countfile.tsv",1e-3)
counts<-read.delim(args[1], header=T)
fdr_thresh<-args[2]

#rownames(counts)<-sub("(.+)\\..*", "\\1", counts[,1], perl=TRUE)        
rownames(counts)<-counts[,1]
counts<-counts[,2:ncol(counts)]

#splitCounts <- split(as.data.frame(counts), grepl("ERCC", rownames(counts)))
#samplesCounts <- splitCounts[["FALSE"]]

#validCells<-colSums(samplesCounts>=20) > 1000
#samplesCounts<-samplesCounts[,validCells]

ensembl_mouse<-rownames(counts)[grepl("ENSMUSG", rownames(counts))]   
rownames(counts)<-sub("(ENSMUSG\\d+|ERCC.+)\\.*\\d*", "\\1", rownames(counts), perl=TRUE)        

## Winsorize matrix ##
#models<-knn.error.models(counts, k=ncol(counts)/4, n.cores = MAX_CORES, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
#trim<-3/ncol(counts)
#
#fpm <- t((t(log(counts))-models$corr.b)/models$corr.a)
#tfpm <- winsorize.matrix(fpm, trim)
#rn <- rownames(counts)
#cn <- colnames(counts)
#counts <- round(exp(t(t(tfpm)*models$corr.a+models$corr.b)))
#counts[counts<0] <- 0
#rownames(counts) <- rn
#colnames(counts) <- cn
#rm(fpm, tfpm)
#counts <- counts[rowSums(counts) > 0, ] # omit genes without any data after Winsorization


#spikeCounts <- splitCounts[["TRUE"]]
#spikeCounts <- spikeCounts[,validCells]

numCells <- ncol(counts)
batch<-c(rep("JUNE5",numCells))
june29_ids<-grep("June_29", colnames(counts))
batch[june29_ids]<-"JUNE29"

pd<-data.frame(batch=batch)
pd<-new('AnnotatedDataFrame', data=pd)
rownames(pd)<-colnames(counts)
fd<-data.frame(allgenes=rownames(counts))
rownames(fd)<-rownames(counts)
fd<-new('AnnotatedDataFrame', data=fd)

#samplesCounts<-as.matrix(samplesCounts)
erccs_raw<-rownames(counts[grepl("ERCC",rownames(counts)),])
mtecdev<-scater::newSCESet( 
                           countData=counts,
                           phenoData=pd
                           )

mtecdev<-scater::calculateQCMetrics(
                                    mtecdev,
                                    feature_controls=list(ERCC=erccs_raw)
                                    )


filter_by_total_counts <- (mtecdev$total_counts > 1000000)

filter_by_expr_features <- (mtecdev$total_features > 3000)

filter_by_ERCC <- (mtecdev$pct_counts_feature_controls_ERCC < 50)

mtecdev$use <- (
                filter_by_total_counts &
                filter_by_expr_features & 
                filter_by_ERCC
                )

table(mtecdev$use)

mtecdev$use_default <- (

                        !mtecdev$filter_on_total_features &
                        !mtecdev$filter_on_total_counts &
                        !mtecdev$filter_on_pct_counts_feature_controls_ERCC &
                        !mtecdev$is_cell_control
                )
table(mtecdev$use_default)

filter_genes <- apply(counts(mtecdev[ , pData(mtecdev)$use]), 1, function(x) length(x[x > 5]) >= 2)

fData(mtecdev)$use <- filter_genes
table(filter_genes)

dim(mtecdev[fData(mtecdev)$use, pData(mtecdev)$use])

set_exprs(mtecdev, "log2_counts") <- log2(counts(mtecdev) + 1)

mtecdev.qc <- mtecdev[fData(mtecdev)$use, pData(mtecdev)$use]
endog_genes <- !fData(mtecdev.qc)$is_feature_control
erccs <- fData(mtecdev.qc)$is_feature_control

#er<-counts(mtecdev.qc[erccs,])

#set.seed(0)
#knner<-knn.error.models(er, k=ncol(er)/4, n.cores = MAX_CORES, min.count.threshold=2, min.nonfailed=5, max.model.plots=10, min.size.entries=20)
#varinfoer<-pagoda.varnorm(knner, counts=er, max.adj.var=5, n.cores=MAX_CORES, plot=FALSE)

#co<-counts(mtecdev.qc[endog_genes,])
#use_cells<-colnames(co)
#knnbatch<-c(rep("JUNE5",length(use_cells)))
#june29_ids<-grep("June_29", use_cells)
#knnbatch[june29_ids]<-"JUNE29"

#set.seed(0)
#knn<-knn.error.models(co, k=ncol(co)/4, n.cores = MAX_CORES, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)

#varinfo<-pagoda.varnorm(knn, counts=co, batch=knnbatch, trim=3/ncol(co), max.adj.var=5, n.cores=MAX_CORES, plot=FALSE)
#varinfo<-pagoda.varnorm(knn, counts=co, batch=knnbatch, max.adj.var=5, n.cores=MAX_CORES, plot=FALSE)
#varinfo<-pagoda.subtract.aspect(varinfo, colSums(co[, rownames(knn)]>0))

#combined<-rbind(varinfo$mat, varinfoer$mat)

#exprs_rows <- rownames(get_exprs(mtecdev.qc, exprs_value="exprs"))
#combined2 <- combined[exprs_rows,,drop=FALSE]

#set_exprs(mtecdev.qc, "varinfo")<-combined2

#all.co<-counts(mtecdev.qc)
#set.seed(0)
#knn.all<-knn.error.models(all.co, k=ncol(all.co)/4, n.cores = MAX_CORES, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
#exprs<-pagoda.varnorm(knn.all, counts=all.co, batch=knnbatch, max.adj.var=5, n.cores=MAX_CORES, plot=FALSE, fit.genes=rownames(mtecdev.qc[erccs,]))
#
#set_exprs(mtecdev.qc, "exprs")<-exprs$mat

#ruvg <- RUVg(counts(mtecdev.qc), erccs, k = 1)
#set_exprs(mtecdev.qc, "ruvg1") <- ruvg$normalizedCounts
#set_exprs(mtecdev.qc, "log2_ruvg1") <- log2(ruvg$normalizedCounts + 1)
#ruvg <- RUVg(counts(mtecdev.qc), erccs, k = 2)
#set_exprs(mtecdev.qc, "ruvg2") <- ruvg$normalizedCounts
#set_exprs(mtecdev.qc, "log2_ruvg2") <- log2(ruvg$normalizedCounts + 1)

#co<-counts(mtecdev.qc)

#ruvg <- RUVg(co, erccs, k = 1)
#set_exprs(mtecdev.qc, "ruvg1") <- ruvg$normalizedCounts
#set_exprs(mtecdev.qc, "log2_ruvg1") <- log2(ruvg$normalizedCounts + 1)
#ruvg <- RUVg(co, erccs, k = 2)
#set_exprs(mtecdev.qc, "ruvg2") <- ruvg$normalizedCounts
#set_exprs(mtecdev.qc, "log2_ruvg2") <- log2(ruvg$normalizedCounts + 1)


setSpike(mtecdev.qc)<-"ERCC"
mtecdev.qc <- computeSpikeFactors(mtecdev.qc, type="ERCC", general.use=FALSE)
mtecdev.qc <- computeSumFactors(mtecdev.qc, sizes=seq(20, 60, 5))
#summary(sizeFactors(mtecdev.qc))


#mtecdev.new <- normaliseExprs(
#                         mtecdev.qc,
#                         method = "RLE", 
#                         feature_set = endog_genes
#                             )

mtecdev.new<-normalize(mtecdev.qc)

## HVG ANALYSIS ##
var.fit<-trendVar(mtecdev.qc[endog_genes,], assay="exprs", trend="loess", use.spikes=FALSE, span=0.2)
var.out<-decomposeVar(mtecdev.qc[endog_genes,], var.fit, assay="exprs")
hvg.out<-var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out<-hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)

## correlated analysis ##
set.seed(0)
var.cor <- correlatePairs(mtecdev.qc, assay="exprs", subset.row=rownames(hvg.out))
head(var.cor)
sig.cor <- var.cor$FDR <= 0.05

chosen <- unique(c(var.cor$gene1[sig.cor],var.cor$gene2[sig.cor]))
exprs <- get_exprs(mtecdev.qc, exprs_value="exprs")[chosen,,drop=FALSE]
heat.vals <- exprs - rowMeans(exprs)



pdf("scater_hist.pdf")
hist(mtecdev$total_counts, breaks=20)
hist(mtecdev$total_features, breaks=20)
scater::plotPhenoData(
                      mtecdev,
                      aes_string(
                                 x = "total_features",
                                 y = "pct_counts_feature_controls_ERCC",
                                 colour = "batch"
                                 )
                      )

scater::plotQC(mtecdev, type = "highest-expression")

scater::plotPCA(mtecdev.qc[endog_genes, ],
               ntop = 500,
               colour_by = "batch",
               size_by = "total_features",
               exprs_values = "exprs")

scater::plotTSNE(mtecdev.qc[endog_genes, ],
                 ntop = 500,
#                 perplexity = 10,
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "exprs",
                 rand_seed = 0)
scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "find-pcs",
               variable="total_features",
               exprs_values="exprs"
               )

scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "expl",
               exprs_values="exprs",
               variables=c("total_features",
                           "total_counts",
                           "batch",
                           "pct_counts_feature_controls_ERCC"
                           )
               )
scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "find-pcs",
               variable="total_features",
               exprs_values="exprs"
               )

scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "expl",
               exprs_values="exprs",
               variables=c("total_features",
                           "total_counts",
                           "batch",
                           "pct_counts_feature_controls_ERCC"
                           )
               )

scater::plotPCA(mtecdev.qc[endog_genes, ],
               ntop = 500,
               colour_by = "batch",
               size_by = "total_features",
               exprs_values = "exprs")

scater::plotTSNE(mtecdev.qc[endog_genes, ],
                 ntop = 500,
#                 perplexity = 10,
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "exprs",
                 rand_seed = 0)

plotExpression(mtecdev.qc, head(rownames(hvg.out))) 

scater::plotPCA(mtecdev.qc[endog_genes, ],
                feature_set = rownames(hvg.out),
                colour_by = "batch",
                size_by = "total_features",
                exprs_values = "exprs")

scater::plotTSNE(mtecdev.qc[endog_genes, ],
                 feature_set = rownames(hvg.out),
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "exprs",
                 rand_seed = 0)

heat.out <- heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6)
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "find-pcs",
#               variable="total_features",
#               exprs_values="exprs"
#               )
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "find-pcs",
#               variable="total_features",
#               exprs_values="ruvg1"
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "find-pcs",
#               variable="total_features",
#               exprs_values="log2_ruvg1"
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "find-pcs",
#               variable="total_features",
#               exprs_values="ruvg2"
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "find-pcs",
#               variable="total_features",
#               exprs_values="log2_ruvg2"
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "expl",
#               exprs_values="ruvg1",
#               variables=c("total_features",
#                           "total_counts",
#                           "batch",
#                           "pct_counts_feature_controls_ERCC"
#                           )
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "expl",
#               exprs_values="log2_ruvg1",
#               variables=c("total_features",
#                           "total_counts",
#                           "batch",
#                           "pct_counts_feature_controls_ERCC"
#                           )
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "expl",
#               exprs_values="ruvg2",
#               variables=c("total_features",
#                           "total_counts",
#                           "batch",
#                           "pct_counts_feature_controls_ERCC"
#                           )
#               )
#
#scater::plotQC(mtecdev.qc[endog_genes, ],
#               type= "expl",
#               exprs_values="log2_ruvg2",
#               variables=c("total_features",
#                           "total_counts",
#                           "batch",
#                           "pct_counts_feature_controls_ERCC"
#                           )
#               )
#
#scater::plotPCA(mtecdev.qc[endog_genes, ],
#               ntop = 500,
#               colour_by = "batch",
#               size_by = "total_features",
#               exprs_values = "ruvg1")
#
#scater::plotTSNE(mtecdev.qc[endog_genes, ],
#                 ntop = 500,
##                 perplexity = 10,
#                 colour_by = "batch",
#                 size_by = "total_features",
#                 exprs_values = "ruvg1",
#                 rand_seed = 0)
#
#scater::plotPCA(mtecdev.qc[endog_genes, ],
#               ntop = 500,
#               colour_by = "batch",
#               size_by = "total_features",
#               exprs_values = "log2_ruvg1")
#
#scater::plotTSNE(mtecdev.qc[endog_genes, ],
#                 ntop = 500,
##                 perplexity = 10,
#                 colour_by = "batch",
#                 size_by = "total_features",
#                 exprs_values = "log2_ruvg1",
#                 rand_seed = 0)
#
#scater::plotPCA(mtecdev.qc[endog_genes, ],
#               ntop = 500,
#               colour_by = "batch",
#               size_by = "total_features",
#               exprs_values = "ruvg2")
#
#scater::plotTSNE(mtecdev.qc[endog_genes, ],
#                 ntop = 500,
##                 perplexity = 10,
#                 colour_by = "batch",
#                 size_by = "total_features",
#                 exprs_values = "ruvg2",
#                 rand_seed = 0)
#
#scater::plotPCA(mtecdev.qc[endog_genes, ],
#               ntop = 500,
#               colour_by = "batch",
#               size_by = "total_features",
#               exprs_values = "log2_ruvg2")

#scater::plotTSNE(mtecdev.qc[endog_genes, ],
#                 ntop = 500,
#                 perplexity = 10,
#                 colour_by = "batch",
#                 size_by = "total_features",
#                 exprs_values = "log2_ruvg2",
#                 rand_seed = 0)

dev.off()
