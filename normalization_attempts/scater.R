###########################################################
# scater.R
# Author: Aparna Rajpurkar
# single cell RNAseq normalization and analysis pipeline
# Depends on scater, pagoda, scran
###########################################################

#############
# libraries #
#############

library(scater)
library(scran)
library(scde)
library(org.Mm.eg.db)
library(gplots) # heatmap.2
library(RBGL) # community detection
library(matrixStats)

####################
# user-set globals #
####################

MAX_CORES = 35
MIN_READS = 1000000
MIN_FEATURES = 3000
MAX_ERCC_PERC = 50
RUN_NAME = "v3_" # for file name 

plotFigs = TRUE
commDetect = FALSE

#############
# functions #
#############

getPCA <- function(sceset, ntop = 500, exprs = "varinfo.ercc") {
    exprs_mat <- get_exprs(sceset, exprs_value = exprs)
    rv <- matrixStats::rowVars(exprs_mat)
    feature_set <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    exprs_to_plot <- exprs_mat[feature_set,, drop = FALSE]
    exprs_to_plot <- scale(t(exprs_to_plot), scale = TRUE)
    keep_feature <- (matrixStats::colVars(exprs_to_plot) > 0.001)
    keep_feature[is.na(keep_feature)] <- FALSE
    exprs_to_plot <- exprs_to_plot[, keep_feature]
    pca <- prcomp(exprs_to_plot)
    return(pca)
}

communityDetect<-function(var.cor, sig.cor) {
    g <- ftM2graphNEL(cbind(var.cor$gene1, var.cor$gene2)[sig.cor,], 
                          W=NULL, V=NULL, edgemode="undirected")
    cl <- highlyConnSG(g)$clusters
    cl <- cl[order(lengths(cl), decreasing=TRUE)]
    head(cl)

    cl.df<-do.call(rbind, lapply(seq_along(cl), function(i){
                                 data.frame(CLUSTER=i, cl[[i]])
    }))

    write.table(cl.df, file=paste0(RUN_NAME, "communities.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    return(cl.df)
}

#############
# arguments #
#############

args <- c("countfile.tsv")
counts<-read.delim(args[1], header=T)

##############
# clean data #
##############

# clean ENSEMBL gene names
rownames(counts)<-counts[,1]
counts<-counts[,2:ncol(counts)]

# get total number of cells from raw counts
numCells <- ncol(counts)

# set batches
batch<-c(rep("JUNE5",numCells))
june29_ids<-grep("June_29", colnames(counts))
batch[june29_ids]<-"JUNE29"

##########################################################################

# deal with the ENSEMBL gene name issue
ensembl_mouse<-rownames(counts)[grepl("ENSMUSG", rownames(counts))]   
rownames(counts)<-sub("(ENSMUSG\\d+|ERCC.+)\\.*\\d*", "\\1", rownames(counts), perl=TRUE)        

# get ercc names from raw counts
erccs_raw<-rownames(counts[grepl("ERCC",rownames(counts)),])

# set phenotype data (describe cells)
pd<-data.frame(batch=batch)
pd<-new('AnnotatedDataFrame', data=pd)
rownames(pd)<-colnames(counts)

# create SCESet to hold all data
mtecdev<-scater::newSCESet( 
                           countData=counts,
                           phenoData=pd
                           )

# calculate basic QC, setting controls as ERCCs 
mtecdev<-scater::calculateQCMetrics(
                                    mtecdev,
                                    feature_controls=list(ERCC=erccs_raw)
                                    )


# set filters
filter_by_total_counts <- (mtecdev$total_counts > MIN_READS)
filter_by_expr_features <- (mtecdev$total_features > MIN_FEATURES)
filter_by_ERCC <- (mtecdev$pct_counts_feature_controls_ERCC < MAX_ERCC_PERC)

# filter cells
mtecdev$use <- (
                filter_by_total_counts &
                filter_by_expr_features & 
                filter_by_ERCC
                )

# view number of cells kept and dropped
table(mtecdev$use)

# create very basic gene filter
filter_genes <- apply(counts(mtecdev[ , pData(mtecdev)$use]), 1, function(x) length(x[x > 5]) >= 2)

# apply filter to data
fData(mtecdev)$use <- filter_genes

# view number of filtered genes
table(filter_genes)

# view dimensions of dataset (cells and genes)
dim(mtecdev[fData(mtecdev)$use, pData(mtecdev)$use])

# calculate log2 raw counts (DO NOT USE FOR ANYTHING OTHER THAN QC)
set_exprs(mtecdev, "log2_counts") <- log2(counts(mtecdev) + 1)

# create new SCESet with QC'd, filtered data
mtecdev.qc <- mtecdev[fData(mtecdev)$use, pData(mtecdev)$use]

# get all non-control genes (boolean)
endog_genes <- !fData(mtecdev.qc)$is_feature_control
# get all control genes (boolean)
erccs <- fData(mtecdev.qc)$is_feature_control

# get counts for all useable cells and genes
all.co<-counts(mtecdev.qc)
# get names for usable cells
use_cells<-colnames(all.co)

#################
# normalization #
#################

# create a new batch vector for only the usable cells
# to be used with pagoda.varnorm()
knnbatch<-c(rep("JUNE5",length(use_cells)))
june29_ids<-grep("June_29", use_cells)
knnbatch[june29_ids]<-"JUNE29"

# set seed before running knn.error.models
set.seed(0)
# get the error model from pagoda: note MAX_CORES
knn.all<-knn.error.models(all.co, k=ncol(all.co)/4, n.cores = MAX_CORES, min.count.threshold=2, min.nonfailed=5, max.model.plots=10)
# standardize variance to follow chi-square distribution by ERCC norm
varinfo.ercc<-pagoda.varnorm(knn.all, counts=all.co, batch=knnbatch, max.adj.var=5, n.cores=MAX_CORES, plot=FALSE, fit.genes=rownames(mtecdev.qc[erccs,]))

# create a new exprs matrix in the SCESet object
set_exprs(mtecdev.qc, "varinfo.ercc")<-varinfo.ercc$mat

###############
# PCA results #
###############

# since plotPCA() from scater doesn't give back an object, 
# calcuate it separately
pca<-getPCA(mtecdev.qc, ntop=500, exprs="varinfo.ercc")

################
# HVG analysis #
################

# model variance of data with a loess fit 
var.fit<-trendVar(mtecdev.qc[endog_genes,], assay="varinfo.ercc", trend="loess", use.spikes=FALSE, span=0.2)

# decompose variance
var.out<-decomposeVar(mtecdev.qc[endog_genes,], var.fit, assay="varinfo.ercc")

# get HVGs at a cutoff of FDR 5%
hvg.out<-var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]

# order by bio dispersion
hvg.out<-hvg.out[order(hvg.out$bio, decreasing=TRUE),]

# display number of HVGs
nrow(hvg.out)

#############################
# gene correlation analysis #
#############################

# set seed
set.seed(0)
# find correlated gene pairs, using HVGs as subset
var.cor <- correlatePairs(mtecdev.qc, assay="varinfo.ercc", subset.row=rownames(hvg.out))

# see top 10 highly corrleated pairs
head(var.cor)

# only take highly correlated gene pairs of FDR 5%
sig.cor <- var.cor$FDR <= 0.05

# find unique gene names of correlated pairs
chosen <- unique(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor])

# view number of these genes
length(chosen)

# get varinfo for these gene 
norm.exprs <- get_exprs(mtecdev.qc, exprs_value="varinfo.ercc")[chosen,,drop=FALSE]
heat.vals <- norm.exprs - rowMeans(norm.exprs)

cl <- NULL
if (commDetect) {
    cl <- communityDetect(var.cor, sig.cor)
}

#########
# plots #
#########

stopifnot(plotFigs)

pdf(paste0(RUN_NAME, "scater.pdf"))

# histogram of raw total counts
hist(mtecdev$total_counts, breaks=20)
	
# histogram of raw total features
hist(mtecdev$total_features, breaks=20)
	
# scatterplot to view raw ERCC percent as a function of total features
scater::plotPhenoData(
                      mtecdev,
                      aes_string(
                                 x = "total_features",
                                 y = "pct_counts_feature_controls_ERCC",
                                 colour = "batch"
                                 )
                      )
	
# raw count plot of highest-expressed features
scater::plotQC(mtecdev, type = "highest-expression")

# qc'd PCA plot
scater::plotPCA(mtecdev.qc[endog_genes, ],
               ntop = 500,
               colour_by = "batch",
               size_by = "total_features",
               exprs_values = "exprs")

# qc'd tSNE plot
scater::plotTSNE(mtecdev.qc[endog_genes, ],
                 ntop = 500,
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "exprs",
                 rand_seed = 0)

# qc'd plot of top contribution of total features to variance in any PCs
scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "find-pcs",
               variable="total_features",
               exprs_values="exprs"
               )

# qc'd density plot of top likely technical confounders
scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "expl",
               exprs_values="exprs",
               variables=c("total_features",
                           "total_counts",
                           "batch",
                           "pct_counts_feature_controls_ERCC"
                           )
               )

# post-normalization QC contribution plot
scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "find-pcs",
               variable="total_features",
               exprs_values="varinfo.ercc"
               )

# post-normalization density plot
scater::plotQC(mtecdev.qc[endog_genes, ],
               type= "expl",
               exprs_values="varinfo.ercc",
               variables=c("total_features",
                           "total_counts",
                           "batch",
                           "pct_counts_feature_controls_ERCC"
                           )
               )

# post-normalization PCA plot
scater::plotPCA(mtecdev.qc[endog_genes, ],
               ntop = 500,
               colour_by = "batch",
               size_by = "total_features",
               exprs_values = "varinfo.ercc")

# post-normalization tSNE plot
scater::plotTSNE(mtecdev.qc[endog_genes, ],
                 ntop = 500,
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "varinfo.ercc",
                 rand_seed = 0)

# violin plots of top highly variable genes
plotExpression(mtecdev.qc, rownames(hvg.out)[1:25]) 

# PCA plot based on top highly variable genes
scater::plotPCA(mtecdev.qc[endog_genes, ],
                feature_set = chosen,
                colour_by = "batch",
                size_by = "total_features",
                exprs_values = "varinfo.ercc")

# tSNE plot based on top highly variable genes
scater::plotTSNE(mtecdev.qc[endog_genes, ],
                 feature_set = chosen,
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "varinfo.ercc",
                 rand_seed = 0)

# heatmap of top highly variable correlated genes
heat.out <- heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6)

dev.off()
