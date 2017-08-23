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
RUN_NAME = "v3_fc_" # for file name 

plotFigs = TRUE
commDetect = TRUE

#############
# functions #
#############

scdeExprPseudo <- function(models, counts) {
        if(!all(rownames(models) %in% colnames(counts))) { stop("ERROR: provided count data does not cover all of the cells specified in the model matrix") }
    t((t(log(counts[, rownames(models), drop = FALSE] + 1))-models$corr.b)/models$corr.a)
}

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

fitTechnicalNoise=function(nCountsEndo,nCountsERCC=NULL,use_ERCC = TRUE,fit_type = "counts",plot=TRUE){
	
	meansEndo <- rowMeans( nCountsEndo )
	varsEndo <- rowVars( nCountsEndo )
	cv2Endo <- varsEndo / meansEndo^2
	
	meansERCC <- rowMeans( nCountsERCC )
	varsERCC <- rowVars( nCountsERCC )
	cv2ERCC <- varsERCC / meansERCC^2
	
	#Do fitting of technical noise
#	if(!is.null(fit_opts)){
#		if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
#		if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
#	}else{
	mincv2 = 0.3
	quan = 0.8
#	}
	
	#normalised counts (with size factor)
	minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > mincv2 ) ], quan ) )
	useForFitA <- meansERCC >= minMeanForFitA
	fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
			cv2ERCC[useForFitA] )
	
	#4. Transform to log-space and propagate error
	eps=1
	LogNcountsEndo=log10(nCountsEndo+eps)
	dLogNcountsEndo=1/((meansEndo+eps)*log(10))
	var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
	LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
	
	if(plot==TRUE){
		#plot fit
		pdf(file.path(".", "techNoise.pdf"),width=12,height=6)
		
		par(mfrow=c(1,2))
		plot( meansERCC, cv2ERCC, log="xy", col=1+2*useForFitA)
		xg <- 10^seq( -3, 5, length.out=100 )
		lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue' )
		segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
				meansERCC[useForFitA], fitA$fitted.values, col="gray" )
		legend('bottomleft',c('ERCCs used for fit', 'Fit technical noise'),pch=c(1, NA),lty =c(NA,1),col=c('green','blue'),cex=0.8)
		title('Mean-CV2 fit using ERCCs')
		
		#plot fot with all genes
		plot( meansEndo, cv2Endo, log="xy", col=1, xlab = 'Means of human genes', ylab = 'CV2')
		points(meansERCC, cv2ERCC, col='blue', pch=15)
		xg <- 10^seq( -3, 5, length.out=100 )
		lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue',lwd=2 )
		legend('bottomleft',c('Endogenous genes','ERCCs', 'Fit technical noise'),pch=c(1,15, NA),lty =c(NA,NA,1),col=c('black','blue','blue'),cex=0.8)        
		title('Mean-CV2 relationship')
		dev.off()
		
		par(mfrow=c(1,1))
	}
	res = list()
	res$fit = fitA
	res$techNoiseLog = LogVar_techEndo
	res$meansERCC=meansERCC
	return(res)
}



getVariableGenes <- function(nCountsEndo, techNoise, method = "fit", threshold = 0.1, fit_type=NULL,sfEndo=NULL, sfERCC=NULL, plot=T){
	
#	nCountsEndo=nCountsHms
#	
#	method = "fdr" 
#	threshold = 0.01 
#	fit_type="counts"
#	sfEndo=sfHms
#	sfERCC=sfERCC

	fit=techNoise$fit
        meansERCC=techNoise$meansERCC

	if(method=='fdr' & fit_type!="counts"){stop("method='fdr', can only be used with fit_type 'counts'")}
	if(method=='fdr' & (is.null(sfERCC) | is.null(sfEndo))){stop("Please specify sfERCC and sfEndo when using method='fdr'")}
	
	if(method=='fdr'){
		meansEndo <- rowMeans( nCountsEndo )
		varsEndo <- rowVars( nCountsEndo )
		cv2Endo <- varsEndo / meansEndo^2
		
		minBiolDisp <- .5^2
		xi <- mean( 1 / sfERCC )
		m <- ncol(nCountsEndo)

		psia1thetaA <- mean( 1 / sfERCC ) + ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfEndo )

		cv2thA <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
		testDenomA <- ( meansEndo * psia1thetaA + meansEndo^2 * cv2thA ) / ( 1 + cv2thA/m )
		
		pA <- 1 - pchisq( varsEndo * (m-1) / testDenomA, m-1 )
		padjA <- p.adjust( pA, "BH" )
		print( table( padjA < threshold ))
		is_het =  padjA < threshold
		is_het[is.na(is_het)] = FALSE
		
		if(plot==TRUE){
			#cairo_pdf('./tech_noise_genes_sfEndo.pdf',width=4.5,height=4.5)
			
			
			pdf(file.path(".", "variableGenes.pdf"),width=6,height=6)
			plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,250), xlab='Mean Counts', ylab='CV2 Counts')
			xg <- 10^seq( -3, 5, length.out=100 )
			lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='green' )      
			try(points( meansERCC, cv2ERCC, pch=20, cex=1, col="blue" ))
			legend('bottomleft',c('Endo. genes','Var. genes','ERCCs',"Fit"),pch=c(1,1,20,NA),lty = c(NA,NA,NA,1),col=c('black','red','blue', 'green'),cex=0.7)   
			dev.off()
			
			
		}
		
	}
	
	is_het
}


#############
# arguments #
#############

args <- c("all_fc_counts.counts")
counts<-read.delim(args[1], header=T)

##############
# clean data #
##############

# clean ENSEMBL gene names
rownames(counts)<-counts[,4]
counts<-counts[,6:ncol(counts)]

neg_control<-"June29_KW.96.featureCounts"
neg<-grepl(neg_control, colnames(counts))
counts<-counts[,!neg]
weird_cell<-"June5_KW.6.featureCounts"
weird<-grepl(weird_cell, colnames(counts))
counts<-counts[,!weird]

# get total number of cells from raw counts
numCells <- ncol(counts)

# set batches
batch<-c(rep("JUNE5",numCells))
june29_ids<-grep("June29", colnames(counts))
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
#set.seed(0)
# get the error model from pagoda: note MAX_CORES
#knn.all<-knn.error.models(all.co, k=ncol(all.co)/4, n.cores = MAX_CORES, min.count.threshold=10, min.nonfailed=5, max.model.plots=10)
# standardize variance to follow chi-square distribution by ERCC norm
#varinfo.ercc<-pagoda.varnorm(knn.all, counts=all.co, batch=knnbatch, max.adj.var=5, n.cores=MAX_CORES, plot=FALSE, fit.genes=rownames(mtecdev.qc[erccs,]))


#knn.all<-scde.error.models(all.co,n.cores = MAX_CORES, min.count.threshold=10, min.nonfailed=5)
#scde.mat <- scde.expression.magnitude(knn.all, counts = all.co)


print("before norm")

library(DESeq2)
library(statmod)
co<-counts(mtecdev.qc)
splitCounts <- split(as.data.frame(co), grepl("ERCC", rownames(co)))
samplesCounts <- splitCounts[["FALSE"]]

spikeCounts <- splitCounts[["TRUE"]]

#samplesNF <- estimateSizeFactorsForMatrix( samplesCounts )
#spikeNF <- estimateSizeFactorsForMatrix( spikeCounts )

samplesNF <- estimateSizeFactorsForMatrix( co[endog_genes,] )
spikeNF <- estimateSizeFactorsForMatrix( co[erccs,] )
nCo<-co
nCo[endog_genes,] <- t( t(co[endog_genes,]) / samplesNF )
nCo[erccs,] <- t( t(co[erccs,]) / spikeNF )

lognCo<-log(as.matrix(nCo+1))

set_exprs(mtecdev.qc, "varinfo.ercc") <- lognCo

#nCountsSamples <- t( t(samplesCounts) / samplesNF )
#nCountsSpike <- t( t(spikeCounts) / spikeNF )

print("before fitTechnicalNoise")
#my_techNoise <- fitTechnicalNoise(nCountsEndo=nCountsSamples, nCountsERCC=nCountsSpike, fit_type="counts")
my_techNoise <- fitTechnicalNoise(nCountsEndo=nCo[endog_genes,], nCountsERCC=nCo[erccs,], fit_type="counts")

# this command errors with "object "meansERCC" not found". Ignore the error, it's a plotting bug
print("before getVariableGenes")
#varGenes <- getVariableGenes(nCountsSamples, my_techNoise, method="fdr", threshold=1e-3, fit_type="counts", sfEndo=samplesNF, sfERCC=spikeNF)
varGenes <- getVariableGenes(nCo[endog_genes,], my_techNoise, method="fdr", threshold=1e-3, fit_type="counts", sfEndo=samplesNF, sfERCC=spikeNF)
print("after getVariableGenes()")

table(varGenes)
endog_nCo<-nCo[endog_genes,]
varGenesNames<-rownames(endog_nCo[varGenes,])
#rownames(varGenesNames)<-varGenesNames
#varGenesVals <- log(as.matrix(nCountsSamples+1))[varGenes,]








#setSpike(mtecdev.qc)<-"ERCC"
#mtecdev.qc<-computeSpikeFactors(mtecdev.qc, type="ERCC", general.use=FALSE)
#mtecdev.qc<-computeSumFactors(mtecdev.qc, sizes=seq(20,60,5))
#mtecdev.qc<-normalize(mtecdev.qc)
#
#co<-get_exprs(mtecdev.qc, exprs_value="norm_exprs")
#library(SCnorm)
#cond<-rep(c(1),ncol(counts(mtecdev.qc)))
#norm<-SCnorm(Data=counts(mtecdev.qc), Conditions=cond, FilterCellNum=10, useSpikes=TRUE, NCores=20)

# remove batch effect


# create a new exprs matrix in the SCESet object
#set_exprs(mtecdev.qc, "varinfo.ercc")<-ruvg$normalizedCounts

#set_exprs(mtecdev.qc, "varinfo.ercc")<-log2(results(norm))
#set_exprs(mtecdev.qc, "varinfo.ercc")<-get_exprs(mtecdev.qc, exprs_value="exprs")

#co<-counts(mtecdev.qc) + as.integer(1) 
#error.mod<-scde.error.models(counts=co, n.cores=35, threshold.segmentation=TRUE, save.crossfit.plots=FALSE, save.model.plots=FALSE, verbose=1)
#
#o.fpm<-scde.expression.magnitude(error.mod, counts=co)
#set_exprs(mtecdev.qc, "varinfo.ercc")<-o.fpm



print("after norm")

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
#var.fit<-trendVar(mtecdev.qc[endog_genes,], assay="varinfo.ercc", trend="loess", use.spikes=FALSE, span=0.2)

# decompose variance
#var.out<-decomposeVar(mtecdev.qc[endog_genes,], var.fit, assay="varinfo.ercc")

# get HVGs at a cutoff of FDR 5%
#hvg.out<-var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]

# order by bio dispersion
#hvg.out<-hvg.out[order(hvg.out$bio, decreasing=TRUE),]

# display number of HVGs
#nrow(hvg.out)

#############################
# gene correlation analysis #
#############################

# set seed
set.seed(0)
# find correlated gene pairs, using HVGs as subset
#var.cor <- correlatePairs(mtecdev.qc, assay="varinfo.ercc", subset.row=rownames(hvg.out))


var.cor <- correlatePairs(mtecdev.qc, assay="varinfo.ercc", subset.row=varGenesNames, iters=1e6)

# see top 10 highly corrleated pairs
head(var.cor)

# only take highly correlated gene pairs of FDR 5%
sig.cor <- var.cor$FDR <= 0.05

# find unique gene names of correlated pairs
chosen <- unique(c(var.cor$gene1[sig.cor], var.cor$gene2[sig.cor]))
ercc_chosen<-grepl("ERCC",chosen)
chosen<-chosen[!ercc_chosen]

# view number of these genes
length(chosen)

# get varinfo for these gene 
norm.exprs <- get_exprs(mtecdev.qc, exprs_value="varinfo.ercc")[chosen,,drop=FALSE]
heat.vals <- norm.exprs - rowMeans(norm.exprs)

hc <- hclust(dist(t(heat.vals)), method="ward.D2")

cell_clusters <- cutree(hc,2)
cell_cl1<-names(cell_clusters[cell_clusters==1])
cell_cl2<-names(cell_clusters[cell_clusters==2])

norm.exprs.cl1<-norm.exprs[,cell_cl1]
norm.exprs.cl2<-norm.exprs[,cell_cl2]

cl <- NULL
if (commDetect) {
    cl <- communityDetect(var.cor, sig.cor)
    write.table(cl, file=paste0(RUN_NAME, "communities.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

###############
# DE analysis #
###############

# wilcox-test
pVals <- apply(norm.exprs[chosen,], 1, function(x) {
               wilcox.test(x[cell_clusters==1], x[cell_clusters==2], alternative="two.sided")$p.value
                })
pVals <- p.adjust(pVals, method = "fdr")

chosen_p1e2<-names(pVals[pVals<1e-2])
fold<-rowMedians(norm.exprs.cl1[chosen_p1e2,]) - rowMedians(norm.exprs.cl2[chosen_p1e2,])
genes_fold<-cbind.data.frame(chosen_p1e2,fold)
fold2<-genes_fold[genes_fold$fold < -2 | genes_fold$fold > 2,]
fold2<-fold2$chosen_p1e2
write.table(fold2,file="chosen_p1e2_de_clus12.txt",row.names=FALSE, col.names=FALSE,quote=FALSE)


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
#plotExpression(mtecdev.qc, rownames(hvg.out[chosen,])[1:25], exprs_values="varinfo.ercc") 
plotExpression(mtecdev.qc, varGenesNames[1:25], exprs_values="varinfo.ercc") 

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
