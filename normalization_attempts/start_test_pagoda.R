####################################################
# Author: Aparna Rajpurkar
#
# note: this script is hardcoded to use 35 cores!
#
# This is heavily taken from the official Pagoda/SCDE tutorial
# with modifications for our data
#
# data format: columns = samples/cells, rows = genes, matrix = counts (integer)
# header line for samples, and 1st col is gene names (ENSEMBL)
####################################################

# required libraries
library("scde")
library(Rtsne)
library(org.Mm.eg.db)
library(fpc)
library(cluster)

# set seed
set.seed(0)

# set max cores
MAX_CORES <- 35

# get data
mtec<-read.delim("countfile.tsv", header=T)

# rename rows with 1st col: gene IDs (ENSEMBL)
rownames(mtec) <- mtec[,1]

# remove first column, leaving only count data
mtec<-mtec[,2:ncol(mtec)]

# clean counts
#cd_mtec <- clean.counts(mtec, min.lib.size=1000, min.reads = 1, min.detected = 1)

# remove and save ERCCs
#cd<-cd_mtec[-grep("ERCC", rownames(cd_mtec)), ]
cd<-mtec[-grep("ERCC", rownames(mtec)), ]
erccs<-mtec[grep("ERCC", rownames(mtec)), ]

print("clean.counts()")
cd <- clean.counts(cd, min.lib.size=1000, min.reads = 1, min.detected = 1)

# Deal with the ENSEMBLE ID issue where the GTF contains IDs with ####.## format. Remove everything after and including the period
rownames(cd)<-sub("(.+)\\..*", "\\1", rownames(cd), perl=TRUE)

# get ENTREZ gene ids that correspond the the ENSEMBL IDs
ids <- unlist(lapply(mget(rownames(cd), org.Mm.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]))

# rename rows by ENTREZ IDs
rownames(cd)<-names(ids)

# run KNN error model
print("knn.error.models()")
set.seed(0)
##REPRODUCIBLE
knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = MAX_CORES, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)

# batch effect removal
numcols<-length(colnames(cd))
batch<-c(rep("JUNE5", 126))
june29_ids<-grep("June_29", colnames(cd))
batch[june29_ids]<-"JUNE29"

# get variances using the error model
print("pagoda.varnorm()")
set.seed(0)
# REPRODUCIBLE
varinfo <- pagoda.varnorm(knn, counts = cd, batch=batch, trim = 3/ncol(cd), max.adj.var = 5, n.cores = MAX_CORES, plot = TRUE)

# handle sequencing depth correction
print("pagoda.subtract.aspect()")
# REPRODUCIBLE
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

#genes<-sub("(.+)\\..*", "\\1", rownames(cd_mtec), perl=TRUE)
#ids <- unlist(lapply(mget(genes, org.Mm.egENSEMBL2EG, ifnotfound = NA), function(x) x[1]))

# reverse association of IDs
rids <- names(ids)

# name by IDs
names(rids) <- ids

# get GO pathway terms
go.env <- eapply(org.Mm.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
go.env <- clean.gos(go.env) # Remove GOs with too few or too many genes
go.env <- list2env(go.env)  # Convert to an environment

# using PCA, find main aspects of heterogeneity
print("pagoda.pathway.wPCA()")
#set.seed(0)
## NOT REPRODUCIBLE, even with seed!
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = MAX_CORES)
