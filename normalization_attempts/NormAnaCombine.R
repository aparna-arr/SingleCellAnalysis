# TODO: Add comment
# 
# Author: wave
###############################################################################


library("Biostrings")
library("ShortRead")

library(genefilter)
library(statmod)
require(ggplot2)
library(gplots)
require(DESeq2)

#serverDir="/g/steinmetz/wuwei/"
serverDir="/home/wuwei5/"

#projectDir=file.path(serverDir,"projects/dmel/")

projectDir=file.path(serverDir,"projects/neuron3d/2016-08-24")

dataDir=file.path(projectDir,"data/")
resDir=file.path(projectDir,"res/combinded")


###############################
#
# combine data..
#

combHTCS=function(){
	
	exps=c("2016-08-04","2016-08-24")
	
	infDat=lapply(exps,function(x){
				
				print(paste(x,"is loading"))
				
				exprsDataDir=paste(projectDir,"..",x,"data/rdata",sep="/")
				
				
				load(file.path(exprsDataDir,"sampleInfo.rda"))
				sampleInfo$projectDir=projectDir
				sampleInfo$info$sc=sampleInfo$info$cellNumber %in% c("Single Cell","Single cell")
				
				sinfos=sampleInfo$info
				
				sinfos$seqTime=x
				
				
				load(file.path(exprsDataDir,"htcs.rda"))
				
				return(list(info=sinfos,htcs=htcs))
				
			})
	
	
	lapply(infDat,function(x) colnames(x$info))
	
	infos=do.call(rbind,lapply(infDat,function(x) subset(x$info)))
	
	
	clbind<-function(x) Reduce(cbind,x)
	htcs=clbind(lapply(infDat,function(x) x$htcs))
	
	table(colnames(htcs)==infos$sample)
	
	save(infos,htcs,file=file.path(dataDir,"rdata/htcsCombine.rda"))
	
}

############################
#
# change rownames to gene names..
#

id2genenames=function(){
	
	
	load(file.path(dataDir,"rdata/htcsCombine.rda"))
	
	dat=htcs[1:(nrow(htcs)-5),]
	
	geneTypes <- factor( c( ENSG="ENSG", ERCC="ERCC" )[
					substr( rownames(dat), 1, 4 ) ] )
	
	load(file.path(serverDir,"projects/referenceSequences/rdata/human/gtf.rda"))
	genenames=rownames(dat)
	hgIds=geneTypes=="ENSG"
	genenames[hgIds]=as.character(gtf$gene_name[match(rownames(dat)[hgIds],gtf$gene_id)])
	
	#fix the problem where two genes have identical names
	which(genenames=="TUBB3")
	which(genenames=="RBFOX3")
	genenames[18079]="TUBB3-2"
	genenames[58392]="RBFOX3-2"

	dIds=duplicated(genenames)
	
	ngn=genenames
	ngn[dIds]=paste(ngn[dIds],"--d",1:length(dIds),sep="")
	
	gdat=dat
	rownames(gdat)=ngn
	
	save(infos,dat,gdat,file=file.path(dataDir,"rdata/gdats.rda"))
	
}

###############################
#
# load combined ht-count files...
#

load(file.path(dataDir,"rdata/htcsCombine.rda"))


dat=htcs[1:(nrow(htcs)-5),]

geneTypes <- factor( c( ENSG="ENSG", ERCC="ERCC" )[
				substr( rownames(dat), 1, 4 ) ] )

###############
#
cc=dat[ which( geneTypes=="ENSG" ),]
ii=which (colSums(cc>=20) <1000 & infos$sc)
colSums(cc>=20)[ii]

#filter out empty wells; only single cells
scIds= colSums(cc>=20) >1000 & infos$sc

# calculate normalisation for counts
countsHms <- dat[ which( geneTypes=="ENSG" ), scIds]
countsERCC <- dat[ which( geneTypes=="ERCC" ), scIds]
lengthshms <- dat[ which( geneTypes=="ENSG" ), 1 ]
lengthsERCC <- dat[ which( geneTypes=="ERCC" ), 1 ]


sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfHms <- estimateSizeFactorsForMatrix( countsHms )

#sfHms <- sfERCC #also use ERCC size factor for endogenous genes

#normalise read counts
nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsHms <- t( t(countsHms) / sfHms )




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
		pdf(file.path(resDir, "techNoise.pdf"),width=12,height=6)
		
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
			
			
			pdf(file.path(resDir, "variableGenes.pdf"),width=6,height=6)
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


techNoise = fitTechnicalNoise(nCountsHms,nCountsERCC=nCountsERCC, fit_type = 'counts')  
#is_het = getVariableGenes(nCountsHms, techNoise, method = "fdr", threshold = 0.01, fit_type="counts",sfEndo=sfHms, sfERCC=sfERCC)
is_het = getVariableGenes(nCountsHms, techNoise, method = "fdr", threshold = 1e-3, fit_type="counts",sfEndo=sfHms, sfERCC=sfERCC)

table(is_het)

mat <- log10(as.matrix(nCountsHms)+1)[is_het,]


save(mat,is_het,)



breaks()


#base.pca <- prcomp(t(nCountsHms))
base.pca <- prcomp(t(mat))

sg=factor(infos$type[match(colnames(mat),infos$sample)])
sgCol <- rainbow(length(levels(sg)))[sg]

sg2=factor(infos$seqTime[match(colnames(mat),infos$sample)])
sgCol2 <- rainbow(length(levels(sg2)))[sg2]
pch2=as.numeric(sg2)+15

pdf(file.path(resDir, "sampleSeqtimepca.pdf"),width=6,height=6)
plot(base.pca$x[,1], base.pca$x[,2], col=sgCol, pch=pch2, main='PCA',xlab="PCA1",ylab="PCA2")
legend('bottomleft',c('P cells','SP cells'),pch=c(16,16),col=c(sgCol[1],sgCol[length(sgCol)]))  
legend('bottomright',c('first seq','second seq'),pch=c(16,17),col=1)  
dev.off()


pdf(file.path(resDir, "samplepca.pdf"),width=6,height=6)
plot(base.pca$x[,1], base.pca$x[,2], col=sgCol, pch=16, main='PCA',xlab="PCA1",ylab="PCA2")
legend('bottomleft',c('P cells','SP cells'),pch=c(16,16),col=c(sgCol[1],sgCol[length(sgCol)]))  
#legend('bottomright',c('first seq','second seq'),pch=c(16,17),col=1)  
dev.off()


#v <- apply(mat, 1, var)
#vi <- names(sort(v)[1:1000])
#hc <- hclust(dist(t(mat[vi,])))
## visualize as heatmap
#heatmap(mat[vi,], Rowv=NA, Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=colorRampPalette(c('blue', 'white', 'red'))(100))

library("RColorBrewer")
co=colorRampPalette(brewer.pal(9, "YlGnBu"))(256)


png(file.path(resDir, "sampleclusterCorBoth.png"),width=1024,height=1024)
heatmap(mat, ColSideColors = sgCol,  col=co)
dev.off()

#########
#
# top variable genes..

vIds=order(apply(mat,1,var),decreasing=TRUE)[1:500]
heatmap(mat[vIds,], ColSideColors = sgCol,  col=co)

png(file.path(resDir, "sampleclusterCorTop500.png"),width=1024,height=1024)
heatmap(mat[vIds,], ColSideColors = sgCol,  col=co)
dev.off()



hc <- hclust(dist(t(mat)))
hr <- hclust(dist((mat)))

#heatmap(mat, Rowv=NA, Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=colorRampPalette(c('blue', 'white', 'red'))(100))
pdf(file.path(resDir, "sampleclusterCorRow.pdf"),width=6,height=6)
heatmap(mat, Rowv=NA, Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=co)
dev.off()

pdf(file.path(resDir, "sampleclusterCorBoth.pdf"),width=6,height=6)
heatmap(mat, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=co)
dev.off()


png(file.path(resDir, "sampleclusterCorRow.png"),width=512,height=512)
heatmap(mat, Rowv=NA, Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=co)
dev.off()

pdf(file.path(resDir, "sampleclusterCorNone.pdf"),width=6,height=6)
heatmap(mat, Rowv=NA, Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=co,scale="none")
dev.off()

hr <- hclust(dist((mat)))
heatmap(mat, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), ColSideColors = sgCol,  col=co, keep.dendro=FALSE)

heatmap.2(mat,dendrogram='none', Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),trace='none')


##############################
#
# gene list
#

load(file.path(serverDir,"projects/referenceSequences/rdata/human/gtf.rda"))

genenames=as.character(gtf$gene_name[match(rownames(nCountsHms),gtf$gene_id)])



#fix the problem where two genes have identical names
which(genenames=="TUBB3")
which(genenames=="RBFOX3")
genenames[18079]="TUBB3-2"
genenames[58392]="RBFOX3-2"


geneListFile=file.path(dataDir,"../2016-08-04/geneList.txt")

geneListR=read.delim(geneListFile,stringsAsFactors =F)

#addGenes=data.frame(gene=c("MT-RNR2","CDH8","PDCD6IP"),type="")

geneList=rbind(geneListR,addGenes)

geneList$gene[!geneList$gene %in% genenames]

mgIds=match(geneList$gene,genenames)

is_gl=genenames %in% geneList$gene

nch=nCountsHms
rownames(nch)=genenames
#matt=log2(as.matrix(nch)+1)[is_gl,]


#mgIds=match(c(geneList$gene,"CDH8","MT-RNR2"),genenames)
mgIds=match(c(geneList$gene),genenames)


sg=factor(infos$type[match(colnames(mat),infos$sample)])

matt <- log2(as.matrix(nch)+1)[mgIds,order(sg)]
#rownames(matt)=c(geneList$gene,"CDH8","MT-RNR2")
rownames(matt)=c(geneList$gene)

sgN=factor(infos$type[match(colnames(matt),infos$sample)])
sgColN <- rainbow(length(levels(sgN)))[sgN]


heatmap.2((as.matrix(matt)),dendrogram="none",col = co,ColSideColors = sgColN,trace="none")


#hc <- hclust(dist(t(matt)))
hr = hclust(dist(matt))
#heatmap(matt, Rowv=as.dendrogram(hr), Colv=NA, ColSideColors = sgCol,  col=colorRampPalette(c('blue', 'white', 'red'))(100))
#heatmap(matt, Rowv=as.dendrogram(hr), Colv=NA, ColSideColors = sgCol,  col=rev(topo.colors(75)))

#heatmap(matt, Rowv=as.dendrogram(hr), Colv=NA, ColSideColors = sgCol,  col=heat.colors(100))

#co=colorRampPalette(brewer.pal(9, "RdYlGn"))(256)


co=colorRampPalette(brewer.pal(9, "YlGnBu"))(256)
heatmap(matt, Rowv=as.dendrogram(hr), Colv=NA, ColSideColors = sgColN,  col=co)


pdf(file.path(resDir, "knowGenesClust.pdf"),width=12,height=12)
co=colorRampPalette(brewer.pal(9, "YlGnBu"))(256)
heatmap(matt, Rowv=as.dendrogram(hr), Colv=NA, ColSideColors = sgColN,  col=co)
dev.off()


pdf(file.path(resDir, "knowGenesClust_ht2.pdf"),width=10,height=10)
co=colorRampPalette(brewer.pal(9, "YlGnBu"))(256)
heatmap.2((as.matrix(matt)),col = co,ColSideColors = sgColN,trace="none", density.info="none")
dev.off()




heatmap.2((as.matrix(matt)),Colv=F,Rowv=F,col = co,ColSideColors = sgColN,trace="none")

heatmap.2((as.matrix(matt)),Colv="NA",col = co,ColSideColors = sgColN,trace="none")

#heatmap.2((as.matrix(matt)),Colv=F,Rowv=F,col = rev(topo.colors(75)),trace="none")
#heatmap(nch[mgIds,], Rowv=NA, Colv=NA, ColSideColors = sgColN,  col=colorRampPalette(c('blue', 'white', 'red'))(100))
#heatmap(log10(nch+1)[mgIds,], Rowv=NA, Colv=NA, ColSideColors = sgColN,  col=colorRampPalette(c('blue', 'white', 'red'))(100))









##########################
#
# most different genes between P and SP
#
#

meanp=apply(mat[,sg=="2242_P"],1,mean)
meansp=apply(mat[,sg=="2242_SP"],1,mean)

mdif=meanp-meansp
hist(mdif)

dIds= mdif>0.3 |mdif< -0.3

heatmap.2(mat[dIds,order(sg)],Colv="NA",col = co,ColSideColors = sgColN,trace="none")


############################
# correlation for all samples...
#
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- abs(cor(x, y,use="complete.obs"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex * r)
}

panel.smoothN=function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
		cex = 1, col.smooth = "red", span = 2/3, iter = 3,...) 
{
#	ids=which(x<1000 & y<1000)
	points(x, y, pch = pch, col = col, bg = bg, cex = cex)
	ok <- is.finite(x) & is.finite(y)
#	if (any(ok)) 
#		lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
}
#load(file.path(dataDir,"rdata/htcs.rda"))

png(paste(resDir, "ERCCcor.png",sep="/"),width=1024,height=1024)
pairs(log(countsERCC[,1:10]+1),lower.panel=panel.smoothN,upper.panel=panel.cor)
dev.off()

png(paste(resDir, "hgsCor.png",sep="/"),width=1024,height=1024)
pairs(log(countsHms[,1:10]+1),lower.panel=panel.smoothN,upper.panel=panel.cor)
dev.off()

#png(paste(resDir, "Neuroncor_p10_100.png",sep="/"),width=1024,height=1024)
#pairs(log(log10(countsHms+1)[,84:90]+1),lower.panel=panel.smooth,upper.panel=panel.cor)
#dev.off()


