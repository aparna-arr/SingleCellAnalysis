library(ggplot2)
library(DESeq2)
library(genefilter)
library("Biostrings")
#library("ShortRead")
library(statmod)
library(gplots)


## read in a count file. Genes as the columns, cells as the rows
#geneCountFile <- read.table("/Users/kristen/Documents/Stanford/Steinmetz_lab/mTECs/Analysis/RAnalysis/HTSeq_gene_file.txt",
#                              sep="\t", header=TRUE, row.names = 1)
geneCountFile <- read.table("countfile.tsv", sep="\t", header=TRUE, row.names = 1)
resDir <- "./"
#resDir <- "/Users/kristen/Documents/Stanford/Steinmetz_lab/mTECs/Analysis/RAnalysis/"

#load("/Users/kristen/Documents/Stanford/Steinmetz_lab/mTECs/Analysis/RAnalysis/sampleInfo.rda")

# establish if each gene is ERCC or mouse
geneTypes <- factor(c(ENSM="ENSM", ERCC="ERCC") [
  substr(rownames(geneCountFile), 1, 4)])

# pull out mouse genes
mouseGenes <- geneCountFile[which(geneTypes== "ENSM"),]

# determines which cells have fewer than 1000 genes with counts 20 or higher
ii=which(colSums(mouseGenes>=20)<1000)
colSums(mouseGenes>=20)[ii]

# boolean of which cells have greater than 1000 genes with counts 20 or higher
single_cells <- colSums(mouseGenes>=20) >1000

# only mouse genes and cells that pass the treshold
countsMouse <- geneCountFile[ which( geneTypes=="ENSM" ), single_cells]

# only ERCCs and cells that pass the threshold
countsERCC <- geneCountFile[ which( geneTypes=="ERCC" ), single_cells]
lengthsMouse <- geneCountFile[ which( geneTypes=="ENSM" ), 1 ]
lengthsERCC <- geneCountFile[ which( geneTypes=="ERCC" ), 1 ]

# estimateSizeFactors from DEseq
sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMouse <- estimateSizeFactorsForMatrix( countsMouse )

# normalized counts
nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMouse <- t( t(countsMouse) / sfMouse )


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
    plot( meansEndo, cv2Endo, log="xy", col=1, xlab = 'Means of mouse genes', ylab = 'CV2')
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


techNoise = fitTechnicalNoise(nCountsMouse,nCountsERCC=nCountsERCC, fit_type = 'counts')  
#is_het = getVariableGenes(nCountsHms, techNoise, method = "fdr", threshold = 0.01, fit_type="counts",sfEndo=sfHms, sfERCC=sfERCC)
is_het = getVariableGenes(nCountsMouse, techNoise, method = "fdr", threshold = 1e-3, fit_type="counts",sfEndo=sfMouse, sfERCC=sfERCC)

table(is_het)

mat <- log10(as.matrix(nCountsMouse)+1)[is_het,]



breaks()

base.pca <- prcomp(t(mat))

