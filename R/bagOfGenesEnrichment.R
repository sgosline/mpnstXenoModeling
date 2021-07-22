


#' Old plot using clusterProfiler
#' @name plotOldGSEA
#' @param genes.with.values data frame of gene names and values
#' @param prot.univ total proteins
#' @param prefix used to create file
#' @export 
#' @import BiocManager
#' 
plotOldGSEA<-function(genes.with.values,prot.univ,prefix){
  
  if(!require(org.Hs.eg.db)){
    BiocManager::install('Biobase')
    require(org.Hs.eg.db)
  }
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  genes.with.values<-genes.with.values%>%
   # dplyr::left_join(mapping,by='Gene')%>%
    arrange(desc(value))
  
  genelist=genes.with.values$value
  names(genelist)=genes.with.values$Gene
  
  # symbs<-names(genelist)[!is.na(genelist)]
  # xx <- as.list(org.Hs.egALIAS2EG)
  # ents<-unlist(sapply(intersect(names(xx),symbs), function(x) xx[[x]]))
  # print(ents)
  if(!require(clusterProfiler)){
    BiocManager::install('clusterProfiler')
    library(clusterProfiler)
  }
  gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH')#,eps=1e-10)
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("KEGG Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_KEGG.pdf'),width=10,height=10)
  
  return(gr)
}

#' Runs regular bag of genes enrichment
#' @name doRegularGo
#' @description Performs GO enrichment
#' @export 
#' @import BiocManager
#' 
doRegularGo<-function(genes,bg=NULL){
  if(!require(org.Hs.eg.db)){
    BiocManager::install('Biobase')
    require(org.Hs.eg.db)
  }
  #genes<-unique(as.character(genes.df$Gene))
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  eg<-subset(mapping,Gene%in%genes)
  if(!require(clusterProfiler)){
    BiocManager::install('clusterProfiler')
    require(clusterProfiler)
  }
  res<-clusterProfiler::enrichGO(eg$gene_id,'org.Hs.eg.db',keyType='ENTREZID',ont='BP')
    #sprint(res)
  ret=as.data.frame(res)%>%
    dplyr::select(ID,Description,pvalue,p.adjust)
  return(ret)
  
  
}

#'
#'limmaTwoFactorDEAnalysis
#'@name limmaTwoFactorDEAnalysis
#'@description Runs limma on two groups
#'@author Osama
#'@import BiocManager
#'@export
#'@param data matrix
#'@param group1 ids
#'@param group2 ids
limmaTwoFactorDEAnalysis <- function(dat, sampleIDs.group1, sampleIDs.group2) {
  # Conduct DE expression analysis using limma from the expression matrix dat (group2 vs group1, group1 is reference)
  #
  # Args:
  #   dat: Expression data matrix, rows are genes, columns are samples
  #   sampleIDs.group1: Vector with ids of samples in reference group (eg. normal samples)
  #   sampleIDs.group2: Vector with ids of samples in interest group (eg. tumor samples) 
  #
  # Returns:
  #   limma Differential Expression results.
  #
  #http://www.biostat.jhsph.edu/~kkammers/software/CVproteomics/R_guide.html
  #http://genomicsclass.github.io/book/pages/using_limma.html
  #https://wiki.bits.vib.be/index.php/Tutorial:_Testing_for_differential_expression_I
  if(!require('limma')){
    BiocManager::install('limma')
    library(limma)
  }
  fac <- factor(rep(c(2,1), c(length(sampleIDs.group2), length(sampleIDs.group1))))
  design <- model.matrix(~fac)
  fit <- lmFit(dat[,c(sampleIDs.group2, sampleIDs.group1)], design)
  fit <- eBayes(fit)
  print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="none")
  res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  return(arrange(res,P.Value))
}
