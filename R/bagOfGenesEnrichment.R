


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

#'plotTopGenesHeatmap
#'@name plotTopGenesHeatmap
#'@description Filters and plots expression matrix
#'@author Jess
#'@import BiocManager
#'@import reticulate
#'@export
#'@param data matrix
#'@param txid, genename identifieres
#'@param str var
#'@param df of var
plotTopGenesHeatmap <- function(de.out, counts, identifiers, myvar, var.ID, adjpval=0.5, upload=FALSE, path='.', parentID=NULL) {
  # Downfilter DE expression table by Adjusted P Value and generate pheatmap
  #
  # Args:
  #   de.out: Expression data matrix, rows are genes, columns are samples
  #   counts: Count data matrix, columns are sampleID.Count or sampleID.TPM, rows are TXID
  #   identifiers: Dataframe, columns include 'TXID', 'GENENAME', merged output from do_ensembl_match
  #   myvar: Str, variable tested for differential expression
  #   var.ID: Dataframe, rows are sample IDs and columns are variables used in DE
  #   upload: Bool, should filtered gene dataframe be written to file and uploaded to synapse
  #   path: Str, path to filewrite location
  # Returns:
  #   pheatmap results.
  library(reticulate)
  if(!require('pheatmap')){
    BiocManager::install('pheatmap')
    library(pheatmap)
  }
  if(!require('edgeR')){
    BiocManager::install('edgeR')
    library(edgeR)
  }
  if(!require('tibble')){
    BiocManager::install('tibble')
    library(tibble)
  }
  synapse=reticulate::import('synapseclient')
  sync=synapse$login()

  names(de.out)[names(de.out) == "featureID"] <- "TXID"
  dge <- DGEList(counts)
  dge <- calcNormFactors(dge)

  #gather read counts that are normalized to 1 via calcnormfactors()
  norm.counts <- dge$counts
  norm.counts <- as.data.frame(norm.counts)
  norm.counts <- tibble::rownames_to_column(norm.counts,var='TXID')

  #combined differentially expressed txids, genenames, and normalized counts
  de.df <- left_join(de.out,identifiers,by="TXID")
  de.df <- left_join(de.df,norm.counts,by="TXID")
  de.df <- de.df[de.df$adj.P.Val < adjpval,]
  if (isTRUE(upload)) {
    write.csv(de.df, file.path(path, paste0(myvar,'_topgenes_adjpval_',adjpval,'.csv')))
    synapseStore(file.path(path,paste0(myvar,'_topgenes_adjpval_',adjpval,'.csv')),parentId=parentID)
  }
  if (dim(de.df)[1] == 0) {
    stop("No top genes within specified adj.p.val threshold to make heatmap")
  }

  de.table <- de.df %>% select(GENENAME,contains(c('WU', 'JHU', 'MN')))
  de.table <- de.table[!is.na(de.table$GENENAME),]
  de.table <- de.table[!duplicated(de.table$GENENAME), ]
  de.table <- tibble::remove_rownames(de.table)
  de.table <- tibble::column_to_rownames(de.table,var = 'GENENAME')

  h.t <- sapply(de.table,as.numeric)
  h.t <- as.matrix(h.t)
  rownames(h.t) = rownames(de.table)
  # remove genes with no counts
  h.t <- h.t[rowSums(h.t)>0,]

  # Convert counts to zScore of relevant genes
  cal_z_score <- function(x){(x - mean(x)) / sd(x)}
  data_subset_norm <- t(apply(h.t, 1, cal_z_score))

  # Pairwise correlation between samples (columns)
  cols.cor <- cor(h.t, use = "pairwise.complete.obs", method = "pearson")
  # Pairwise correlation between rows (genes)
  rows.cor <- cor(t(h.t), use = "pairwise.complete.obs", method = "pearson")

  # Plot pheatmap
  heatmap <- pheatmap(data_subset_norm,
                     show_rownames=TRUE,
                     cellheight=10,
                     annotation_col=var.ID,
                     clustering_distance_cols = as.dist(1 - cols.cor),
                     clustering_distance_rows = as.dist(1 - rows.cor),
                     filename=file.path(path, paste0(myvar,'_DE_heatmap_adjpval',adjpval,'.png'))
                    )
  if (isTRUE(upload)) {
    synapseStore(file.path(path, paste0(myvar,'_DE_heatmap_adjpval',adjpval,'.png')),parentId=parentID)
  }
  return(heatmap)

}
