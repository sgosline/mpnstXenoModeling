
#' Used to make reversed logarithmic scales
#' @import scales
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}


#'ploGenesetResults
#'Internal function to plot results of cluster profile
#'@param res
#'@param prefix
#'@param pathway.plot.size
#'@param order.by
#'@param clean.names
#'@param width
#'@param height
plotGenesetResults<-function(res,prefix,pathway.plot.size=3,
                             order.by='NES',clean.names=F,width=11,height=8.5){
  
  all.gsea<-res %>% 
    dplyr::rename(pathway = 'Description') %>% 
    dplyr::rename(tosort = order.by) %>% 
    dplyr::mutate(status = case_when(tosort > 0 ~ "Up", tosort < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    group_by(status) %>% 
    top_n(20, wt = abs(tosort) )%>% 
    ungroup()
  
  p.NES <- ggplot(all.gsea, aes(x = tosort, y = reorder(pathway, tosort))) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c(Up = "firebrick2", Down = "dodgerblue3")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 11),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") + 
    labs(x = order.by) +
    ggtitle("Normalized Enrichment Score")
  
  p.Pval <- ggplot(all.gsea, aes(x = p.adjust, y = reorder(pathway, tosort))) +
    scale_x_continuous(trans = reverselog_trans(10)) +  
    theme_minimal() +
    geom_bar(stat = "identity") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18), 
          axis.title.x = element_text(size = 16), 
          axis.text.x = element_text(size = 12), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_blank(), 
          legend.position = "none") + 
    labs(x = "Adjusted p-value") + 
    ggtitle("Significance")
  
  arrange_matrix <- t(as.matrix(c(rep(1,pathway.plot.size),2)))
  p.both <- grid.arrange(p.NES,p.Pval, layout_matrix = arrange_matrix)
  
  try(ggsave(paste0("sig-included", prefix,"-gsea-plot.png"), p.both, 
             height = height, width = width, units = "in"))
  
  #df<-as.data.frame(res)%>%mutate(Condition=prefix)
  
  return(p.both)
}

#' Old plot using clusterProfiler
#' @name doGSEA
#' @param genes.with.values data frame of gene names and values
#' @param prot.univ total proteins
#' @param prefix used to create file
#' @export
#' @import BiocManager
#'
doGSEA<-function(genes.with.values,prot.univ,prefix,useEns=FALSE,pathway.plot.size=3,
                      order.by='NES',clean.names=F,width=11,height=8.5,  gsea_FDR=0.05){

  if(!require('org.Hs.eg.db')){
    BiocManager::install('org.Hs.eg.db')
    library(org.Hs.eg.db)
  }
  
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')

  emapping <-as.data.frame(org.Hs.egENSEMBLTRANS2EG)%>%
    dplyr::rename(Gene='trans_id')
  
  if(useEns)
    mapping <- emapping
  
  genes.with.values<-genes.with.values%>%
    dplyr::left_join(mapping,by='Gene')%>%
    arrange(desc(value))%>%
    subset(!is.na(gene_id))

  genelist=genes.with.values$value
  names(genelist)=as.character(genes.with.values$gene_id)

 # print(head(genelist))
  # symbs<-names(genelist)[!is.na(genelist)]
  # xx <- as.list(org.Hs.egALIAS2EG)
  # ents<-unlist(sapply(intersect(names(xx),symbs), function(x) xx[[x]]))
  # print(ents)
  if(!require(clusterProfiler)){
    BiocManager::install('clusterProfiler')
    library(clusterProfiler)
  }
  gr<-NULL
  try(gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="ENTREZID",
                               OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH'))
  res<-filter(as.data.frame(gr),p.adjust<gsea_FDR)
  
  if(nrow(res)==0){
    return(gr)
  }else{
    return(plotGenesetResults(res,prefix=prefix,pathway.plot.size=pathway.plot.size,
                                        order.by='NES',clean.names=F,width=width,height=height))
  }
  
  
}

#' Runs regular bag of genes enrichment
#' @name doRegularGo
#' @description Performs GO enrichment
#' @export
#' @import BiocManager
#'
doRegularGo<-function(genes,bg=NULL,prefix='',gsea_FDR=0.05,pathway.plot.size=3,width=11,height=8.5){
  if(!require(org.Hs.eg.db)){
    BiocManager::install('org.Hs.eg.db')
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
#  print(head(eg))
  res<-clusterProfiler::enrichGO(eg$gene_id,'org.Hs.eg.db',keyType='ENTREZID',ont='BP',qvalueCutoff =0.2)
    #sprint(res)
  ret=as.data.frame(res)%>%
    dplyr::select(ID,Description,pvalue,p.adjust,Count)
  
  res<-filter(as.data.frame(res),p.adjust<gsea_FDR)
  
  if(nrow(res)==0){
    return(res)
  }else{
    return(plotGenesetResults(res,prefix=prefix,pathway.plot.size=pathway.plot.size,
                              order.by='Count',clean.names=F,width=width,height=height))
  }
  return(res)


}



#' ds2FactorDE
#' @name ds2FactorDE
#' @author Sara
#' @import BiocManager
#' @param dds DESeq object
#' @param ids1 Sample ids
#' @param ids2 Other sample ids
#' @param name for condition
#' @export
#' @import BiocManager
ds2FactorDE<-function(dds,ids1,ids2,name){
  if(!require('DESeq2')){
    BiocManager::install('DESeq2')
    library(DESeq2)
  }
  ##create an additional column

  tcd<-colData(dds)

  ##chek namess
  t.ids1<-intersect(ids1,rownames(tcd))
  t.ids2<-intersect(ids2,rownames(tcd))
  
  message(paste("Found",length(t.ids1),'samples that overlap with expression out of',length(ids1)))
  message(paste("Found",length(t.ids2),'samples that overlap with expression out of',length(ids2)))
  
  tcd$newvar<-rep(NA,nrow(tcd))
  tcd[t.ids1,]$newvar<-TRUE
  tcd[t.ids2,]$newvar<-FALSE
  tcd<-subset(tcd,!is.na(newvar))#%>%
    #dplyr::rename(newvar=name)

  new.counts<-counts(dds)[,rownames(tcd)]
 
  new.dds<-DESeq2::DESeqDataSetFromMatrix(new.counts,design=~newvar,
                                          colData=tcd)
  ###re run dds  
  #design(dds)<-~newvar
  new.dds <- DESeq(new.dds) ##rerun
  res <- results(new.dds)#,contrasts=c("newvar","TRUE","FALSE"))
#  print(summary(res))  
  as.data.frame(results(new.dds))%>%arrange(pvalue)

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
  sampleIDs.group1 <- intersect(sampleIDs.group1,colnames(dat))
  sampleIDs.group2 <- intersect(sampleIDs.group2,colnames(dat))
  
  fac <- factor(rep(c(2,1), c(length(sampleIDs.group2), length(sampleIDs.group1))))
  design <- model.matrix(~fac)
  fit <- lmFit(dat[,c(sampleIDs.group2, sampleIDs.group1)], design)
  fit <- eBayes(fit)
#  print(topTable(fit, coef=2))
  res <- topTable(fit, coef=2, number=Inf, sort.by="none")
  res <- data.frame(featureID=rownames(res), res, stringsAsFactors = F)
  return(arrange(res,P.Value))
}


#' geneIdToSymbolMatrix
#' This takes a symbol wtih gene ids and maps them to symbols
#' @param gene.mat
#' @param identifiers data table of identifiers
#' @export
geneIdToSymbolMatrix<-function(gene.mat,identifiers){
  count.mat<-gene.mat%>%
    as.data.frame()%>%
    tibble::rownames_to_column("GENEID")%>%
    tidyr::pivot_longer(cols=c(-GENEID),names_to='Sample',values_to='counts')%>%
    left_join(tibble::rownames_to_column(identifiers,'GENEID'))%>%
    group_by(GENENAME)%>%
    dplyr::select(Sample,counts,GENENAME)%>%distinct()%>%
    subset(GENENAME!="")%>%
    tidyr::pivot_wider(names_from=Sample,values_from=counts,values_fn=list(counts=mean))%>%
    tibble::column_to_rownames('GENENAME')%>%as.matrix()
  
  return(count.mat)
}

#'plotTopGenesHeatmap
#'@name plotTopGenesHeatmap
#'@description Filters and plots expression matrix
#'@author Jess
#'@import BiocManager
#'@import reticulate
#'@import pheatmap
#'@export
#'@param de.out diffex results
#'@param dds, DESEq object
#'@param identifiers mapping to gene name
#'@param myvar is name of variable
plotTopGenesHeatmap <- function(de.out, dds, identifiers, myvar, patients=NULL, adjpval=0.5, 
                                upload=FALSE, path='.', parentID=NULL,newVar="") {
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
  if(!require('tibble')){
    BiocManager::install('tibble')
    library(tibble)
  }
  synapse=reticulate::import('synapseclient')
  sync=synapse$login()

  de.out<-cbind(de.out,identifiers[rownames(de.out),])%>%
    subset(!is.na('GENENAME'))%>%
    subset(GENENAME!="")
  
  if(is.null(patients))
    patients <- rownames(colData(dds))
  else
    patients <- intersect(patients,rownames(colData(dds)))
  
  print(paste0('plotting expression across ',length(patients),' samples'))
  if(length(patients)==0)
    return(NULL)
  #  names(de.out)[names(de.out) == "featureID"] <- "TXID"
  
  #combined differentially expressed txids, genenames, and normalized counts
 # de.df <- left_join(de.out,identifiers,by="TXID")
#  de.df <- left_join(de.df,norm.counts,by="TXID")
#  de.df <- de.df[de.df$adj.P.Val < adjpval,]
  if (isTRUE(upload)) {
    write.csv(de.out, file.path(path, paste0(myvar,'_topgenes_adjpval_',adjpval,'.csv')))
    synapseStore(file.path(path,paste0(myvar,'_topgenes_adjpval_',adjpval,'.csv')),parentId=parentID)
  }
  if (dim(de.out)[1] == 0) {
    print("No top genes within specified adj.p.val threshold to make heatmap")
    return(NULL)
  }

  
  sigs <-subset(de.out,padj<adjpval)%>%dplyr::select(GENENAME)
  
  if(newVar!="")
    all.vars <- c('Sex','MicroTissueQuality','Clinical Status','Age',newVar)
  else
    all.vars <-c('Sex','MicroTissueQuality','Clinical Status','Age')
  
  var.ID<-colData(dds)[,all.vars]%>%
    as.data.frame()%>%
    mutate(MicroTissueQuality=unlist(MicroTissueQuality))

  count.mat <- geneIdToSymbolMatrix(counts(dds,normalized=TRUE)[rownames(sigs),],identifiers)
 
  count.mat<-count.mat[,patients]
  var.ID <- var.ID[patients,]
  library(pheatmap)
    heatmap <- pheatmap(log10(0.01+count.mat),
                     cellheight=10,
                     annotation_col=var.ID,
                     filename=file.path(path, paste0(myvar,'_DE_heatmap_adjpval',adjpval,'.png')))

  if (isTRUE(upload)) {
    synapseStore(file.path(path, paste0(myvar,'_DE_heatmap_adjpval',adjpval,'.png')),parentId=parentID)
  }
  return(heatmap)

}

#' Plot using correlation enrichment from leapR package. 
#' A single plot is saved to the working directory
#' @export 
#' @import ggplot2
#' @import gridExtra
#' @import scales
#' @import dplyr
#' @import remotes
#' @param exprs A matrix of intensities with accessions as row names, along with samples in the columns.
#' @param prefix string, used for naming the saved plots.
#' @param order.by This determines how the pathways are sorted. Default is pathway correlation of "Ingroup mean", but can also use "BH_pvalue" to sort by significance of the pathways.
#' @param geneset Pathway/Kinase database, eg ncipid, msigdb, both of which are included in leapr.
#' @param clean.names Boolean, if TRUE removes the "_pathway" ending in pathway names, making the plot easier to read.
plotCorrelationEnrichment <- function(exprs, geneset, fdr.cutoff = 0.05, 
                                      corr.cutoff = 0.1, prefix, width = 11, 
                                      height = 8.5, order.by = "Ingroup mean", 
                                      clean.names = FALSE, pathway.plot.size = 3) {
  if(!require('leapr')){
    remotes::install_github('biodataganache/leapr')
    library(leapr)
  }
  corr.enrichment <- leapr::leapR(geneset, 
                           enrichment_method = "correlation_enrichment",
                           datamatrix = exprs) 
  corr.enrichment <- corr.enrichment %>%
    mutate(Pathway = rownames(.)) %>%
    rename(`Ingroup mean` = ingroup_mean,
           `Outgroup mean` = outgroup_mean) %>%
    mutate(Status = case_when(`Ingroup mean` > 0 ~ "Positively Correlated", 
                              `Ingroup mean` < 0 ~ "Negatively Correlated")) %>%
    select(Pathway, `Ingroup mean`, `Outgroup mean`, 
           ingroup_n, outgroup_n, pvalue, BH_pvalue, Status)
  
  corr.enrichment.filtered <- corr.enrichment %>%
    filter(BH_pvalue < fdr.cutoff & abs(`Ingroup mean`) > corr.cutoff) %>%
    mutate(BH_pvalue = case_when(BH_pvalue > 1e-10 ~ BH_pvalue,
                                 BH_pvalue < 1e-10 ~ 1e-10))
  
  if (clean.names) {
    corr.enrichment.filtered$Pathway <- sub("_pathway$", "", 
                                            corr.enrichment.filtered$Pathway)
  }
  
  #sprint(res)
 # ret=as.data.frame(res)%>%
#    dplyr::select(ID,Description,pvalue,p.adjust,Count)
  
#  res<-filter(as.data.frame(res),p.adjust<gsea_FDR)
  
#  if(nrow(res)==0){
#    return(res)
#  }else{
#    return(plotGenesetResults(res,prefix=prefix,pathway.plot.size=pathway.plot.size,
#                              order.by='Count',clean.names=F,width=width,height=height))
#  }
  
  p.corr <- ggplot(corr.enrichment.filtered, aes(x = `Ingroup mean`, 
                                                 y = reorder(Pathway, `Ingroup mean`))) +
    geom_bar(stat='identity', aes(fill = Status)) +
    scale_fill_manual(values = c("Positively Correlated" = "mediumturquoise", 
                                 "Negatively Correlated" = "firebrick2")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 9),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") + 
    labs(x = "Average correlation") +
    ggtitle("Correlation Enrichment")
  
  p.pval <- ggplot(corr.enrichment.filtered, aes(x = BH_pvalue, 
                                                 y = reorder(Pathway, `Ingroup mean`))) +
    geom_bar(stat='identity') +
    scale_x_continuous(trans = reverselog_trans(10)) + 
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_line(color = "black"),
          legend.position = "none") + 
    labs(x = "Adjusted p-value") +
    ggtitle("Significance")
  
  arrange_matrix <- t(as.matrix(c(rep(1,pathway.plot.size),2)))
  p.both <- grid.arrange(p.corr,p.pval,layout_matrix = arrange_matrix)
  
  try(ggsave(paste0("sig-included-", prefix,"-correlation-enrichment-plot.png"), p.both, 
         height = height, width = width, units = "in"))
  
  return(p.both)
}
