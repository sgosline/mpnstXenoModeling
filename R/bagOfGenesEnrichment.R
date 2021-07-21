
#' compute gene set enrichment - osama's code wrapped in package.
#' @export
#' @import ggplot2
#' @author Osama 
#' @param genes.with.values of genes and difference values
#' @param prot.univ the space of all proteins we are considering
#' @return gSEA output type stuff
#' computeGSEA<-function(genes.with.values,prefix,gsea_FDR=0.01){
#'   
#'   library(WebGestaltR)
#'   library(ggplot2)
#'   inputdfforWebGestaltR <- genes.with.values%>%
#'     dplyr::rename(genes='Gene',scores='value')%>%
#'     dplyr::arrange(scores)
#'   
#'   
#'   #' * GSEA using gene ontology biological process gene sets
#'   
#'   go.bp.res.WebGestaltR <- WebGestaltR(enrichMethod = "GSEA", 
#'                                        organism="hsapiens", 
#'                                        enrichDatabase="geneontology_Biological_Process", 
#'                                        interestGene=inputdfforWebGestaltR, 
#'                                        interestGeneType="genesymbol", 
#'                                        collapseMethod="mean", perNum = 1000,
#'                                        fdrThr = gsea_FDR, nThreads = 2, isOutput = F)
#'   write.table(go.bp.res.WebGestaltR, paste0("proteomics_", prefix, "_gseaGO_result.txt"), sep="\t", row.names=FALSE, quote = F)
#'   
#'   top_gseaGO <- go.bp.res.WebGestaltR %>% 
#'     filter(FDR < gsea_FDR) %>% 
#'     dplyr::rename(pathway = description, NES = normalizedEnrichmentScore) %>% 
#'     arrange(desc(NES)) %>% 
#'     dplyr::mutate(status = case_when(NES > 0 ~ "Up",
#'                                      NES < 0 ~ "Down"),
#'                   status = factor(status, levels = c("Up", "Down"))) %>% 
#'     #\group_by(status) %>% 
#'     top_n(30, wt = NES) %>% 
#'     ungroup() %>% 
#'     ggplot2::ggplot(aes(x=reorder(pathway, NES), y=NES)) +
#'     geom_bar(stat='identity', aes(fill=status)) +
#'     scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
#'     coord_flip() +
#'     theme_minimal() +
#'     theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
#'           axis.title.x = element_text(size=16),
#'           axis.title.y = element_blank(), 
#'           axis.text.x = element_text(size = 14),
#'           axis.text.y=element_text(size = 14),
#'           axis.line.y = element_blank(),
#'           axis.ticks.y = element_blank(),
#'           legend.position = "none") +
#'     labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
#'     ggtitle(paste('Up-regulated',prefix))
#'   ggsave(paste0("upRegProts_", prefix,"_gseaGO_plot.pdf"), top_gseaGO, height = 8.5, width = 11, units = "in")
#'   
#'   
#'   all_gseaGO <- go.bp.res.WebGestaltR %>% 
#'     filter(FDR < gsea_FDR) %>% 
#'     dplyr::rename(pathway = description, NES = normalizedEnrichmentScore) %>% 
#'     arrange(NES) %>% 
#'     dplyr::mutate(status = case_when(NES > 0 ~ "Up",
#'                                      NES < 0 ~ "Down"),
#'                   status = factor(status, levels = c("Up", "Down"))) %>% 
#'     group_by(status) %>% 
#'     top_n(20, wt = abs(NES)) %>% 
#'     ungroup() %>% 
#'     ggplot2::ggplot(aes(x=reorder(pathway, NES), y=NES)) +
#'     geom_bar(stat='identity', aes(fill=status)) +
#'     scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
#'     coord_flip() +
#'     theme_minimal() +
#'     theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
#'           axis.title.x = element_text(size=16),
#'           axis.title.y = element_blank(), 
#'           axis.text.x = element_text(size = 14),
#'           axis.text.y=element_text(size = 14),
#'           axis.line.y = element_blank(),
#'           axis.ticks.y = element_blank(),
#'           legend.position = "none") +
#'     labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
#'     ggtitle(paste('All',prefix))
#'   ggsave(paste0("allRegProts_", prefix,"_gseaGO_plot.pdf"), all_gseaGO, height = 8.5, width = 11, units = "in")
#'   
#'   
#'   bot_gseaGO <- go.bp.res.WebGestaltR %>% 
#'     filter(FDR < gsea_FDR) %>% 
#'     dplyr::rename(pathway = description, NES = normalizedEnrichmentScore) %>% 
#'     arrange(NES) %>% 
#'     dplyr::mutate(status = case_when(NES > 0 ~ "Up",
#'                                      NES < 0 ~ "Down"),
#'                   status = factor(status, levels = c("Up", "Down"))) %>% 
#'     #group_by(status) %>% 
#'     top_n(40, wt = rev(NES)) %>% 
#'     ungroup() %>% 
#'     ggplot2::ggplot(aes(x=reorder(pathway, rev(NES)), y=NES)) +
#'     geom_bar(stat='identity', aes(fill=status)) +
#'     scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
#'     coord_flip() +
#'     theme_minimal() +
#'     theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
#'           axis.title.x = element_text(size=16),
#'           axis.title.y = element_blank(), 
#'           axis.text.x = element_text(size = 14),
#'           axis.text.y=element_text(size = 14),
#'           axis.line.y = element_blank(),
#'           axis.ticks.y = element_blank(),
#'           legend.position = "none") +
#'     labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
#'     ggtitle(paste('Down-regulated',prefix))
#'   ggsave(paste0("downRegProts_", prefix,"_gseaGO_plot.pdf"), bot_gseaGO, height = 8.5, width = 11, units = "in")
#'   
#'   return(go.bp.res.WebGestaltR) 
#' }



#' Old plot using clusterProfiler
#' @export 
#' @import BiocManager
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
  #gr<-clusterProfiler::gseKEGG(genelist[!is.na(genelist)],organism='hsa',keyType="kegg",
  #OrgDb=org.Hs.eg.db,
  #                           pAdjustMethod = 'BH')#,eps=1e-10)
  
  # if(nrow(as.data.frame(gr))==0){
  #    gr<-clusterProfiler::gseGO(genelist[!is.na(genelist)],ont="BP",keyType="SYMBOL",
  #                             OrgDb=org.Hs.eg.db,pAdjustMethod = 'BH',pvalueCutoff = 0.1)#,eps=1e-10)
  #  }
  
  enrichplot::ridgeplot(gr,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("KEGG Terms for ",prefix))
  ggplot2::ggsave(paste0(prefix,'_KEGG.pdf'),width=10,height=10)
  
  
  return(gr)
}

#'Runs regular bag of
#'@name doRegularGo
#'@description Performs GO enrichment
#'@export 
#'@import BiocManager
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
#'uses Osama's code to compute de from limma
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
