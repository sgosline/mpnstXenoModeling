


#' Old plot using clusterProfiler
#' @export 
#' @require org.Hs.eg.db
#' @import clusterProfiler
plotOldGSEA<-function(genes.with.values,prot.univ,prefix){
  require(org.Hs.eg.db)
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
#'@export 
#'@require org.Hs.eg.db
#'@import clusterProfiler
doRegularGo<-function(genes,bg=NULL){
  require(org.Hs.eg.db)
  #genes<-unique(as.character(genes.df$Gene))
  mapping<-as.data.frame(org.Hs.egALIAS2EG)%>%
    dplyr::rename(Gene='alias_symbol')
  
  eg<-subset(mapping,Gene%in%genes)
  
  res<-clusterProfiler::enrichGO(eg$gene_id,'org.Hs.eg.db',keyType='ENTREZID',ont='BP')
    #sprint(res)
  ret=as.data.frame(res)%>%
    dplyr::select(ID,Description,pvalue,p.adjust)
  return(ret)
  
  
}
