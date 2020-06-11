###biomarker stats
## collecting various evidence of biomarkers in drugs


#' drugTranscriptCor
#' Evaluates correlation 
#' @param drug.tab
#' @param metric
#' @param rna.tab
#' @export
drugTranscriptCor<-function(drug.tab,metric,rna.tab){
  full.res<-rna.tab%>%
    inner_join(drug.tab,by=c('individualID'))%>%
    group_by(drug,Symbol)%>%
    purrr::group_map(~cor(totalCounts,!!metric))
  
    dplyr::summarize(drugGene=cor(totalCounts,!!metric))
}