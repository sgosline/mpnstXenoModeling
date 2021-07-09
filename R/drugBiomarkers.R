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


#' drugMutationCor
#' Evaluates correlation betwen allele depth and various mutational metrics
#' @param pat.drug output of `getAllDrugStats`
#' @param mutData output of ...
#' @import dplyr
#' @import tidyr
#' @return table of drug, gene, metric, and correlation value
#' @export
drugMutationsCor<-function(pat.drug,mutData){
  #reshape drug data
  drug.long<-pat.drug%>%
    mutate(drug=tolower(drug))%>%
    tidyr::pivot_longer(c(3,4,5),values_to='Value',names_to='Metric')
  
  drug.mut<-mutData%>%
      inner_join(drug.long)
  
  drugGeneCors<-drug.mut%>%
    group_by(Symbol,drug,Metric)%>%
    summarize(corVal=cor(AD,Value,method='pearson',use='pairwise.complete.obs'),.groups='keep')%>%
    right_join(drug.mut)%>%
    subset(!is.na(corVal))
  
  return(drugGeneCors)
}