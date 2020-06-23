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
    tidyr::pivot_longer(c(3,4,5),values_to='Value',names_to='Metric')
  
  drug.mut<-mutData%>%
    mutate(AD=as.numeric(AD))%>%
    tidyr::separate(specimenID,into=c('patient','sample'),remove=TRUE)%>%
    filter(sample=='xenograft')%>%
    dplyr::select(-c(patient,sample,tranche))%>%
    tidyr::pivot_wider(names_from=individualID,values_from=AD,values_fn=list(AD=mean),
                       values_fill=list(AD=0.0))%>%
    tidyr::pivot_longer(-Symbol,names_to="individualID",values_to="AD")%>%
    inner_join(drug.long)
  
  drugGeneCors<-drug.mut%>%
    group_by(Symbol,drug,Metric)%>%
    summarize(corVal=cor(AD,Value,method='spearman'))%>%right_join(drug.mut)
  
  return(drugGeneCors)
}