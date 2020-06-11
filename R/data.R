##get datasets from figshare and whatever drug screening data we have. 

#' loadPDXData gets data from synapse and stores them as global variables
#' @export
#' @import reticulate
#' @import dplyr
#' @import tidyr
loadPDXData<-function(){
  library(reticulate)
  library(dplyr)
  syn<-reticulate::import('synapseclient')$login()
  
  newMutData<-getLatestVariantData(syn)%>%select(Symbol='Gene',AD,specimenID,individualID)
  mutData<-getOldVariantData(syn)%>%mutate(AD=as.numeric(AD))
  
  mutData<<-rbind(data.frame(mutData,tranche='oldData'),data.frame(newMutData,tranche='newData'))
  
  
  csvs = syn$tableQuery(paste0("SELECT id,individualID FROM syn11678418 WHERE \"dataType\" = 'drugScreen'"))$asDataFrame()
  ids<-csvs$id
  indiv<-csvs$individualID
  names(indiv)<-ids
  res=do.call(rbind,lapply(names(indiv),function(x)
  { 
    read.csv(syn$get(x)$path)%>%
      dplyr::select(model.id='individual_id',specimen_id,drug='compound_name',volume='assay_value',time='experimental_time_point')%>%
      dplyr::mutate(individualID=indiv[[x]])
  }))
  res$drug <-sapply(res$drug,function(x) gsub('Doxorubinsin','doxorubicin',gsub('N/A','vehicle',x)))
  nas=which(res$model.id=="")
  res$individualID<-sapply(res$individualID,function(x) gsub('2-','JHU',x))
  if(length(nas)>0)
    res<-res[-nas,]
  drugData<<-res
    
  #now get RNA-Seq
  rnaSeq<<-getPdxRNAseqData(syn)%>%
    dplyr::select(totalCounts,Symbol,zScore,specimenID,individualID,sex,species,experimentalCondition)
  
 
}



#'getPdxRNAseqData gets all rna seq counts for xenografts
#'#'@export
getPdxRNAseqData<-function(syn){
  wu.rnaSeq=syn$tableQuery("SELECT * FROM syn21054125 where transplantationType='xenograft'")$asDataFrame()
  jh.rnaSeq=syn$tableQuery("SELECT * FROM syn20812185 where transplantationType='xenograft'")$asDataFrame()

  com.cols=intersect(colnames(wu.rnaSeq),colnames(jh.rnaSeq))%>%setdiff(c("ROW_ID","ROW_VERSION"))
  count.tab=rbind(wu.rnaSeq[,com.cols],jh.rnaSeq[,com.cols])
  count.tab$individualID<-sapply(count.tab$individualID,function(x) gsub('2-','JHU',x))
  count.tab$specimenID<-sapply(count.tab$specimenID,function(x) gsub('2-','JHU',x))

  
  return(count.tab)
  
}

#' getAllNF1Expression collects all data from NF1 processed data in synapse
getAllNF1Expression<-function(syn){
  tabs<-syn$tableQuery('select * from syn21221980')$asDataFrame()
  
  allDat<-lapply(tabs$tableId,function(y){
    syn$tableQuery(paste('select Symbol,totalCounts,zScore,specimenID,diagnosis,tumorType,studyName,species,isCellLine,transplantationType from ',y))$asDataFrame()
  })
  full.dat<-do.call(rbind,allDat)
  full.dat$tumorType<-sapply(full.dat$tumorType,function(x) gsub('Malignant peripheral nerve sheath tumor','Malignant Peripheral Nerve Sheath Tumor',x))
  return(full.dat) 
  
}


#germline CSV calls
#are we using thesee?
getGermlineCsv<-function(syn,fileid,specimen){
  tab<- read.csv2(syn$get(fileid)$path,sep='\t')%>%
    select(gene_name,Tumor_VAF)%>%
    mutate(specimenID=specimen)
  return(tab)
}

##out of date xls data
#DEPRACATED
#' @param syn
#' @param fileid
#' @param indId
#' @import dplyr
#' @import readxl
processMergedXls<-function(syn,fileid,indId){
  library(readxl)
  library(dplyr)
#  id_list=list(BI3686='',)
  tab<- readxl::read_excel(syn$get(fileid)$path)
  tab<- tab%>%
    select(gene_name,ends_with("_var_count"))%>%
    tidyr::pivot_longer(cols=ends_with("_var_count"),names_to='Value',values_to='ADs')%>%
    tidyr::separate(2,sep='_',into=c('individualID','CountType'))%>%
    mutate(individualID=indId)%>%
    rowwise()%>%
    mutate(specimenID=paste(individualID,CountType,sep='_'))%>%
    ungroup()%>%
    subset(CountType!='RNA')%>% ##rna type gets lost in this parsing
    select(Symbol='gene_name',individualID,specimenID,AD='ADs')
  
  return(tab)
}

#' 
#' james said this is the format for future data
#' @import dplyr
#' @import tidyr
#' @import biomaRt
getNewSomaticCalls<-function(syn,fileid,specimen){
    library(dplyr)
    library(biomaRt)
#  print(fileid)
  tab<-read.csv2(syn$get(fileid)$path,sep='\t')%>%
    tidyr::separate(HGVSc,into=c('trans_id','var'))%>%
    mutate(trans_id=stringr::str_replace(trans_id,'\\.[0-9]+',''))%>%
    tidyr::separate(HGVSp,into=c('prot_id','pvar'))%>%
    mutate(prot_id=stringr::str_replace(prot_id,'\\.[0-9]+',''))
  
  ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")

  pmap<-biomaRt::getBM(mart=ensembl,attributes=c('ensembl_transcript_id','ensembl_peptide_id','hgnc_symbol'))
  ftab<-tab%>%rename(ensembl_transcript_id='trans_id')%>%left_join(pmap)
  res<-ftab%>%
    dplyr::select(c('hgnc_symbol',colnames(ftab)[grep('.AD$',colnames(ftab))]))%>%
    distinct()%>%
    mutate(specimenID=specimen)
  colnames(res)<-c('Gene','AD','nAD','specimenID')
  return(res)
}

#' getOldVariantData
#'@import dplyr
#'@import purrr

getOldVariantData<-function(syn){
  library(dplyr)
  files<-syn$tableQuery("SELECT * FROM syn11678418 WHERE \"name\" like '%merged.xlsx'")$asDataFrame()%>%
    dplyr::select(id,name,specimenID,individualID)%>%
    subset(is.na(specimenID))

  
  som.tab<-files%>%select(fileid='id',indId='individualID')%>%
    purrr::pmap_df(processMergedXls,syn)
  #purrr::map2_df(.f=processMergedXls,.x=files$id,.y=files$individualID)
  
  return(som.tab)
  ##create one giant table of variant allele frequences
  
}

#' getLatestVariantdata
#' @import dplyr
#' @import purrr
#' @export
getLatestVariantData<-function(syn){
  files<-syn$tableQuery("SELECT id,name,specimenID,individualID FROM syn21993642 WHERE ( (\"dataType\" = 'genomicVariants' ) )")$asDataFrame()%>%
    dplyr::select(id,name,specimenID,individualID)
  samps<-files%>%select(specimenID,individualID)%>%distinct()
  
  som.tab<-files%>%select(fileid='id',specimen='specimenID')%>%
    purrr::pmap_df(getNewSomaticCalls,syn)
  som.tab<-som.tab%>%subset(!is.na(Gene))%>%subset(Gene!="")%>%left_join(samps)
  return(som.tab)
}

#' until we have fully processed new data, grab missing samples from the old data 
mergeMutData<-function(mutData,newMutData){
  
}

