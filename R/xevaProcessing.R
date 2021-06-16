##Xeva-specific analyses


#' This function takes a tidied table of data and makes into expresion dataset
#' can be used for RNASeq
#' TODO: expand for genomic data and latent variables
#'@import Biobase
#'@import tibble 
#'@import tidyr
#'@import dplyr
#'@export
#'
tidiedTableToExpressionSet<-function(tidied.tb,featureVar='totalCounts',
                                     nonFeatures=c('Symbol','zScore','parent')){
  require(Biobase)
  
  vfn=list(mean)
  names(vfn)<-featureVar
  
  vfi=list(0.0)
  names(vfi)<-featureVar
  
  assayData = tidied.tb%>%
    dplyr::select(Symbol,specimenID,!!featureVar)%>%
    subset(!is.na(Symbol))%>%
    tidyr::pivot_wider(names_from='specimenID',values_from=featureVar,values_fn=vfn,values_fill=vfi)%>%
    tibble::column_to_rownames('Symbol')%>%
    as.matrix()
  
  cvars<-setdiff(names(tidied.tb),c(featureVar,nonFeatures))
  print(cvars)
  rownames(tidied.tb)<-c()
  phenoData<-tidied.tb[,cvars]%>%
    distinct()%>%
    tibble::column_to_rownames('specimenID')%>%
    AnnotatedDataFrame()
  
  return(ExpressionSet(assayData,phenoData=phenoData))
}

#' formatDataToXeva
#'@import Xeva
#'@import dplyr
#'@export
formatDataToXeva<-function(){
  require(Xeva)
  require(dplyr)

  drugDat<-drugData
  controls=c('vehicle1','vehicle2','vehicle3')
  drugD = subset(drugDat,!drug%in%controls)
  contD = subset(drugDat,drug%in%controls)
  #contD$drug=rep('control',nrow(contD))
  
  
  expDesign =lapply(unique(drugD$batch),function(x){
    pat_drug=unlist(strsplit(x,split='_',fixed=T))
    treats = subset(drugD,batch==x)%>%
      dplyr::select(model.id)%>%
      unique()
    conts = subset(contD,Sample==pat_drug[1])%>%
      dplyr::select(model.id)%>%
      unique()
    
    list(batch.name=x,
         treatment=as.character(treats$model.id),
         control=as.character(conts$model.id))})
  
  model=drugDat%>%dplyr::select(model.id,patient.id='Sample')%>%
    distinct()%>%
    dplyr::mutate(tissue='MPNST')%>%
    dplyr::mutate(tissue.name='Malignant Peripheral Nerve Sheath Tumor')%>%
    ungroup()
  
  experiment=rbind(drugD,contD)%>%
    #mutate()%>%
    dplyr::select(-c(batch,Sample))
  
  drug = dplyr::select(experiment,drug)%>%distinct()%>%
    dplyr::mutate(treatment.type='single')%>%
    dplyr::rename(drug.id='drug')%>%
    dplyr::mutate(standard.name=drug.id)%>%
    dplyr::mutate(treatment.target='None')
  
  rnaDat = rnaSeq%>%distinct()%>%
    tidiedTableToExpressionSet()
  
  mutDat = varData%>%tidiedTableToExpressionSet(.,featureVar='AD')
  
  wesMap<-varData%>%
    dplyr::select(sample='individualID',biobase.id='specimenID')%>%
    distinct()%>%
    mutate(mDataType='mutation')
  
  rnaMap<-rnaSeq%>%
    dplyr::select(sample='individualID',biobase.id='specimenID')%>%
    distinct()%>%
    mutate(mDataType='RNASeq')
  
  modToBiobaseMap<-drugDat%>%
    dplyr::select(model.id,sample='Sample')%>%
    distinct()%>%
    left_join(rbind(rnaMap,wesMap),by='sample')%>%
    subset(!is.na(biobase.id))
  
  print(modToBiobaseMap)
  xeva.set = createXevaSet(name="MPNST PDX Data", 
                           model=as.data.frame(model), drug=as.data.frame(drug),
                           experiment=as.data.frame(experiment), expDesign=expDesign,
                           molecularProfiles=list(RNASeq = rnaDat,
                                                  mutation=mutDat),
                           modToBiobaseMap = as.data.frame(modToBiobaseMap))
  return(xeva.set)
}
