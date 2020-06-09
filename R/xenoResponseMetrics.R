##compute response metrics using various types of claculations

###calculate metrics for each
#https://link.springer.com/article/10.1208/s12248-018-0284-8
#' @description SPI survival prolongation index calculates 
#' ratio of time it takes treated to tget to volume compared to control
#' @param tvThreshold - thresholl at which to measure
#' @param treatedTab - table of treated values
#' @param contTab - table of control values
computeSPI<-function(treatedTab,contTab,tvThreshold){
  ttime=subset(ttab,volume<=tvThresholdl)%>%
      group_by(model.id)%>%summarize(time=min(as.numeric(time)))
  ctime=subset(ctab,volume<=tvThresholdl)%>%group_by(model.id)%>%summarize(time=min(as.numeric(time)))
  res=ttime$time/ctime$time 
  return(data.frame(metric='SPI',param='volumeThreshold',value=tvThreshold,result=res))
}

#' TGI - tumor growth index
computeTGI<-function(treatedTab,contTab,finalTimePoint){
  
}

#' AUC - area under the curve
#' computes difference between AUC of treated vs control normalized by treatment
#' @import Xeva
computeAUC<-function(treatedTab,contTab){
  #https://link.springer.com/article/10.1208/s12248-018-0284-8
  library(Xeva)
  
  tauc=treatedTab%>%mutate(volume=as.numeric(volume))%>%
         group_by(model.id)%>%
        group_map(~ Xeva::AUC(.x$time,.x$volume),.keep=TRUE)
  cauc=contTab%>%mutate(volume=as.numeric(volume))%>%
    group_by(model.id)%>%
    select(time,volume)%>%
    group_map(~ unlist(Xeva::AUC(.x$time,.x$volume)),.keep=TRUE)
  
  
  }

#' GRI - growth rate inhibition
#' 
computeGRI<-function(){
  #https://link.springer.com/article/10.1208/s12248-018-0284-8
}

#' statsForDrugPatient
#' @description Function that calculates all stats for each drug/patient combo
#' @export
statsForDrugPatient<-function(indivId='WU545',treat='doxorubicin'){
  ptab<-subset(drugData,individualID==indivId)
  
  ttab<-subset(ptab,drug==treat)
  ctab<-subset(ptab,drug=='vehicle')
  
  minVol=max(min(as.numeric(ttab$volume),na.rm=T),
    min(as.numeric(ctab$volume),na.rm=T))
  maxTime=min(max(as.numeric(ttab$time),na.rm=T),
              max(as.numeric(ctab$time),na.rm=T))
  
  

  
}
