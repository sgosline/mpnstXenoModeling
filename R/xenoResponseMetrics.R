##compute response metrics using various types of claculations

###calculate metrics for each
#https://link.springer.com/article/10.1208/s12248-018-0284-8
#' @description SPI survival prolongation index calculates 
#' ratio of time it takes treated to tget to volume compared to control
#' @param tvThreshold - thresholl at which to measure
#' @param treatedTab - table of treated values
#' @param contTab - table of control values
computeSPI<-function(treatedTab,contTab,tvThreshold){
  ttime=subset(treatedTab,volume<=tvThreshold)%>%
    subset(time>0)%>%
      group_by(model.id)%>%summarize(time=min(as.numeric(time)))
  ctime=subset(contTab,volume<=tvThreshold)%>%
    subset(time>0)%>%
      group_by(model.id)%>%
      summarize(time=min(as.numeric(time)))
  res= mean(ttime$time,na.rm=T)/mean(ctime$time,na.rm=T)#ttime$time/ctime$time 
  #return(data.frame(metric='SPI',param='volumeThreshold',value=tvThreshold,result=res))
  return(res)
  }

#' TGI - tumor growth index
computeTGI<-function(treatedTab,contTab,finalTimePoint){
  tvt=subset(treatedTab,time==finalTimePoint)%>%
    group_by(model.id)%>%summarize(vol=max(as.numeric(volume)))
  
  tvc=subset(contTab,time==finalTimePoint)%>%
    group_by(model.id)%>%summarize(vol=max(as.numeric(volume)))
  tv0=subset(contTab,time==0)%>%
    group_by(model.id)%>%summarize(vol=max(as.numeric(volume)))
  
  res=(mean(tvc$vol)-mean(tvt$vol))/(mean(tvc$vol)-mean(tv0$vol))
  
#  return(data.frame(metric='TGI',param='finalTimePoint',value=finalTimePoint,result=res))
 return(res) 
}

#' AUC - area under the curve
#' computes difference between AUC of treated vs control normalized by treatment
#' @import Xeva
computeAUC<-function(treatedTab,contTab){
  #https://link.springer.com/article/10.1208/s12248-018-0284-8
  library(Xeva)
  
  tauc=treatedTab%>%mutate(volume=as.numeric(volume))%>%
         group_by(model.id)%>%
        group_map(~ unlist(Xeva::AUC(.x$time,.x$volume))[['value']],.keep=TRUE)
  cauc=contTab%>%mutate(volume=as.numeric(volume))%>%
    group_by(model.id)%>%
    select(time,volume)%>%
    group_map(~ unlist(Xeva::AUC(.x$time,.x$volume))[['value']],.keep=TRUE)
  #sprint(tauc)
  ret= (mean(as.numeric(unlist(cauc)),na.rm=T)-mean(as.numeric(unlist(tauc)),na.rm=T))/mean(as.numeric(unlist(cauc)),na.rm=T)
  return(ret)
  }

#' GRI - growth rate inhibition
#' 
computeGRI<-function(){
  #https://link.springer.com/article/10.1208/s12248-018-0284-8
}

#' statsForDrugPatient
#' @description Function that calculates all stats for each drug/patient combo
#' @export
statsForDrugPatient<-function(indivId,treat){
  ptab<-subset(drugData,individualID==indivId)
  
  ttab<-subset(ptab,drug==treat)
  ctab<-subset(ptab,drug=='vehicle')
  
  nzt<-subset(ttab,time>10)%>%subset(volume>0)
  nzc<-subset(ctab,time>10)%>%subset(volume>0)
  
  minVol=max(min(as.numeric(nzt$volume),na.rm=T),
    min(as.numeric(nzt$volume),na.rm=T))
  maxTime=min(max(as.numeric(ttab$time),na.rm=T),
              max(as.numeric(ctab$time),na.rm=T))
  
  
  return(list(AUC=computeAUC(ttab,ctab),
              SPI=computeSPI(ttab,ctab,minVol),
              TGI=computeTGI(ttab,ctab,maxTime)))
  
}
