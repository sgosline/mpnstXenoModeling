#these are a few helper functions that should help with data visualization
#the goal is to visualize expression of various genes

#'@export
#'@import ggplot2
#'@import cowplot
plotSingleGene<-function(gene='CD274'){
  all.dat<-getAllNF1Expression()
  dotplot<-subset(all.dat,Symbol==gene)%>%ggplot()+geom_point(aes(x=totalCounts,y=zScore,col=tumorType,shape=studyName))
  barplot<-subset(all.dat,Symbol==gene)%>%ggplot()+geom_boxplot(aes(x=studyName,y=totalCounts,fill=tumorType))
  cowplot::plot_grid(dotplot,barplot,nrow=2)
  ggsave(paste(gene,'expression.png',sep=''))
}


#' plotHistogram of drug data
#' @name plotDrugData
#' @export
#' @import ggridges
#' @import ggplot2
#' @param drugData table of drug data to be plotted
plotDrugData<-function(drugData){
  library(ggplot2)
  library(ggridges)
  p<-ggplot(drugData,aes(x=as.numeric(volume),fill=drug,y=individualID))+
    geom_density_ridges(alpha=0.5)+  scale_y_discrete(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_fill_viridis_d
  print(p)
  ggsave("allDrugVolumeScores.png")
  p
}



#' plotPDXTreatmentBySample
#' @param dt is the drugDatatable
#' @export
plotPDXTreatmentBySample<-function(dt){
  #print(dt$Sample)
  sample=dt$Sample[1]
  batch=dt$batch[1]
  #print(sample)
  tt<-dt%>%dplyr::select(time,volume,drug)%>%distinct()
  tm <-tt%>%group_by(time,drug)%>%summarize(mvol=median(volume),minVol=median(volume)-sd(volume),maxVol=median(volume)+sd(volume))%>%distinct()
  
  p<-ggplot(tm,aes(x=time,y=mvol,ymin=minVol,ymax=maxVol,col=drug,fill=drug))+geom_line()+ggtitle(sample)+geom_ribbon(alpha=0.25)
  #ggsave(paste0(sample,'PDXmodeling.png'),p)
  print(p)
  return(p)
}



#'plotPDXTreatmentByDrug
#'@param samps
#'@param drugs
#'@export
plotPDXTreatmentByDrug<-function(samps,drugs){
  
  
}

#'plotMTTreatmentByDrug
#'@param samps list of sample identifiers
#'@param drugs list of drug names
#'@export
plotMTTreatmentByDrug<-function(samps,drugs){
  
  ##first filter
  res<-subset(mt.meta,experimentalCondition%in%drugs)%>%
    subset(individualID%in%samps)
  #out = res[order(res$Conc),]
  
  #then get drug data
  res2 = mpnstXenoModeling::getMicroTissueDrugData(syn,res)
  
  ##then plot for each drug
  dplots<-lapply(drugs,function(drug){ 
    print(drug);
    generate_DR_plots(res2,drug)}
    )

}

#' plotTumorGrowthCorrelations
#' @param drugGeneCors the output of `drugMutationsCor`
#' @param minCor absolutely correlation to plot
#' @export
#' @import ggplot2
#' @import ggridges
#' 
plotTumorGrowthCorrelations<-function(drugGeneCors,minCor=0.8){
  library(ggplot2)
  library(ggridges)
  

  plotGeneDrug<-function(tab){
    sym=tab$Symbol[1]
    drug=tab$drug[1]
    met=tab$Metric[1]
    ggplot(tab,aes(x=AD,y=Value))+
      geom_point()+
      geom_text(aes(label=individualID))+
      ggtitle(paste(sym,'by',drug,met))
  }
  
 
  plots<-drugGeneCors%>%
    subset(abs(corVal)>minCor)%>%
    group_by(Symbol,drug,Metric)%>%
    group_map(~ plotGeneDrug(.x),.keep=TRUE)
  plots
  }


plotWaterfall<-function(drug.tab,mut.tab,gene,drug){
 
  
}