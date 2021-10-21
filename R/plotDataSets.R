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
  pal<-c(nationalparkcolors::park_palette('GeneralGrant'),nationalparkcolors::park_palette('CraterLake'))
  
  sample=dt$Sample[1]
  batch=dt$batch[1]
  #print(sample)
  tt<-dt%>%dplyr::select(time,volume,drug)%>%distinct()
  tm <-tt%>%group_by(time,drug)%>%summarize(mvol=median(volume),
                                            minVol=median(volume)-sd(volume),
                                            maxVol=median(volume)+sd(volume))%>%distinct()
  
  p<-ggplot(tm,aes(x=time,y=mvol,ymin=minVol,ymax=maxVol,col=drug,fill=drug))+
    geom_line()+ggtitle(sample)+geom_ribbon(alpha=0.25)+scale_color_manual(values=pal)+scale_fill_manual(values=pal)
  #ggsave(paste0(sample,'PDXmodeling.png'),p)
  print(p)
  return(p)
}


#'plotMTTreatmentByDrug
#'@param dt mt data table
#'@export
plotMTTreatmentByDrug<-function(dt,sample){
  
  pal<-c(nationalparkcolors::park_palette('GeneralGrant'),nationalparkcolors::park_palette('CraterLake'))
  dt<-subset(dt,CellLine==sample)
  #batch=dt$batch[1]
  #print(sample)
  tt<-dt%>%dplyr::select(Conc,Viabilities,DrugCol)%>%distinct()
  tm <-tt%>%group_by(Conc,DrugCol)%>%summarize(mvol=median(Viabilities,na.rm=T),
                                            minVol=median(Viabilities,na.rm=T)-sd(Viabilities,na.rm=T),
                                            maxVol=median(Viabilities,na.rm=T)+sd(Viabilities,na.rm=T))%>%distinct()
  #print(tm)  
  p<-ggplot(tm,aes(x=Conc,y=mvol,ymin=minVol,ymax=maxVol,col=DrugCol,fill=DrugCol))+
    geom_line()+ggtitle(sample)+geom_ribbon(alpha=0.25)+scale_color_manual(values=pal)+scale_fill_manual(values=pal)
  #ggsave(paste0(sample,'PDXmodeling.png'),p)
  print(p)
  
  return(p)
  
}

#'plotMTTreatmentBySample
#'@param dt mt data table
#'@export
plotMTTreatmentBySample<-function(dt,drug){
  
  pal<-c(nationalparkcolors::park_palette('GeneralGrant'),nationalparkcolors::park_palette('CraterLake'))
  dt<-subset(dt,DrugCol==drug)
  #batch=dt$batch[1]
  #print(sample)
  tt<-dt%>%dplyr::select(Conc,Viabilities,CellLine)%>%distinct()
  tm <-tt%>%group_by(Conc,CellLine)%>%summarize(mvol=median(Viabilities,na.rm=T),
                                                minVol=median(Viabilities,na.rm=T)-sd(Viabilities,na.rm=T),
                                                maxVol=median(Viabilities,na.rm=T)+sd(Viabilities,na.rm=T))%>%
            distinct()
  
  p<-ggplot(tm,aes(x=Conc,y=mvol,ymin=minVol,ymax=maxVol,col=CellLine,fill=CellLine))+
    geom_line()+ggtitle(drug)+geom_ribbon(alpha=0.25)+scale_color_manual(values=pal)+scale_fill_manual(values=pal)
  #ggsave(paste0(sample,'PDXmodeling.png'),p)
  print(p)
  
  return(p)

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