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
  ggsave(paste(gene,'expresion.png',sep=''))
}


plotExpressionByResponse<-function(genelist,drugresponse){}

#' plotHistogram of drug data
#' @export
#' @param drugData
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



#' plotGrowthForBatch
#' uses Xeva pacage to plot growth of treatment vs. control
#' @export
#' @import Xeva
plotGrowthForBatch<-function(xeva.obj,batch){
  plotPDX(xeva.obj,batch=batch,SE.plot='ribbon',
          vol.normal=FALSE,
          title=batch,
          control.name='vehicle')
}



plotWaterfall<-function(drug.tab,mut.tab,gene,drug){
 
  
}