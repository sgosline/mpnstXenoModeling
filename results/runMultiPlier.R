##run multiPLIER

require(mpnstXenoModeling)

library(PLIER)
library(dplyr)
library(tidyr)
#rc2path = "https://ndownloader.figshare.com/files/10881866"

##download these behemoth separately
#https://figshare.com/articles/recount_rpkm_RData/5716033/4

#then we want to looad this file
plier.results <- readRDS("recount_PLIER_model.RDS")


loadPDXData()



exprs.mat<-rnaSeq%>%
  dplyr::select(specimenID,'Symbol','totalCounts')%>%
  pivot_wider(values_from=totalCounts, names_from=specimenID,
              values_fn=list(totalCounts=mean),
              values_fill=list(totalCounts=0.01))%>% #no zeroes allowed!
  tibble::column_to_rownames('Symbol')%>%as.matrix()+0.01

zvar<-which(apply(exprs.mat,1,var)==0)
if(length(zvar)>0)
  exprs.mat<-exprs.mat[-zvar,]
source("../../multi-plier/util/plier_util.R")
pat.recount.b <- GetNewDataB(exprs.mat = exprs.mat,
                             plier.model = plier.results)


###now that all the data are loaded we can compute the correlation, regression, and random forest values

lv.df<-pat.recount.b%>%
  as.data.frame()%>%
  tibble::rownames_to_column("Latent Variable")%>%
  pivot_longer(-`Latent Variable`,names_to='specimenID',values_to="Loading")
