##quick script to move mutational data to table


library(reticulate)
library(dplyr)
library(mpnstXenoModeling)




getLatestVariantData <- function(syn){
  library(dplyr)
  library(purrr)
   files<-syn$tableQuery("SELECT id,name,specimenID,individualID FROM syn21993642 WHERE ( (\"dataType\" = 'genomicVariants' ) )")$asDataFrame()%>%
     dplyr::select(id,name,specimenID,individualID)

   samps<-files%>%dplyr::select(specimenID,individualID)%>%
     distinct()

   som.tab<-files%>%
     dplyr::select(fileid='id',specimen='specimenID')%>%
     purrr::pmap_df(getNewSomaticCalls,syn)
   som.tab<-som.tab%>%subset(!is.na(Gene))%>%subset(Gene!="")%>%
     left_join(samps)
   return(som.tab)
}

syn<-reticulate::import('synapseclient')$login()


data<-mgetLatestVariantData(syn)

res = synTableStore(data,'Updated Somatic Variant Data','syn21984813')