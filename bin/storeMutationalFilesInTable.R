##quick script to move mutational data to table


library(reticulate)
library(dplyr)
library(mpnstXenoModeling)
syn<-reticulate::import('synapseclient')$login()


data<-getLatestVariantData(syn)

res = synTableStore(data,'Updated Somatic Variant Data','syn21984813')