##quick script to move mutational data to table


library(reticulate)
library(dplyr)
library(MXM)
syn<-reticulate::import('synapseclient')$login()


data<-MXM::getLatestVariantData(syn)

res = synTableStore(data,'Updated Somatic Variant Data','syn21984813')