library(mpnstXenoModeling)

syn<-synapseLogin()

tab<-syn$tableQuery('SELECT name,id,readPair,readStrandOrigin FROM syn21993642 where "assay"=\'rnaSeq\' and "specimenID"=\'JH-2-002-CL\' and "fileFormat"=\'fastq\'')$asDataFrame()

library(dplyr)
library(tidyr)

csv <- tab%>%mutate(sample=stringr::str_replace(name,'_R[1|2].fastq.gz',''))

r1<-subset(csv,readPair==1)%>%select(sample,id,readStrandOrigin)%>%
  dplyr::rename(strandedness='readStrandOrigin')%>%
  mutate(fastq_1=stringr::str_c('syn://',id))%>%
  select(sample,fastq_1,strandedness)


r2<-subset(csv,readPair==2)%>%select(sample,id,readStrandOrigin)%>%
  dplyr::rename(strandedness='readStrandOrigin')%>%
  mutate(fastq_2=stringr::str_c('syn://',id))%>%
  select(sample,fastq_2,strandedness)

csv<-r1%>%full_join(r2)