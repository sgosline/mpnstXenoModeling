###fix identifiers

require(synapser)
synLogin()
require(dplyr)

key.map<-synapser::synTableQuery("SELECT * FROM syn22088018")$asDataFrame()%>%
  select(oldId='private_id_3',individualID='public_id')

wu.rnaSeq=synapser::synTableQuery("SELECT * FROM syn21054125 where transplantationType='xenograft'")$asDataFrame()%>%
  rename(oldId='individualID')

new.tab<-wu.rnaSeq%>%left_join(key.map,by='oldId')%>%select(-oldId)%>%rowwise()%>%
  mutate(specimenID=paste(individualID,'tumor',sep='_'))


#now delete rows from the table

#and replace
allRows=synapser::synTableQuery("SELECT * FROM syn21054125")
deleted<-synapser::synDelete(allRows)
synapser::synStore(synapser::Table("syn21054125",select(new.tab,-c(ROW_ID,ROW_VERSION))))
