
## This package uses the synapse web platform at http://synapse.org
# you must first register and update your config

##we use the synapse python client and reticulate

##

#####Synapse python client wrapper files



#condaenv="C:\\Users\\gosl241\\OneDrive - PNNL\\Documents\\GitHub\\amlresistancenetworks\\renv\\python\\r-reticulate\\"

#' Logs into Synapse using local information
#' @import reticulate
#' @return Synapse login python entity
#' @export
synapseLogin<-function(){
  library(reticulate)
 # reticulate::use_condaenv(condaenv)
  syn=reticulate::import('synapseclient')
  sync=syn$login()
}

#' Synapse store 
#' @import reticulate
#' @param path to file
#' @param parentId of folder to store
#' @export
synapseStore<-function(path,parentId){
  library(reticulate)
 # reticulate::use_condaenv(condaenv)
  
  synapse=reticulate::import('synapseclient')
  sync=synapse$login()
  sync$store(synapse$File(path,parentId=parentId))
}


#' Synapse table store
#' REQUIRES PANDAS ON YOUR PYTHON PATH
#' @param table to store on synapse
#' @param tabname name of table
#' @param parentId id of project
#' @import reticulate
#' @export
synTableStore<-function(tab,tabname,parentId='syn22128879'){
  #we have to first write the table to a file, then build it and store it
  library(reticulate)
  print(head(tab))
  fpath=write.table(tab,file='tmp.csv',sep=',',row.names = FALSE,quote=FALSE)
 # reticulate::use_condaenv(condaenv)
  synapse=reticulate::import('synapseclient')
  
  tab<-synapse$build_table(tabname,parentId,'tmp.csv')
  sync=synapse$login()
  sync$store(tab)
  
}


#' query synapse table
#' This is how you get data from the project
#' @param tableid
#' @export
querySynapseTable<-function(tableid){
  syn=synapseLogin()
  res<-syn$tableQuery(paste('select * from',tableid))$asDataFrame()
  return(res)
}

