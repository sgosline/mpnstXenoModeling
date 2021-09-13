##get datasets from figshare and whatever drug screening data we have.

#'parses a list of synapse ids, like form `syn24215021`
#'@param list.column Name of column to select and unlist
#'@return list of lists
parseSynidListColumn<-function(list.column){

  lc <- unlist(list.column)
  res<-lapply(lc,function(x){
    unlist(strsplit(gsub('"','',gsub('[','',gsub(']','',x,fixed=T),fixed=T)),split=','))
  })
  return(res)

}

#' do_deseq2 import
#' Uses tximport functionality to properly re-count salmon based estimates and returns
#' counts
#' TODO: update to include more than just counts
#' @param file in sf format
#' @import BiocManager
#' @export 
do_deseq_import <- function(file) {
  
  if(!require('DESeq2')){
    BiocManager::install("DESeq2")
    library(DESeq2)
  }
  
  if(!require('tximportData')){
    BiocManager::install("tximportData")
    library(tximportData)
  }

  if(!require('tximport')){
    BiocManager::install("tximport")
    library(tximport)
  }

  dir <- system.file("extdata", package="tximportData")
  tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))

  qnt.table <- data.frame(counts=tximport::tximport(file,type='salmon',
                                                    tx2gene = tx2gene)$counts)%>%
    tibble::rownames_to_column('GENEID')

  return(qnt.table)
}

 

#' gets a list of synapse ids and binds them together
#' @param tab table of MPNST samples
#' @param syn synapse login client
#' @param colname name of column to select
#' @import readxl
#' @import tidyr
#' @import dplyr
#' @export
dataFromSynTable<-function(tab,syn,colname){


  samps <- tab$Sample
  print(colname)
  synids<-parseSynidListColumn(tab[,colname])
  names(synids)<-samps
  #print(head(tab))
  ##get the columns from the csvs
  schemas=list(`PDX Drug Data`=c('individual_id','specimen_id','compound_name','dose','dose_unit','dose_frequency',
              'experimental_time_point','experimental_time_point_unit',
              'assay_value','assay_units'),
              `Somatic Mutations`=c('Symbol','individualID','specimenID','AD'),
              #RNASeq=c('TXID','Symbol','TPM','NumReads'),
              RNASeq=c('GENEID','counts'),
              `Microtissue Drug Data`=c())

  ##the columns in the table we need
  othercols<-c('Sample','Age','Sex','MicroTissueQuality','Location','Size','Clinical Status')

  res<-lapply(samps,function(y){
    other.vals<-subset(tab,Sample==y)%>%
      dplyr::select(othercols)
    synds = synids[[y]]
    if(synds[1]=='NaN')
      return(NULL)

    full.tab<-do.call(rbind,lapply(synds,function(x){
      #print(x)
      x=unlist(x)
      path <-syn$get(x)$path
      fend<-unlist(strsplit(basename(path),split='.',fixed=T))
      fend <- fend[length(fend)]
      if(fend=='csv')
        tab<-read.csv(path,fileEncoding = 'UTF-8-BOM')
      else if(fend=='xlsx'){
        tab<-readxl::read_xlsx(path)
        if(colname=='Somatic Mutations')
          tab<- tab%>%
               dplyr::select(gene_name,ends_with("_var_count"))%>%
               tidyr::pivot_longer(cols=ends_with("_var_count"),names_to='Value',values_to='ADs')%>%
               tidyr::separate(2,sep='_',into=c('individualID','CountType'))%>%
               mutate(individualID=y)%>%
               rowwise()%>%
               mutate(specimenID=paste(individualID,CountType,sep='_'))%>%
               ungroup()%>%
               subset(CountType!='RNA')%>% ##rna type gets lost in this parsing
               dplyr::select(Symbol='gene_name',individualID,specimenID,AD='ADs')

      }else if(fend=='tsv'){
        tab<-getNewSomaticCalls(read.csv2(path,sep='\t',header=T),y)
      }
      else if(fend=='sf'){##add check from rnaseq data
        tab<-do_deseq_import(path)#%>%
        # dplyr::rename(counts=x)
       # print(head(tab))
      }else
        tab<-readxl::read_xls(path)
     # print(head(tab))
      tab<-tab[,schemas[[colname]]]%>%
        mutate(synid=x)
      #print(head(tab))
      data.frame(cbind(tab,other.vals))
    }))
    return(full.tab)
  })
  res <- do.call(rbind,res)
  return(res)
}

#' @name fixDrugData
#' @param drugData data frame of drug data to harmonize
#'@export
fixDrugData<-function(drugData){
  drugDat = #subset(drugData,individualID==specimenId)%>%
    drugData%>%
    dplyr::select(model.id='individual_id',Sample,drug,volume,time,specimen_id,batch='synid')%>%
    #dplyr::select(individual_id,drug='compound_name',volume='assay_value',time='experimental_time_point')%>%
    # dplyr::mutate(model.id=as.character(model.id))%>%
    # dplyr::mutate(model.id=stringr::str_replace(model.id,'_[0-9]',''))%>%
    #rowwise()%>%
    mutate(drug=tolower(drug))%>%
   # rename(batch=synid)%>%#paste(Sample,drug,sep='_'))%>%ungroup()%>%
    mutate(volume=as.numeric(volume))%>%
    subset(!is.na(volume))

  ##these files are supposed to be in the same batch -each file is doxo, everlo, vehicle
  b1<-drugDat$batch%in%c("syn22024434", "syn22024435","syn22024436")
  b2<-drugDat$batch%in%c('syn22024430','syn22024431','syn22024432')
  b3<-drugDat$batch%in%c('syn22024439','syn22024440','syn22024441')
  b4<-drugDat$batch%in%c('syn22024442','syn22024443','syn22024444')
  b5<-drugDat$batch%in%c("syn22018368","syn22018369","syn22018370")
  b6<-drugDat$batch%in%c('syn22024461','syn22024462','syn22024463')

  drugDat$batch[b1]<-'batch1'
  drugDat$batch[b2]<-'batch2'
  drugDat$batch[b3]<-'batch3'
  drugDat$batch[b4]<-'batch4'
  drugDat$batch[b5]<-'batch5'
  drugDat$batch[b6]<-'batch6'
  ##update the control drug name
  #%in%c(NA,'N/A','control')
  # inds0 = grep('vehicle',drugDat$drug)
  # print(inds0)
  # drugDat$drug[inds0]<-'vehicle0'
  # inds1 = which(is.na(drugDat$drug))
  # drugDat$drug[inds1]='vehicle1'
  # inds2 = grep('n/a',drugDat$drug)
  # drugDat$drug[inds2]='vehicle2'
  # inds3 = grep("control",drugDat$drug)
  # drugDat$drug[inds3]='vehicle3'
  #
  # ai<-c(inds0,inds1,inds2,inds3)
  # drugDat$model.id[c(ai)]<-paste(drugDat$model.id[ai],drugDat$drug[ai],sep='_')
  # ##first lets add a replicate to those without
  # inds = grep('_',drugDat$model.id)
  # #  print(inds)
  # inds<-union(inds,ai)
  # drugDat$model.id[-inds]<-paste(drugDat$model.id[-inds],drugDat$drug[-inds],'1',sep='_')
  ai<-drugDat$inds%in%c('vehicle','N/A','control',NA)
  drugDat$drug[ai]<-rep('control',length(ai))

  drugDat$time[which(drugDat$time<0)]<-0

  return(drugDat)
}

#' loadPDXData gets data from synapse and stores them as global variables
#' @export
#' @import reticulate
#' @import dplyr
#' @import tidyr
loadPDXData<-function(){
  library(reticulate)
  library(dplyr)
  syn<-reticulate::import('synapseclient')$login()

  ##updated to use harmonized data table
  data.tab<-syn$tableQuery('select * from syn24215021')$asDataFrame()


  varData<<-dataFromSynTable(data.tab,syn,'Somatic Mutations')


  drugData<<-dataFromSynTable(data.tab,syn,'PDX Drug Data')%>%
    rename(drug='compound_name',time='experimental_time_point',volume='assay_value')%>%
    fixDrugData()

  #now get RNA-Seq
  #update to use `RNAseq` column
  rnaSeq<<-dataFromSynTable(data.tab,syn,'RNASeq')
   
  #query microtissue drug data
  mt.meta <- syn$tableQuery('SELECT id,individualID,experimentalCondition,parentId FROM syn21993642 WHERE "dataType" = \'drugScreen\' AND "assay" = \'cellViabilityAssay\'')$asDataFrame()
  mt.meta<<- mt.meta[!(mt.meta$parentId == 'syn25791480' | mt.meta$parentId == 'syn25791505'),]

}



#' #'getPdxRNAseqData gets all rna seq counts for xenografts
#' #'#'@export
#' #'DEPRACATED
#' #'@param syn synapse item from
#' getPdxRNAseqData<-function(syn){
#' #  wu.rnaSeq = syn$tableQuery("SELECT * FROM syn21054125 where transplantationType='xenograft'")$asDataFrame()
#'   jh.rnaSeq = syn$tableQuery("SELECT * FROM syn20812185 where transplantationType='xenograft'")$asDataFrame()%>%
#'     subset(individualID=='2-002')
#'   jh.rnaSeq$individualID<-'JHU 2-002'
#'   #updated 6/8
#'   new.rnaSeq = syn$tableQuery("SELECT * from syn23667380 where transplantationType='xenograft'")$asDataFrame()
#' 
#'   com.cols=intersect(colnames(new.rnaSeq),colnames(jh.rnaSeq))%>%
#'     setdiff(c("ROW_ID","ROW_VERSION"))
#'   count.tab=rbind(jh.rnaSeq[,com.cols],new.rnaSeq[,com.cols])
#'  # count.tab$individualID<-sapply(count.tab$individualID,function(x) gsub('2-','JHU',x))
#' #  count.tab$specimenID<-sapply(count.tab$specimenID,function(x) gsub('2-','JHU',x))
#' 
#' 
#'   return(count.tab)
#' 
#' }

#' getAllNF1Expression 
#' @title getAllNF1Expression
#' collects all data from NF1 processed data in synapse
#' @export
#' @param syn object
getAllNF1Expression<-function(syn){
  tabs<-syn$tableQuery('select * from syn21221980')$asDataFrame()

  allDat<-lapply(tabs$tableId,function(y){
    syn$tableQuery(paste('select Symbol,totalCounts,zScore,specimenID,diagnosis,tumorType,studyName,species,isCellLine,transplantationType from ',y))$asDataFrame()
  })
  full.dat<-do.call(rbind,allDat)
  full.dat$tumorType<-sapply(full.dat$tumorType,function(x) gsub('Malignant peripheral nerve sheath tumor','Malignant Peripheral Nerve Sheath Tumor',x))
  return(full.dat)

}


# #germline CSV calls
# #are we using thesee?
# getGermlineCsv<-function(syn,fileid,specimen){
#   tab<- read.csv2(syn$get(fileid)$path,sep='\t')%>%
#     dplyr::select(gene_name,Tumor_VAF)%>%
#     mutate(specimenID=specimen)
#   return(tab)
# }



#' deseq2NormFilter
#' Utilizes the geometric mean to normalize the entire table to generate a normalized
#' list of RNAseq counts filtered for low abundance
#' Only does column normalization and returns DESeq2 object
#' @param data.table
#' #' @import BiocManager
#' @export 
deseq2NormFilter<-function(data.table){
  library(dplyr)
  
  if(!require('DESeq2')){
    BiocManager::install('DESeq2')
    library(DESeq2)
  }
  
  counts <- data.table%>%
    dplyr::select(GENEID,counts,Sample)%>%
    tidyr::pivot_wider(values_from='counts',names_from='Sample')%>%
    tibble::column_to_rownames('GENEID')%>%
    round()
  
  coldata<-data.table%>%
    dplyr::select(Sample,Sex,MicroTissueQuality,Location,Size,Age,`Clinical Status`)%>%
    distinct()%>%
  #  mutate(Clinical.Status=unlist(Clinical.Status))%>%
    tibble::column_to_rownames('Sample')
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ Sex)# + Clinical.Status)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds<-DESeq2::DESeq(dds)
  
  return(dds)

}

#' normDiffEx
#' Normalizes differential data using limma and voom
normDiffEx<-function(data.table){
  library(limma)
      # create data structure to hold counts and subtype information for each sample.
  counts <- data.table%>%
    dplyr::select(TXID,NumReads,Sample)%>%
    tidyr::pivot_wider(values_from='NumReads',names_from='Sample')%>%
    tibble::column_to_rownames('TXID')

    dge <- limma::voom(counts)

    # Specify a design matrix without an intercept term
    #design <- model.matrix(~ Clinical.Status + MicroTissueQuality + Age + Sex)
  
    #row.names(design) <- colnames(counts)[as.numeric(rownames(design))] #sampleIDs[names]
    #attr(design, "row.names") <- subset.names
    
    # Limma voom model fitting
   # v <- voom(dge[,row.names(design)],design,plot=TRUE)
    
    # Limma fit analysis
    fit <- limma::lmFit(dge)
    #fit <- contrasts.fit(fit, contrasts=contr.matrix)
    fit <- limma::eBayes(fit,trend=TRUE)

}

#' ##out of date xls data
#' #DEPRACATED
#' #' @param syn
#' #' @param fileid
#' #' @param indId
#' #' @import dplyr
#' #' @import readxl
#' processMergedXls<-function(syn,fileid,indId){
#'   library(readxl)
#'   library(dplyr)
#'   tab<- readxl::read_excel(syn$get(fileid)$path)
#'   tab<- tab%>%
#'     dplyr::select(gene_name,ends_with("_var_count"))%>%
#'     tidyr::pivot_longer(cols=ends_with("_var_count"),names_to='Value',values_to='ADs')%>%
#'     tidyr::separate(2,sep='_',into=c('individualID','CountType'))%>%
#'     mutate(individualID=indId)%>%
#'     rowwise()%>%
#'     mutate(specimenID=paste(individualID,CountType,sep='_'))%>%
#'     ungroup()%>%
#'     subset(CountType!='RNA')%>% ##rna type gets lost in this parsing
#'     dplyr::select(Symbol='gene_name',individualID,specimenID,AD='ADs')
#' 
#'   return(tab)
#' }

#' @title getNewSomaticCalls
#' james said this is the format for future data
#' @import dplyr
#' @import tidyr
#' @import BiocManager
getNewSomaticCalls<-function(tab,specimen){
    library(dplyr)
  if(!require('EnsDb.Hsapiens.v86')){
    BiocManager::install('EnsDb.Hsapiens.v86')
    library(EnsDb.Hsapiens.v86)
  }else{
    library('EnsDb.Hsapiens.v86')
  }
    library(ensembldb)
  
#  print(fileid)
  tab<-tab%>%#read.csv2(syn$get(fileid)$path,sep='\t')%>%
    tidyr::separate(HGVSc,into=c('trans_id','var'))%>%
    mutate(trans_id=stringr::str_replace(trans_id,'\\.[0-9]+',''))%>%
    tidyr::separate(HGVSp,into=c('prot_id','pvar'))%>%
    mutate(prot_id=stringr::str_replace(prot_id,'\\.[0-9]+',''))
  
  database <- EnsDb.Hsapiens.v86
  pmap <- ensembldb::select(database, keys=tab$trans_id, keytype = "TXNAME", columns = c("GENENAME"))
  
  ftab<-tab%>%dplyr::rename(TXNAME='trans_id')%>%left_join(pmap)

  res<-ftab%>%
    dplyr::select(c('GENENAME',colnames(ftab)[grep('.AD$',colnames(ftab))]))%>%
    distinct()%>%
    mutate(specimenID=specimen)%>%
    mutate(individualID=specimen)
  colnames(res)<-c('Symbol','AD','nAD','specimenID','individualID')
  return(res)
}

#' getOldVariantData
#'@import dplyr
#'@import purrr

getOldVariantData<-function(syn){
  library(dplyr)
  files<-syn$tableQuery("SELECT * FROM syn11678418 WHERE \"name\" like '%merged.xlsx'")$asDataFrame()%>%
    dplyr::select(id,name,specimenID,individualID)%>%
    subset(is.na(specimenID))

  som.tab<-files%>%select(fileid='id',indId='individualID')%>%
    purrr::pmap_df(processMergedXls,syn)
  #purrr::map2_df(.f=processMergedXls,.x=files$id,.y=files$individualID)

  return(som.tab)
  ##create one giant table of variant allele frequences

}

#' getLatestVariantdata
#' @import dplyr
#' @import purrr
#' @export
getLatestVariantData<-function(syn){
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

#' until we have fully processed new data, grab missing samples from the old data
mergeMutData<-function(mutData,newMutData){

}

#' getMicroTissueDrugData
#' @name getMicroTissueDrugData
#' @param syn
#' @param mtd
#' @import dplyr
#' @import tidyr
#' @export
getMicroTissueDrugData <- function(syn, mtd) {
  library(dplyr)
  library(tidyr)

  #ids is list of synapse ids
  ids<-mtd$id

  #indiv is list of patient IDs
  indiv<-mtd$individualID

  #sets filenames to names of ids
  names(indiv)<-ids

  res=do.call(rbind,lapply(names(indiv),function(x)
  {
   # print(x)
    tab<-read.csv(syn$get(x)$path,fileEncoding = 'UTF-8-BOM')
   # print(head(tab))
    tab%>%
    dplyr::select(DrugCol='compound_name', CellLine='model_system_name', Conc='dosage',
                  Resp='response', RespType='response_type', ConcUnit='dosage_unit') %>%
    tidyr::pivot_wider(names_from=RespType, names_sep='.', values_from=Resp) %>%
    dplyr::rename(Viabilities='percent viability')
  }))
  #res$`total cell count`=unlist(res$`total cell count`)
  #res$`live cell count`=unlist(res$`live cell count`)
  return(res)
}
