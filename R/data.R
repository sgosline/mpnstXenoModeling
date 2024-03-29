# 
# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Loading data.R")
# }
# 
# .onLoad <- function(libname, pkgname) {
#   library(reticulate)
#   have_synapse <- reticulate::py_module_available("synapseclient")
#   if (!have_synapse)
#     reticulate::py_install("synapseclient")
#   
#   syn_client <<-
#     reticulate::import("synapseclient", delay_load = TRUE)$login()
# }
# 
# ##get datasets from figshare and whatever drug screening data we have.
# 
# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Loading data.R")
# }
# 
# syn_client <- NULL
# 
# .onLoad <- function(libname, pkgname) {
#   library(reticulate)
#   have_synapse <- reticulate::py_module_available("synapseclient")
#   if (!have_synapse)
#     reticulate::py_install("synapseclient")
#   
#   syn_client <<-
#     reticulate::import("synapseclient", delay_load = TRUE)$login()
# }

##get datasets from figshare and whatever drug screening data we have.
# TODO what is this function doing?
#' parses a list of synapse ids, like form `syn24215021`
#' @param list.column Name of column to select and unlist
#' @return list of lists
parseSynidListColumn <- function(list.column) {
  lc <- unlist(list.column)
  res <- lapply(lc, function(x) {
    unlist(strsplit(gsub('"', '', gsub('[', '', gsub(']', '', x, fixed = T), fixed = T)), split = ','))
  })
  return(res)
}

#' do_deseq2 import
#' Uses tximport functionality to properly re-count salmon based estimates and returns
#' counts
#' TODO: update to include more than just counts
#' @param file in sf format
#' @export
do_deseq_import <- function(file) {
  if (!require('DESeq2')) {
    BiocManager::install("DESeq2")
    library(DESeq2)
  }
  
  if (!require('tximportData')) {
    BiocManager::install("tximportData")
    library(tximportData)
  }
  
  if (!require('tximport')) {
    BiocManager::install("tximport")
    library(tximport)
  }
  
  dir <- system.file("extdata", package = "tximportData")
  tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
  
  qnt.table <-
    data.frame(counts = tximport::tximport(file, type = 'salmon',
                                           tx2gene = tx2gene)$counts) %>%
    tibble::rownames_to_column('GENEID')
  
  return(qnt.table)
}

#' gets a list of synapse ids and binds them together
#' @param tab table of MPNST samples
#' @param colname name of column to select

#' @export
dataFromSynTable <- function(tab, colname) {
  samps <- tab$Sample
  synids <- parseSynidListColumn(tab[, colname])
  names(synids) <- samps
  ##get the columns from the csvs
  schemas = list(
    `PDX Drug Data` = c(
      'individual_id',
      'specimen_id',
      'compound_name',
      'dose',
      'dose_unit',
      'dose_frequency',
      'experimental_time_point',
      'experimental_time_point_unit',
      'assay_value',
      'assay_units'
    ),
    `Somatic Mutations` = c('Symbol', 'individualID', 'specimenID', 'Tumor_AF', 'AD'),
    `RNASeq` = c('GENEID', 'counts'),
    `Microtissue Drug Data` = c(),
    `Incucyte drug Data` = c(
      'model_system_name',
      'compound_name',
      'compound_name_2',
      'dosage',
      'dosage_2...10',
      'dosage_unit',
      'response',
      'response_unit',
      'experimental_time_point',
      'experimental_time_point_unit',
      'replicate'
    )
  )
  ##RNASeq=c('TXID','Symbol','TPM','NumReads'),
  ##the columns in the table we need
  othercols <-
    c('Sample',
      'Age',
      'Sex',
      'MicroTissueQuality',
      'Location',
      'Size',
      'Clinical Status',
      'PRC2 Status')
  res <- lapply(samps, function(y) {
    other.vals <- subset(tab, Sample == y) %>%
      dplyr::select(othercols)
    synds = synids[[y]]
    if (synds[1] == 'NaN')
      return(NULL)
    
    full.tab <- do.call(rbind, lapply(synds, function(x) {
      x = unlist(x)
      path <- syn_client$get(x)$path
      
      # if (is.null(path)){
      #  print(paste("No access to :", x))
      #  return(NULL)
      # }
      
      fend <-
        unlist(strsplit(basename(path), split = '.', fixed = T))
      fend <- fend[length(fend)]
      if (fend == 'csv') {
        # Reading PDX tumor data processed by WU
        # TODO memory footprint optimize this func, specify colClasses as vector
        tab <-
          read.csv(
            path,
            fileEncoding = 'UTF-8-BOM',
            header = TRUE,
            row.names = NULL,
            check.names = TRUE
          )
      }
      else if (fend == 'xlsx') {
        # Reading incucyte data
        if (colname == 'Incucyte drug Data') {
          tab <-
            readxl::read_excel(
              path,
              col_types = c(
                "text",
                "text",
                "text",
                "text",
                "text",
                "text",
                "text",
                "text",
                "numeric",
                "numeric",
                "text",
                "text",
                "numeric",
                "text",
                "text",
                "text",
                "text",
                "text",
                "text",
                "numeric",
                "text",
                "numeric"
              )
            )
        }
        # Reading PDX tumor data not processed by WU
        else if (colname == 'PDX Drug Data') {
          tab <- tryCatch({
            readxl::read_excel(
              path,
              col_types = c(
                "text",
                "text",
                "text",
                "text",
                "text",
                "skip",
                "numeric",
                "text",
                "text",
                "text",
                "text",
                "skip",
                "numeric",
                "text",
                "text",
                "text",
                "numeric",
                "text",
                "text",
                "numeric",
                "text",
                "text"
              )
            )
          },
          error = function(cond) {
            readxl::read_excel(
              path,
              col_types = c(
                "text",
                "text",
                "text",
                "text",
                "text",
                "skip",
                "numeric",
                "text",
                "text",
                "text",
                "numeric",
                "text",
                "text",
                "numeric",
                "text",
                "text"
              )
            )
          })
        }
        else {
          tab <- readxl::read_excel(path)
        }
      }
      else if (fend == 'maf') {
        # Reading in somatic variants
        # TODO memory footprint optimize this func, specify colClasses as vector
        tab <- read.delim(path, header = TRUE, skip = 1)
        tab <- tab %>%
          tidyr::separate(Tumor_Sample_Barcode, c('individualID', NA), sep =
                            '_')
        # TODO check in if we should be subsetting IMPACT
        tab <- tab[tab$IMPACT %in% c('MODIFIER', 'HIGH'),]
        # TODO specimenID is assumed as individualID, will need to be updated
        tab$specimenID <- tab$individualID
        tab$Symbol <- tab$SYMBOL
        tab$AD <- tab$t_depth
        # tab<- tab%>%
        #    dplyr::select(gene_name,ends_with('_var_count'))%>%
        #    tidyr::pivot_longer(cols=ends_with('_var_count'),names_to='Value',values_to='ADs')%>%
        #    tidyr::separate(2,sep='_',into=c('individualID','CountType'))%>%
        #    mutate(individualID=y)%>%
        #    rowwise()%>%
        #    mutate(specimenID=paste(individualID,CountType,sep='_'))%>%
        #    ungroup()%>%
        #    subset(CountType!='RNA')%>% ##rna type gets lost in this parsing
        #    dplyr::select(Symbol='gene_name',individualID,specimenID,AD='ADs')
      }
      else if (fend == 'tsv') {
        # TODO memory footprint optimize this func, specify colClasses as vector
        tab <-
          getNewSomaticCalls(read.delim(path, header = TRUE), y)
      }
      else if (fend == 'sf') {
        # Reading in RNASeq data
        #add check from rnaseq data
        tab <- do_deseq_import(path)#%>%
        # dplyr::rename(counts=x)
      }
      else {
        tab <- readxl::read_xls(path)
      }
      tab <- tab[, schemas[[colname]]] %>%
        mutate(synid = x)
      data.frame(cbind(tab, other.vals))
    }))
    return(full.tab)
  })
  res <- do.call(rbind, res)
  
  return(res)
}

#' fixDrugData
#'
#' @param drugData data frame of drug data to harmonize
#'
fixDrugData <- function(drugData) {
  drugDat = #subset(drugData,individualID==specimenId)%>%
    drugData %>%
    dplyr::select(model.id = 'individual_id',
                  Sample,
                  drug,
                  volume,
                  time,
                  specimen_id,
                  batch = 'synid') %>%
    #dplyr::select(individual_id,drug='compound_name',volume='assay_value',time='experimental_time_point')%>%
    # dplyr::mutate(model.id=as.character(model.id))%>%
    # dplyr::mutate(model.id=stringr::str_replace(model.id,'_[0-9]',''))%>%
    #rowwise()%>%
    mutate(drug = tolower(drug)) %>%
    # rename(batch=synid)%>%#paste(Sample,drug,sep='_'))%>%ungroup()%>%
    mutate(volume = as.numeric(volume)) %>%
    subset(!is.na(volume))
  
  ##these files are supposed to be in the same batch -each file is doxo, everlo, vehicle
  b1 <- drugDat$batch %in% c("syn22024434", "syn22024435", "syn22024436")
  b2 <- drugDat$batch %in% c('syn22024430', 'syn22024431', 'syn22024432')
  b3 <- drugDat$batch %in% c('syn22024439', 'syn22024440', 'syn22024441')
  b4 <- drugDat$batch %in% c('syn22024442', 'syn22024443', 'syn22024444')
  b5 <- drugDat$batch %in% c("syn22018368", "syn22018369", "syn22018370")
  b6 <- drugDat$batch %in% c('syn22024461', 'syn22024462', 'syn22024463')
  
  drugDat$batch[b1] <- 'batch1'
  drugDat$batch[b2] <- 'batch2'
  drugDat$batch[b3] <- 'batch3'
  drugDat$batch[b4] <- 'batch4'
  drugDat$batch[b5] <- 'batch5'
  drugDat$batch[b6] <- 'batch6'
  ##update the control drug name
  # #%in%c(NA,'N/A','control')
  # inds0 = grep('vehicle',drugDat$drug)
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
  # inds<-union(inds,ai)
  # drugDat$model.id[-inds]<-paste(drugDat$model.id[-inds],drugDat$drug[-inds],'1',sep='_')
  ai <- drugDat$inds %in% c('vehicle', 'N/A', 'control', NA)
  
  drugDat$drug[ai] <- rep('control', length(ai))
  
  drugDat$time[which(drugDat$time < 0)] <- 0
  
  return(drugDat)
}

#' loadPDXData
#' gets data from synapse and stores them as global variables
#' @export
loadPDXData <- function() {
  library(dplyr)
  loadSynapse()
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  
  clin.tab <<- data.tab %>%
    dplyr::select(Sample,
                  Age,
                  Sex,
                  MicroTissueQuality,
                  MPNST,
                  Location,
                  `Clinical Status`,
                  Size) %>%
    distinct()
  return(clin.tab)
  
}



#' loadVariantData
#'
#' @return
#' @export
#'
#' @examples
loadVariantData <- function() {
  
  ##updated to use harmonized data table
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  varData <<- dataFromSynTable(data.tab, 'Somatic Mutations')
  return(varData)
}

#' loadRNASeqData
#'
#' @return rnaSeq data frame
#' @export
#'
#' @examples
loadRNASeqData <- function() {
  
  ##updated to use harmonized data table
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  #now get RNA-Seq
  #update to use `RNAseq` column
  rnaSeq <<- dataFromSynTable(data.tab, 'RNASeq') %>%
    mutate(`Clinical Status` = gsub(
      "NED",
      "Alive",
      gsub('Alive with metastatic disease', 'Alive', Clinical.Status)
    )) %>%
    tidyr::separate(GENEID,
                    into = c('GENE', 'VERSION'),
                    remove = FALSE)
  
  return(rnaSeq)
}

#' loadIncucyteData

#'
#' @return
#' @export
#'
#' @examples
loadIncucyteData <- function() {
  
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  
  #get incucyte data
  icyteData <<-
    dataFromSynTable(data.tab, 'Incucyte drug Data') %>%
    rowwise() %>%
    tidyr::unite(
      experimentalCondition,
      compound_name,
      compound_name_2,
      sep = ';',
      na.rm = TRUE,
      remove = TRUE
    ) %>%
    tidyr::unite(
      dosage,
      dosage,
      `dosage_2...10`,
      sep = ';',
      na.rm = TRUE,
      remove = TRUE
    ) %>%
    ungroup()
  return(icyteData)
}

#' loadPDXDrugData
#'
#' @return
#' @export
#'
#' @examples
loadPDXDrugData <- function() {
  ##updated to use harmonized data table
  
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  
  drugData <<- dataFromSynTable(data.tab, 'PDX Drug Data') %>%
    dplyr::rename(drug = compound_name,
                  time = experimental_time_point,
                  volume = assay_value) %>%
    fixDrugData()
  
  pdxDrugStats <<-
    syn_client$tableQuery('select * from syn25955439')$asDataFrame()
  return(list(dose = drugData, summary = pdxDrugStats))
}



#' loadMicrotissueDrugData
#'
#' @return
#' @export
#'
#' @examples
loadMicrotissueDrugData <- function() {
  
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  
  mtDrugData <<-
    syn_client$tableQuery('select * from syn26136282')$asDataFrame() %>%
    subset(CellLine %in% c(data.tab$Sample, 'MN-3'))
  return(mtDrugData)
}


#' loadMicrotissueMetadata loads and stores microtissue metata to global variable
#'
#'
#' @return NULL
#' @export
#'
#' @examples
loadMicrotissueMetadata <- function() {
  
  ##updated to use harmonized data table
  data.tab <<-
    syn_client$tableQuery("SELECT * FROM syn24215021 WHERE ( ( \"Microtissue Manuscript\" = 'true' ) )")$asDataFrame()
  
  ##mt data
  mt.meta <-
    syn_client$tableQuery(
      'SELECT id,specimenID,individualID,modelSystemName,experimentalCondition,experimentId,parentId FROM syn21993642 WHERE "dataType" = \'drugScreen\' AND "assay" = \'3D microtissue viability\' AND "fileFormat" = \'csv\' AND "parentId" not in (\'syn26433454\',\'syn25791480\',\'syn25791505\',\'syn26433485\',\'syn26433524\')'
    )$asDataFrame() %>%
    subset(individualID %in% c(data.tab$Sample, 'MN-3'))
  ##fix CUDC annotations
  cudc <- grep("CUDC", mt.meta$experimentalCondition)
  mt.meta$experimentalCondition[cudc] <- rep("CUDC-907", length(cudc))
  # Alphabetize drug combinations
  mt.meta <- mt.meta %>%
    tidyr::separate(experimentalCondition,
                    c('drug1', 'drug2'),
                    sep = ';',
                    remove = TRUE) %>%
    rowwise() %>%
    dplyr::mutate(experimentalCondition = paste0(sort(c(drug1, drug2)), collapse =
                                                   ';')) %>%
    ungroup() %>%
    dplyr::select(-c('drug1', 'drug2'))
  
  mt.meta <<- mt.meta
  return(mt.meta)
  
}





#' deseq2NormFilter
#' Utilizes the geometric mean to normalize the entire table to generate a normalized
#' list of RNAseq counts filtered for low abundance
#' Only does column normalization and returns DESeq2 object
#' @param data.table
#' @export
deseq2NormFilter <- function(data.table, newVar = NULL) {
  library(dplyr)
  
  if (!require('DESeq2')) {
    BiocManager::install('DESeq2')
    library(DESeq2)
  }
  
  counts <- data.table %>%
    dplyr::select(GENEID, counts, Sample) %>%
    tidyr::pivot_wider(values_from = 'counts', names_from = 'Sample') %>%
    tibble::column_to_rownames('GENEID') %>%
    round()
  
  coldata <- data.table %>%
    dplyr::select(Sample,
                  Sex,
                  MicroTissueQuality,
                  Location,
                  Size,
                  Age,
                  `Clinical Status`,
                  newVar) %>%
    distinct() %>%
    #  mutate(Clinical.Status=unlist(Clinical.Status))%>%
    tibble::column_to_rownames('Sample')
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = coldata,
                                        design = ~ Sex)# + Clinical.Status)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)
  
  return(dds)
}

#' normDiffEx
#' Normalizes differential data using limma and voom
normDiffEx <- function(data.table) {
  library(limma)
  # create data structure to hold counts and subtype information for each sample.
  counts <- data.table %>%
    dplyr::select(TXID, NumReads, Sample) %>%
    tidyr::pivot_wider(values_from = 'NumReads', names_from = 'Sample') %>%
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
  fit <- limma::eBayes(fit, trend = TRUE)
  return(fit)
}

#' ##out of date xls data
#' #DEPRACATED
#' #' @param syn
#' #' @param fileid
#' #' @param indId
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
#'
getNewSomaticCalls <- function(tab, specimen) {
  library(dplyr)
  if (!require('EnsDb.Hsapiens.v86')) {
    BiocManager::install('EnsDb.Hsapiens.v86')
    library(EnsDb.Hsapiens.v86)
  } else{
    library('EnsDb.Hsapiens.v86')
  }
  library(ensembldb)
  
  tab <- tab %>% #read.csv2(syn_clien$get(fileid)$path,sep='\t')%>%
    tidyr::separate(HGVSc, into = c('trans_id', 'var')) %>%
    mutate(trans_id = stringr::str_replace(trans_id, '\\.[0-9]+', '')) %>%
    tidyr::separate(HGVSp, into = c('prot_id', 'pvar')) %>%
    mutate(prot_id = stringr::str_replace(prot_id, '\\.[0-9]+', ''))
  
  database <- EnsDb.Hsapiens.v86
  pmap <-
    ensembldb::select(
      database,
      keys = tab$trans_id,
      keytype = "TXNAME",
      columns = c("GENENAME")
    )
  
  ftab <- tab %>% dplyr::rename(TXNAME = 'trans_id') %>% left_join(pmap)
  
  res <- ftab %>%
    dplyr::select(c('GENENAME', colnames(ftab)[grep('.AD$', colnames(ftab))])) %>%
    distinct() %>%
    mutate(specimenID = specimen) %>%
    mutate(individualID = specimen)
  colnames(res) <- c('Symbol', 'AD', 'nAD', 'specimenID', 'individualID')
  return(res)
}


#' getMicroTissueDrugData
#' @name getMicroTissueDrugData
#' @param mtd
#' @export
getMicroTissueDrugData <- function(mtd) {
  message('Loading MT data')
  library(dplyr)
  library(tidyr)
  
  drugs <- unique(mtd$experimentalCondition)
  res1 <-
    mtd %>% dplyr::select(id, individualID, experimentId, experimentalCondition) %>%
    apply(1, function(x) {
      y = x[['experimentalCondition']]
      is_combo = length(grep(';', y)) > 0
      is_dmso = y == 'DMSO'
      tab <-
        read.csv(syn_client$get(x[['id']])$path, fileEncoding = 'UTF-8-BOM')
      # p  rint(head(tab))
      ##TO  DO get this to work for combo data
      if (is_combo) {
        ttab <- tab %>%
          dplyr::select(
            'compound_name',
            'compound_name_2',
            CellLine = 'model_system_name',
            'dosage',
            'dosage_2',
            Resp = 'response',
            RespType = 'response_type',
            ConcUnit = 'dosage_unit'
          ) %>%
          rowwise() %>%
          mutate(DrugCol = paste(sort(c(
            compound_name, compound_name_2
          )), collapse = ';'),
          Conc = mean(c(dosage, dosage_2), na.rm = TRUE)) %>%
          dplyr::select(-c(compound_name, compound_name_2, dosage, dosage_2)) %>%
          tidyr::pivot_wider(
            names_from = RespType,
            names_sep = '.',
            values_from = Resp
          ) %>%
          dplyr::rename(Viabilities = 'percent viability') %>%
          unnest(cols = c(`total cell count`, `live cell count`, Viabilities))
      } else if (is_dmso) {
        ttab <- tab %>%
          dplyr::select(
            DrugCol = 'compound_name',
            CellLine = 'model_system_name',
            Conc = 'dosage',
            Resp = 'response',
            RespType = 'response_type',
            ConcUnit = 'dosage_unit'#,
       #     MeasID = 'measurement_id'
          ) %>%
          tidyr::pivot_wider(
            names_from = RespType,
            names_sep = '.',
            values_from = Resp
          ) %>%
          dplyr::rename(Viabilities = 'percent viability') %>% 
          unnest(cols = c(`total cell count`, `live cell count`, Viabilities))
      } else{
        ttab <- tab %>%
          dplyr::select(
            DrugCol = 'compound_name',
            CellLine = 'model_system_name',
            Conc = 'dosage',
            Resp = 'response',
            RespType = 'response_type',
            ConcUnit = 'dosage_unit'
          ) %>%
          tidyr::pivot_wider(
            names_from = RespType,
            names_sep = '.',
            values_from = Resp
          ) %>%
          dplyr::rename(Viabilities = 'percent viability') %>%
          unnest(cols = c(`total cell count`, `live cell count`, Viabilities))
      }
      ##override the cellLine in the file with the annotations from the table
      ttab <-
        ttab %>% mutate(CellLine = x[['individualID']], experimentId = x[['experimentId']])
      
      return(ttab)
    })
  res2 = do.call(rbind, res1)
  return(res2)
}
