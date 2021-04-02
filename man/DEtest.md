```R
library('biomaRt')
library(tidyr)
library(dplyr)
library(edgeR)
library(tibble)
library(EnsDb.Hsapiens.v86)
library(viridis)
library(gplots)
library(pheatmap)
library(synapser)
library(synapserutils)
#library(leapr)
library(ggplot2)
library(ggridges)
library(org.Hs.eg.db)
library(clusterProfiler)
# https://shiring.github.io/genome/2016/10/23/AnnotationDbi
```


```R
synLogin("","")
```

    Welcome, jessbade!


    NULL



```R

```


```R
#split df by subsets
group.split <- split(var.ID,var.ID$group)
age.split <- split(var.ID,var.ID$age)
sex.split <- split(var.ID,var.ID$sex)
qual.split <- split(var.ID,var.ID$qual)

```

       ROW_ID ROW_VERSION     Sample Age    Sex      RNASeq
    1       1          14     WU-225  38   Male syn21054117
    2       2          14     WU-356  27   Male syn21054119
    3       3          10     WU-368  49 Female syn21054121
    4       4          10     WU-436  36   Male syn21054123
    5       5          17     WU-505  37 Female syn25152854
    6       6          17     WU-545  52   Male syn25152859
    7       8          18        MN2  67 Female syn25258854
    8      10          14  JHU 2-002   9   Male syn20503817
    9      11          17  JHU 2-023  25   Male syn25152842
    10     12          18  JHU 2-079  18 Female syn25258844
    11     13          18  JHU 2-055  10 Female syn25258842
    12     14          17  JHU 2-031  12   Male syn25152847
    13     16           9     WU-386  34   Male syn21054115
    14     19          18 JHU 2-009a  44   Male syn25258840
    15     21          18     WU-487  42   Male syn25258847
    16     22          18     WU-536  18 Female syn25258849
    17     23          18     WU-561  34 Female syn25258852
                                     PDX Drug Data             Somatic Mutations
    1  ["syn22024434","syn22024435","syn22024436"]               ["syn22092093"]
    2  ["syn22024430","syn22024431","syn22024432"]               ["syn22092094"]
    3                                         <NA>               ["syn22092095"]
    4                                         <NA>               ["syn22092096"]
    5  ["syn22024439","syn22024440","syn22024441"]                          <NA>
    6  ["syn22024442","syn22024443","syn22024444"] ["syn22136423","syn22136424"]
    7                                         <NA>                          <NA>
    8  ["syn22018368","syn22018369","syn22018370"] ["syn22136408","syn22136407"]
    9                                         <NA> ["syn22136412","syn22136411"]
    10                                        <NA>                          <NA>
    11                                        <NA>                          <NA>
    12 ["syn22024461","syn22024462","syn22024463"] ["syn22136405","syn22136404"]
    13                                        <NA>               ["syn22092092"]
    14                                        <NA>                          <NA>
    15                                        <NA>                          <NA>
    16                                        <NA>                          <NA>
    17                                        <NA>                          <NA>
       MicroTissueQuality MicroTissue Drug Data     MPNST NF1 Status
    1            Marginal                  <NA>   Primary        NF1
    2                Poor                  <NA>   Primary        NF1
    3            Marginal                  <NA>   Primary        NF1
    4            Marginal                  <NA>   Primary        NF1
    5                Poor                  <NA>   Primary        NF1
    6                <NA>                  <NA>   Primary        NF1
    7                Good                  <NA>   Primary        NF1
    8                <NA>                  <NA>   Primary        NF1
    9                <NA>                  <NA>   Primary        NF1
    10               <NA>                  <NA> Recurrent        NF1
    11               <NA>                  <NA>   Primary        NF1
    12           Marginal                  <NA>   Primary        NF1
    13               <NA>                  <NA>   Primary        NF1
    14               <NA>                  <NA> Recurrent        NF1
    15               <NA>                  <NA> Recurrent        NF1
    16               <NA>                  <NA> Recurrent        NF1
    17               <NA>                  <NA> Recurrent        NF1
              Location               Clinical Status Size
    1            Thigh                      Deceased   22
    2      Mediastinum                      Deceased   13
    3             Calf                           NED   14
    4            Thigh                           NED   16
    5           Pelvis                      Deceased   11
    6          Humerus                           NED   18
    7  Maxillary Sinus                       Unknown    6
    8           Pelvis                           NED    6
    9       Paraspinal                           NED    6
    10           Thigh                         Alive    5
    11      Scalp/Neck                           NED   NA
    12 Retroperitoneal                      Deceased   10
    13            Neck                           NED    7
    14           Thigh                         Alive   NA
    15       Left Neck                      Deceased    6
    16       Left Neck                           NED    6
    17     Left Pelvis Alive with metastatic disease   10
                  group    sex          age     qual
    WU.225     Deceased   Male      Over.18 Marginal
    WU.356     Deceased   Male      Over.18     Poor
    WU.368        Alive Female      Over.18 Marginal
    WU.436        Alive   Male      Over.18 Marginal
    WU.505     Deceased Female      Over.18     Poor
    WU.545        Alive   Male      Over.18     <NA>
    MN2            <NA> Female      Over.18     Good
    JHU.2.002     Alive   Male Under.and.18     <NA>
    JHU.2.023     Alive   Male      Over.18     <NA>
    JHU.2.079     Alive Female Under.and.18     <NA>
    JHU.2.055     Alive Female Under.and.18     <NA>
    JHU.2.031  Deceased   Male Under.and.18 Marginal
    WU.386        Alive   Male      Over.18     <NA>
    JHU.2.009a    Alive   Male      Over.18     <NA>
    WU.487     Deceased   Male      Over.18     <NA>
    WU.536        Alive Female Under.and.18     <NA>
    WU.561        Alive Female      Over.18     <NA>



```R
synapse_query_and_save <- function(synParams,synTableID,synTableFN) {
    results <- synTableQuery(toString(sprintf('select %s from %s where RNASeq IS NOT NULL',synParams,synTableID)
    syn.df <- as.data.frame(results)
    return(syn.df)
    # OPTIONAL local save
    # write.table(syn.df,file=synTableFN,sep=',')
    # results <- read.table(synTableFN,sep=',',header=TRUE)
    # results <- results[!(is.na(results$RNASeq) | results$RNASeq==""), ]                                      
}
```


```R
syn.df <- synapse_query_and_save('*','syn24215021','MPNST_PDX_MODELS.csv')
```


```R
synapse_table_format <- function(df) {
    #Modify table naming for data subsetting
    names(df)[names(df) == "Clinical Status"] <- "Clinical.Status"
    df$Clinical.Status = gsub("NED","Alive",df$Clinical.Status)
    df$Clinical.Status = gsub("Alive with metastatic disease","Alive",df$Clinical.Status)
    df$Clinical.Status = gsub("Unknown",NA,df$Clinical.Status)
    df$UpdatedTissueQuality <- df$MicroTissueQuality
    df$UpdatedTissueQuality <- gsub("Marginal", "Good",df$UpdatedTissueQuality)
    df$Age <- ifelse(df$Age < 19,  sub('[0-9]*', 'Under.and.18',df$Age), gsub('[0-9]*', 'Over.18',df$Age))
    df$Sample = gsub(" ",".",df$Sample)
    df$Sample = gsub("-",".",df$Sample)
    return(df)
}
```


```R
df <- synapse_table_format(syn.df)
```


```R
#Get IDs from synapse table
synIDs <- df$RNASeq
sampleIDs <- df$Sample

#create factors of clinical variables in synapse table
group <- factor(df$Clinical.Status)
sex <- factor(df$Sex)
age <- factor(df$Age)
qual <- factor(df$UpdatedTissueQuality)

#generate clinical variable partitions for heatmap legend
var.ID <- data.frame(group)
var.ID$sex <- sex
var.ID$age <- age
var.ID$qual <- factor(df$MicroTissueQuality)
rownames(var.ID) <- sampleIDs
```


```R
synget_quants_and_format <- function(df) {
    entity <- lapply(df$RNASeq, function(x) synGet(x))
    for (x in 1:length(sampleIDs)) { 
        file <- entity[[x]]
        qnt.table <- read.table(file$path,header=T)
        fn <- sampleIDs[x]
        edb <- EnsDb.Hsapiens.v86
        ensembl.trans <- qnt.table$Name
        TXNAME <- gsub("\\..*","",ensembl.trans)
        trans.table <- cbind(qnt.table,TXNAME)
        geneIDs <- ensembldb::select(edb, keys=TXNAME, keytype = "TXNAME", columns = c("GENEID", "GENENAME"))
        master.table <- left_join(geneIDs,trans.table,by=c("TXNAME"))%>%
            dplyr::select(TXID,GENEID,GENENAME,TPM,NumReads)# %>%
            #dplyr::rename(x=NumReads)
        names(master.table)[4] <- sprintf("%s.TPM",fn)
        names(master.table)[5] <- sprintf("%s.Count",fn)
        assign(sprintf('mt.%s',x),master.table)
    }
    for (x in 1:length(sampleIDs)) {
       
   }
    return(final.table)
}
```


```R
final.table <- synget_quants_and_format(df)
```


```R
filter_by_low_expressed_genes <- function(final.table) {
    #Filter lowly expressed genes
    sum.table <- final.table %>% mutate(sum = rowSums(.[grep("TPM", names(.))]))# %>% select(-(location:13))
    keep<-sum.table[!(sum.table$sum<1),]
    row.names(keep) <- keep$TXID
    counts <- keep %>% dplyr:: select(contains("Count"))
    names(counts) <- sapply(strsplit(names(counts), ".C"), '[', 1)
    return(counts)
}
```


```R
counts <- filter_by_low_expressed_genes(final.table)
```


```R


#final.table <- cbindX(mt.1,mt.2,mt.3,mt.4,mt.5,mt.6,mt.7,mt.8,mt.9,mt.10)
final.table <- full_join(mt.1,mt.2,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.3,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.4,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.5,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.6,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.7,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.8,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.9,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.10,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.11,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.12,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.13,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.14,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.15,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.16,by=c("TXID","GENEID","GENENAME"))
final.table <- full_join(final.table,mt.17,by=c("TXID","GENEID","GENENAME"))
#final.table[is.na(final.table)] <- 0
#mt.1 %>% map_df(rev)

```


```R

```


```R
design_matrix_and_limma_analysis <- function(counts) {
    # create data structure to hold counts and subtype information for each sample.
    dge <- DGEList(counts)
    dge <- calcNormFactors(dge)

    # Specify a design matrix without an intercept term
    design <- model.matrix(~0 + qual + age + sex)

    # Clean up column names of design matrix
    colnames(design) <- gsub('group','',colnames(design))
    colnames(design) <- gsub('sex','',colnames(design))
    colnames(design) <- gsub('age','',colnames(design))
    colnames(design) <- gsub('qual','',colnames(design))
    names <- row.names(design)
    names <- as.numeric(gsub("[^0-9.]", "", names))
    subset.names <- sampleIDs[names]
    attr(design, "row.names") <- subset.names

    # Limma voom model fitting
    v <- voom(dge[,subset.names],design) #,plot=TRUE

    # Limma fit analysis
    fit <- lmFit(v, design)
    #fit <- contrasts.fit(fit, contrasts=contr.matrix)
    fit <- eBayes(fit,trend=TRUE)
    return(fit)
}
```


    Error in DGEList(counts): could not find function "DGEList"
    Traceback:




```R
fit <- design_matrix_and_limma_analysis(counts)
```


```R
get_topGenes_and_plot_heatmap <- function(fit,coefname,heatmapFN) {
    topGenes <- topTable(fit,number=Inf,coef=c(coefname),p.value=0.01) 
    summary(decideTests(fit))
    de.out <- rownames_to_column(topGenes,var='TXID')
    
    #gather read counts that are normalized to 1 via calcnormfactors()
    norm.counts <- dge$counts[,subset.names]
    norm.counts <- as.data.frame(norm.counts)
    norm.counts <- tibble::rownames_to_column(norm.counts,var='TXID')
    
    #combined differentially expressed txids, genenames, and normalized counts
    de.tx <- select(de.out,TXID)
    gene.id <- gene_logFC %>% select(TXID,GENENAME)
    de.table <- left_join(de.tx,gene.id,by="TXID")
    de.table <- left_join(de.table,norm.counts,by="TXID")

    #downselect by just GENENAME, counts
    de.table <- de.table %>% select(GENENAME,contains('WU'),contains('MN'),contains('JHU'))
    de.table <- de.table[!is.na(de.table$GENENAME),]
    de.table <- de.table[!duplicated(de.table$GENENAME), ]
    de.table <- tibble::remove_rownames(de.table)
    de.table <- tibble::column_to_rownames(de.table,var = 'GENENAME')
    
    h.t <- sapply(de.table,as.numeric)
    h.t <- as.matrix(h.t)

    rownames(h.t) = rownames(de.table)
    h.t <- h.t[rowSums(h.t)>0,]

    # Convert counts to zScore of relevant genes
    cal_z_score <- function(x){(x - mean(x)) / sd(x)}
    data_subset_norm <- t(apply(h.t, 1, cal_z_score))

    # Pairwise correlation between samples (columns)
    cols.cor <- cor(h.t, use = "pairwise.complete.obs", method = "pearson")
    # Pairwise correlation between rows (genes)
    rows.cor <- cor(t(h.t), use = "pairwise.complete.obs", method = "pearson")
    
    # Plot heatmap
    pheatmap(data_subset_norm,
         show_rownames=FALSE, 
         annotation_col=var.ID,
         clustering_distance_cols = as.dist(1 - cols.cor),
         clustering_distance_rows = as.dist(1 - rows.cor),
         filename=heatmapFN)
    dev.off()
    
    return(topGenes)
}
```


```R
topGenes <- get_topGenes_and_plot_heatmap(fit,'good','DE_GoodTissue_sex_age.pdf')
```


```R
do_gseGO_of_topGenes <- function(topGenes, rankedlistFN, gseoutFN) {
    tt <- tibble::rownames_to_column(topGenes, "TXID")

    gene_logFC <- left_join(tt,sum.table,by="TXID")
    gene_logFC <- gene_logFC[!duplicated(gene_logFC$GENENAME), ]
    ranked.gene.list <- gene_logFC %>% 
        dplyr::select(GENENAME,logFC) %>%
        tibble::remove_rownames() %>%
        column_to_rownames(var="GENENAME")

    #write.table(ranked.gene.list,rankedlistFN,sep='\t')

    geneList = sapply(ranked.gene.list,as.numeric)
    geneList = geneList[,1]
    names(geneList) = rownames(ranked.gene.list)
    geneList <- sort(geneList, decreasing = TRUE)

    gse = gseGO(geneList=geneList,ont="ALL",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,eps=0) #pAdjustMethod = 'BH',
    write.table(gse,toString(gseoutFN),sep='\t')
    return(gse)
}
```


```R
gseout <- do_gseGO_of_topGenes(topGenes, 'ranked.gene.list.Male.tsv', 'gsego_Male.tsv')
```


    Error in do_gseGO_of_topGenes(topGenes, "ranked.gene.list.Male.tsv", "gsego_Male.tsv"): could not find function "do_gseGO_of_topGenes"
    Traceback:




```R
plot_gseGO <- function(gseout,ClinVar,) {
    res <- as.data.frame(gse)
    prefix <- ClinVar
    all_gseaGO<-res %>% 
    dplyr::rename(pathway = 'Description') %>% 
    arrange(NES) %>% 
    dplyr::mutate(status = case_when(NES > 0 ~ "Up",
                                     NES < 0 ~ "Down"),
                  status = factor(status, levels = c("Up", "Down"))) %>% 
    group_by(status) %>% 
    top_n(20, wt = abs(NES)) %>% 
    ungroup() %>% 
    ggplot2::ggplot(aes(x=reorder(pathway, NES), y=NES)) +
    geom_bar(stat='identity', aes(fill=status)) +
    scale_fill_manual(values = c("Up" = "darkred", "Down" = "dodgerblue4")) +
    coord_flip() +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14),
          axis.text.y=element_text(size = 14),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none") +
    labs(title = "", y="NES") +#for some reason labs still works with orientation before cord flip so set y
    ggtitle(paste('All',prefix))
    ggsave(paste0(prefix,"_gseGO_plot.pdf"), all_gseaGO, height = 8.5, width = 11, units = "in")

    enrichplot::ridgeplot(gse,showCategory = 50,fill='pvalue')+ggplot2::ggtitle(paste0("GO Terms for ",prefix))
    ggplot2::ggsave(paste0(prefix,'_GO.pdf'),width=10,height=10)

    df<-as.data.frame(gse)%>%mutate(Condition=prefix)
}
```


```R
plot_gseGO(gseout,'Male')
```


```R

#pheatmap(data_subset_norm, annotation_col=group.ID,cellheight=8,color=rev(park_palette('GeneralGrant')),filename='heatmap_gseGO_group.pdf')

# library(Hmisc)
#pdf(file="heatmap_gseGO_group.pdf",width=10,height=30)
# v <- varclus(expr.dat,similarity="spearman")
# library(dendextend)
# dend <- as.dendrogram(v)
#heatmap.2(h.t,key=TRUE,density.info="none",trace="none",margins=c(8,16),cexCol=1,col=viridis,srtCol=45)#, symkey=TRUE,cexRow=1,) #,Colv = dend
#pheatmap(h.t)
#dev.off()
```


<strong>pdf:</strong> 2



```R

```


```R
#plots count value by sample for top Gene in ranked gene list
ix <- match(row.names(de.table)[1], row.names(counts) )
top1 <- ( as.vector( t( counts[ix,]) ) )
names( top1 ) = colnames(counts)
jpeg(file=sprintf('barplot.%s.jpeg',vars))
par(mar=c(12, 4.1, 4.1, 2.1))
barplot( top1 , las=3, xlab="")
dev.off()

```

    Warning message in min(x):
    "no non-missing arguments to min; returning Inf"
    Warning message in max(x):
    "no non-missing arguments to max; returning -Inf"



    Error in plot.window(xlim, ylim, log = log, ...): need finite 'ylim' values
    Traceback:


    1. barplot(top1, las = 3, xlab = "")

    2. barplot.default(top1, las = 3, xlab = "")

    3. plot.window(xlim, ylim, log = log, ...)



```R
#get sampleIDs of subsets
age.ID1 <- rownames(age.split$'18+')
age.ID2 <- rownames(age.split$Under.18)
group.ID1 <- rownames(group.split$Alive)
group.ID2 <-rownames(group.split$Deceased)
sex.ID1 <- rownames(sex.split$Male)
sex.ID2 <-rownames(sex.split$Female)
qual.ID1 <- rownames(qual.split$Poor)
qual.ID2 <-rownames(qual.split$Good)

#fetch count data for subsets to make subset tables
age.ID1.table <- counts[,age.ID1]
age.ID2.table <- counts[,age.ID2]
group.ID1.table <- counts[,group.ID1]
group.ID2.table <- counts[,group.ID2]
sex.ID1.table <- counts[,sex.ID1]
sex.ID2.table <- counts[,sex.ID2]
qual.ID1.table <- counts[,qual.ID1]
qual.ID2.table <- counts[,qual.ID2]
```


```R

```


```R
"""
input <- as.matrix(sapply(ranked.gene.list, as.numeric))
rownames(input) <- rownames(ranked.gene.list)

enrichment.order = leapR(geneset=ncipid, enrichment_method='enrichment_in_order',datamatrix=input,primary_columns='logFC')
print(head(enrichment.order))
"""
```


    Error in parse(text = x, srcfile = src): <text>:1:3: unexpected string constant
    6: print(head(enrichment.order))
    7: "
         ^
    Traceback:




```R

```


```R

```


```R

```
