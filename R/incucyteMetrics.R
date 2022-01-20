#' process_IncucyteDrugData
#'
#' @param master.table 
#' @param upload 
#' @param path 
#' @param parentID 
#'
#' @return
#' @export
#'
#' @examples
process_IncucyteDrugData <- function(master.table, upload=FALSE, path='.', parentID=NULL){
    require(data.table) 
  
   test=master.table%>%dplyr::select(compound_name,dosage)%>%distinct()
   
    print(test)
    master.table$compound_name <- as.factor(master.table$compound_name)
    #master.table$dosage <- as.factor(master.table$dosage)
    drugIDs <- levels(master.table$compound_name)
    drugIDs <- drugIDs[!drugIDs == 'DMSO']
    meta.dt <- data.table(compound_name=character(),model_system_name=character(),dosage=numeric(),Hill=numeric(),ec50=numeric(),MinViability=numeric(),MaxViability=numeric(),ic50=numeric(),auc=numeric())
    plt.list <- list()
    pi=1
    for (drug in drugIDs) {
    dres <- filter_byDrug(master.table, drug)
     dres <- split(dres,dres$dosage)
     for (res in dres){
       dose <- unique(res$dosage)
       plt = mpnstXenoModeling::generate_TR_plots(res,drug)
       meta.dt <- rbind(meta.dt, plt$table)
       plt.list[[pi]] <- cbind(plt$plot,drug)
       pi = pi+1
       
     }
    }
    
   meta.dt<-meta.dt%>%
     tidyr::replace_na(list(ic50=0.0))%>%as.data.frame()
   #fwrite(meta.dt,file="doseResponse_table.csv")
   pdf("doseResponsePlots.pdf")#, onefile = TRUE
   master.plot<-lapply(plt.list, tr_plot)
   print(master.plot)
   #invisible(lapply(master.plot, print))
   dev.off()
   #if (isTRUE(upload)) {
   #  synapseStore(file.path(path, "doseResponsePlots.pdf"),parentId=parentID)
   #  synapseStore(file.path(path, 'doseResponse_table.csv'),parentId=parentID)
      #ile.path(path, "doseResponse_table.csv"),parentId=parentID)
   #}
   return(meta.dt)
}


#' generate_TR_plots
#'
#' @param res 
#' @param drugID 
#'
#' @return
#' @export
#'
#' @examples
generate_TR_plots <- function(res, drugID) {
 #  require(data.table)
   require(dplyr)
   require(tidyr)
   tr.df <- res %>% dplyr::select(experimental_time_point, response, model_system_name, dosage) #%>% setDT()
   tr.df <- tidyr::unnest(tr.df, response)
   tr.df$response<- na_if(tr.df$response, "NA")
   tr.df <- tr.df[complete.cases(tr.df), ]
   tr.dt <- tr.df
   #tr.dt <- as.data.table(tr.df)
   #dt2 <- data.table()
   dat.dt <- data.frame(compound_name=character(),model_system_name=character(),dosage=numeric(),Hill=numeric(),ec50=numeric(),
                        MinViability=numeric(),MaxViability=numeric(),ic50=numeric(),auc=numeric())
   # factor by CellLine
  # tr.dt<-tr.dt%>%mutate(model_system_name=as.factor(tr.dt$model_system_name))
   scale.num <- length(unique(tr.dt$model_system_name))#nlevels(tr.dt$model_system_name)
   # iterate over CellLine: Conc, Viabilities
   for (i in unique(tr.dt$model_system_name)){#1:nlevels(tr.dt$model_system_name)) {
     # subset data.table by CellLine level
     dt <- tr.dt%>%subset(model_system_name==i)
     #dt<-mutate(model_system_name=as.character(dt$model_system_name))# model_system_name := as.character(model_system_name)]
     # TryFit of subsetted data.table
     #dft <- dt[,.(Viabilities = unlist(Viabilities)), by = setdiff(names(dt), 'Viabilities')]
     df <- as.data.frame(dt)
    min_value <- min(df$response)
     max_value <- max(df$response)
     fit.LL4 <- TryFit(df, fixed = c(NA, min_value, max_value, NA), names = c("hill", "min_value", "max_value", "ec_50"), nan.handle = "L4", "model_system_name")
     # hill from fit.LL4
     # ic50 code from https://rstudio-pubs-static.s3.amazonaws.com/378543_5b5bda32bf0541a485ccc49efbea680a.html
     coefs <- setNames(
       c(fit.LL4$coefficients, min_value, max_value),
       c("hill", "ec_50", "min_value","max_value"
       ))
     ic_50 <- with(as.list(coefs),
       exp(
         log(ec_50) + (1 / hill) * log(max_value / (max_value - 2 * min_value)) #log(ec_50)
       )
     )
     #cbind(auc,ED(fit.LL4,50))
     #rbindlist(dt.dat,auc)
     dt<-mutate(Pred=predict(object=fit.LL4))
     sd.lst = dt[, list(SD=sd(response)/4),by=experimental_time_point]
     df <- as.data.frame(dt)
     cells <- unique(df$model_system_name)
     dose <- unique(df$dosage)
     # pROC::auc() https://cran.r-project.org/web/packages/pROC/pROC.pdf
     auc <- auc(df$response,df$experimental_time_point)
     meta2 <- data.table(compound_name=drugID,model_system_name=cells,dosage=dose,Hill=coefs['hill'],ec50=coefs['ec_50'],
                        MinViability=min_value, MaxViability=max_value,ic50=ic_50,auc=as.numeric(auc))
     dat.dt <- rbind(dat.dt,meta2)
     dt = merge(dt,sd.lst)
     dt2 <- rbind(dt2, dt)
   }
   return(list('plot'=dt2,'table'=dat.dt))
 }
 

 # tr_plot<-function(dt2) {
 #     drug=unique(dt2$drug)
 #     dose=unique(dt2$dosage)
 #     plot=ggplot(dt2, aes(x=experimental_time_point, y=response, group=model_system_name, color=model_system_name)) + #, shape=model_system_name
 #             geom_ribbon(aes(y=Pred, ymin=Pred-SD, ymax=Pred+SD, fill=model_system_name), alpha=0.2, color=NA) +
 #             geom_line(aes(y = Pred)) +
 #             scale_shape_manual(values=seq(0, nlevels(dt2$model_system_name))) +
 #             labs(x = "Time", y = "Confluency %", shape="Temp", color="Temp") +
 #             theme_bw() +
 #             ggtitle(paste("Time-response curves for", drug, "at dose", dose, "nM"))
 #     return(plot)
 #     }
  