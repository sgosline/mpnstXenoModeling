##'
##'MT plotting and fitting functions
##'
##'

#' TryFit
#' Takes the dose response data and fits it to various models
#' @param dose_response
#' @param fixed
#' @param names
#' @param nan.handle
#' @param curve_type
#' @return fitted curve
#' @export
TryFit <- function(dose_response, fixed = c(NA, NA, NA, NA), names = c(NA, NA, NA, NA), 
                   nan.handle = c("LL4", "L4"), curve_type){
  if (var(dose_response$Viabilities) == 0) {
    dose_response$Viabilties[nrow(dose_response)] <- dose_response$Viabilities[nrow(dose_response)] + 10^-10
  }
  #dose_response[dose_response == 0] <- 10^-10
  nan.handle <- match.arg(nan.handle)
  #, logDose = exp(10)
  drug.model <- tryCatch({
    drcmod(dose_response, LL.4(fixed = fixed, names = names), curve_type)
  }, warning = function(w) {
    if(nan.handle == "L4"){
      drcmod(dose_response, L.4(fixed = fixed, names = names), curve_type)
    } else {
      drcmod(dose_response, LL.4(fixed = fixed, names = names), curve_type)
    }
  }, error = function(e) {
    drcmod(dose_response, L.4(fixed = fixed, names = names), curve_type)
  })
  return(drug.model = drug.model)
}

#'drcmod
#'Modified dose response curve function
#'@param dose_resp
#'@param fctval
#'@param curve_type
#'@return model
#'@export
drcmod <- function(dose_resp, fctval, curve_type){
  if(!require('drc')){
    BiocManager::install('drc')
    library(drc)
  }
  temp <- drm(formula   = Viabilities ~ Conc
              , curveid   = CellLine
              , data      = dose_resp
              , fct       = fctval
              , na.action = na.omit
              , control   = drmc(errorm = FALSE)
  )
  summary(temp)
  return(temp)
}


#' generate_DR_plots
#' We need to fit the data to a curve to plot
#' @param res
#' @param drugID
#' @export
generate_DR_plots <- function(res, drugID) {
  
  res<-subset(res,DrugCol==drugID)
  if(!require('pROC')){
    BiocManager::install('pROC')
    library(pROC)
  }
  if(!require('data.table')){
    install.packages('data.table')
    library(data.table)
  }
  dr.df <- res %>% dplyr::select(Conc, Viabilities, CellLine) #%>% setDT()
  dr.df <- tidyr::unnest(dr.df, Viabilities)
  dr.df$Viabilities<- na_if(dr.df$Viabilities, "NA")
  dr.df <- dr.df[complete.cases(dr.df), ]
  dr.dt <- data.table(dr.df)
  dt2 <- data.table()
  dat.dt <- data.table(Drug=character(),CellLine=character(),Hill=numeric(),ec50=numeric(),
                       MinViability=numeric(),MaxViability=numeric(),ic50=numeric(),auc=numeric())
  # factor by CellLine
  dr.dt[,CellLine:=as.factor(CellLine)]
  scale.num <- nlevels(dr.dt$CellLine)
  # iterate over CellLine: Conc, Viabilities
  for (i in 1:nlevels(dr.dt$CellLine)) {
    # subset data.table by CellLine level
    dt <- dr.dt[CellLine %in% levels(CellLine)[i]]
    dt[, CellLine := as.character(CellLine)]
    # TryFit of subsetted data.table
    #dft <- dt[,.(Viabilities = unlist(Viabilities)), by = setdiff(names(dt), 'Viabilities')]
    df <- as.data.frame(dt)
    min_value <- min(df$Viabilities)
    max_value <- max(df$Viabilities)
    print(df)
    fit.LL4 <- TryFit(df, fixed = c(NA, min_value, max_value, NA), names = c("hill", "min_value", "max_value", "ec_50"), nan.handle = "L4", "CellLine")
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
    dt[, Pred := predict(object=fit.LL4)]
    sd.lst = dt[, list(SD=sd(Viabilities)/4),by=Conc]
    df <- as.data.frame(dt)
    cells <- unique(df$CellLine)
    # pROC::auc() https://cran.r-project.org/web/packages/pROC/pROC.pdf
    print(head(df$Conc))
    print(head(df$Viabilities))
    auc <- auc(df$Viabilities,df$Conc)
    meta2 <- data.table(Drug=drugID,CellLine=cells,Hill=coefs['hill'],ec50=coefs['ec_50'],
                        MinViability=min_value, MaxViability=max_value,ic50=ic_50,auc=as.numeric(auc))
    dat.dt <- rbind(dat.dt,meta2)
    dt = merge(dt,sd.lst)
    dt2 <- rbind(dt2, dt)
  }
  temp.plot = ggplot(dt2, aes(x=Conc, y=Viabilities, group=CellLine, color=CellLine, shape=CellLine)) +
    geom_ribbon(aes(y=Pred, ymin=Pred-SD, ymax=Pred+SD, fill=CellLine), alpha=0.2, color=NA) +
    geom_line(aes(y = Pred)) +
    scale_shape_manual(values=seq(0, scale.num)) +
    labs(x = "Conc log(M)", y = "Viabilities %", shape="Temp", color="Temp") +
    theme_bw() +
    ggtitle(paste("Dose-response curves for Drug:", drugID))
  return(list('plot'=temp.plot,'table'=dat.dt))
}
