### Script assesses LUMCON allometry equations - seasonality, spatial variation

library(zoo)
library(plyr)

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150826.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150827.R")

# function to apply allometry to other time points (same plot). Generates a residual (observed - predicted)
predictBiomass <- function(allometryData = cwc.params, plotData = cwc, monthYear, plot, quadrat = 0.25, 
                           removeTrainingData = "TRUE", returnData = "FALSE") {
  # allometryData:  object with allometry parameters (cwc.params)
  # plotData:       object with plot data on plant heights and masses (cwc)
  # summaryData:    data with total biomass per plot, per time interval (maybe not needed)
  # monthYear:      month and year to use for estimating biomass. monthYear should be of the form, e.g., "May-13" for May 2013
  # plot:           which plot's data will be used
  # removeTrainingData: determines whether training data should also be modeled
  # returnData:     indicates whether predicted is returned by the function
#   
#   allometryData <- cwc.params
#   plotData      <- cwc
#   summaryData   <- CWC.plots # plot.totals has combined live/dead, CWC.plots distinguishes live/dead
#   
  countsAsTrue <- c("TRUE", "true", "T", "True")
  plotSize     <- quadrat^2
  
  # get allometry params
  a.live <- allometryData$coef.live[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                       ((as.character(allometryData$plot) %in% as.character(plot)))]    # coef
  b.live <- allometryData$exp.live[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                       ((as.character(allometryData$plot) %in% as.character(plot)))]    # exp
  a.dead <- allometryData$coef.dead[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                                     ((as.character(allometryData$plot) %in% as.character(plot)))]    # coef
  b.dead <- allometryData$exp.dead[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                                     ((as.character(allometryData$plot) %in% as.character(plot)))]
  
  # isolate plot x's data (and exclude month(s) used to parameterize model, if desired)
  if (removeTrainingData %in% countsAsTrue) {
    newData.live <- plotData[(as.character(plotData$site) %in% plot) & (as.character(plotData$type) %in% "LIVE") & (!plotData$monthYear %in% monthYear), ]
    newData.dead <- plotData[(as.character(plotData$site) %in% plot) & (as.character(plotData$type) %in% "DEAD") & (!plotData$monthYear %in% monthYear), ]
  } else {
    newData.live <- plotData[(as.character(plotData$site) %in% plot) & (as.character(plotData$type) %in% "LIVE"), ]
    newData.dead <- plotData[(as.character(plotData$site) %in% plot) & (as.character(plotData$type) %in% "DEAD"), ]
  }
  
  # apply allometry from t1 to other time points
  newData.live$predicted <- a.live * exp(b.live * newData.live$hgt)
  newData.dead$predicted <- a.dead * exp(b.dead * newData.dead$hgt)
  
  # error per stem (residual as percent of observed)
  newData.live$stem.err.live <- (newData.live$mass - newData.live$predicted) / newData.live$mass
  newData.dead$stem.err.dead <- (newData.dead$mass - newData.dead$predicted) / newData.dead$mass
  
  # error per plot
  # aggregate data
  live.data <- ddply(newData.live, .(site, time, monthYear, type), summarise,
                  mass.obs          =  sum(mass, na.rm = T) / plotSize,
                  mass.pred         =  sum(predicted, na.rm = T) / plotSize
                  )
  
  dead.data <- ddply(newData.dead, .(site, time, monthYear, type), summarise,
                     dead.mass.obs          =  sum(mass, na.rm = T) / plotSize,
                     dead.mass.pred         =  sum(predicted, na.rm = T) / plotSize
  )
  
  live.data$err.pct <- (live.data$mass.obs - live.data$mass.pred) / live.data$mass.obs
  dead.data$err.pct <- (dead.data$dead.mass.obs - dead.data$dead.mass.pred) / dead.data$dead.mass.obs
  
  # output median, IQR of errors
  output <- data.frame(stem.err.live    = median(newData.live$stem.err.live, na.rm = T), # errors per-stem
                      stem.err.dead     = median(newData.dead$stem.err.dead, na.rm = T),
                      plot.err.live     = median(live.data$err.pct, na.rm = T),
                      plot.err.dead     = median(dead.data$err.pct, na.rm = T),
                      stem.err.live.IQR = IQR(newData.live$stem.err.live, na.rm = T),
                      stem.err.dead.IQR = IQR(newData.dead$stem.err.dead, na.rm = T),
                      plot.err.live.IQR = IQR(live.data$err.pct, na.rm = T),
                      plot.err.dead.IQR = IQR(dead.data$err.pct, na.rm = T)
                      )
  if(returnData %in% countsAsTrue) {
    finalData <- join_all(list(newData.live, newData.dead), by = c(site, monthYear)) # probably won't work 
    rownames(finalData) <- 1:nrow(finalData)
    output <- list(summary = output, data = finalData)
  } else {
    output
  }
}


test <- predictBiomass(monthYear = "Jul-14", plot = "LUM1", removeTrainingData = "FALSE", returnData = "TRUE")

test[[1]]

# allometry data is in cwc.params; actual data is in cwc

### Part 1: use the allometric equation from each time point to estimate 
### biomass for every other time point, calculating error as the different 
### between predicted and observed biomass.


