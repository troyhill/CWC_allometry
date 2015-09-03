### Script assesses LUMCON allometry equations - seasonality, spatial variation

library(zoo)
library(plyr)

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150826.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150827.R")

# function to apply allometry to other time points (same plot). Generates a residual (observed - predicted)
predictBiomass <- function(allometryData, plotData, summaryData, monthYear, plot, quadrat = 0.25) {
  # allometryData:  object with allometry parameters (cwc.params)
  # plotData:       object with plot data on plant heights and masses (cwc)
  # summaryData:    data with total biomass per plot, per time interval (maybe not needed)
  # monthYear:      month and year to use for estimating biomass. monthYear should be of the form, e.g., "May-13" for May 2013
  # plot:           which plot's data will be used
  
  allometryData <- cwc.params
  plotData      <- cwc
  summaryData   <- CWC.plots # plot.totals has combined live/dead, CWC.plots distinguishes live/dead
  
  plotSize <- quadrat^2
  
  # get allometry params
  a.live <- allometryData$coef.live[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                       ((as.character(allometryData$plot) %in% as.character(plot)))]    # coef
  b.live <- allometryData$exp.live[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                       ((as.character(allometryData$plot) %in% as.character(plot)))]    # exp
  a.dead <- allometryData$coef.dead[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                                     ((as.character(allometryData$plot) %in% as.character(plot)))]    # coef
  b.dead <- allometryData$exp.dead[(as.character(allometryData$monthYear) %in% as.character(monthYear))  & 
                                     ((as.character(allometryData$plot) %in% as.character(plot)))]
  
  # isolate plot x's data (and exclude month used to parameterize model)
  newData.live <- plotData[(as.character(plotData$site) %in% plot) & (as.character(plotData$type) %in% "LIVE") & (!plotData$monthYear %in% monthYear), ]
  newData.dead <- plotData[(as.character(plotData$site) %in% plot) & (as.character(plotData$type) %in% "DEAD") & (!plotData$monthYear %in% monthYear), ]
  
  # apply allometry from t1 to other time points
  newData.live$predicted <- a.live * exp(b.live * newData.live$hgt)
  newData.dead$predicted <- a.dead * exp(b.dead * newData.dead$hgt)
  
  # error per stem (residual as percent of observed)
  stem.err.live <- (newData.live$mass - newData.live$predicted) / newData.live$mass
  stem.err.dead <- (newData.dead$mass - newData.dead$predicted) / newData.dead$mass
  median(stem.err.live, na.rm = T) # very skewed distribution
  median(stem.err.dead, na.rm = T)
  
  # error per plot
  # aggregate data
  live.data <- ddply(newData.live, .(site, time, monthYear, type), summarise,
                  mass.obs          =  sum(mass, na.rm = T) / plotSize,
                  mass.pred         =  sum(predicted, na.rm = T) / plotSize
                  )
  
  dead.data <- ddply(newData.dead, .(site, time, monthYear, type), summarise,
                     mass.obs          =  sum(mass, na.rm = T) / plotSize,
                     mass.pred         =  sum(predicted, na.rm = T) / plotSize
  )
  
  plot.err.live <- (live.data$mass.obs - live.data$mass.pred) / live.data$mass.obs
  plot.err.dead <- (dead.data$mass.obs - dead.data$mass.pred) / dead.data$mass.obs
  
  median(plot.err.live, na.rm = T)
  median(plot.err.dead, na.rm = T)
  
  # output median, IQR of errors
  output <- dataframe(stem.err.live = median(stem.err.live, na.rm = T),
                      stem.err.dead = median(stem.err.dead, na.rm = T),
                      plot.err.live = median(plot.err.live, na.rm = T),
                      plot.err.dead = median(plot.err.dead, na.rm = T),
                      stem.err.live.IQR = IQR(stem.err.live, na.rm = T),
                      stem.err.dead.IQR = IQR(stem.err.dead, na.rm = T),
                      plot.err.live.IQR = IQR(plot.err.live, na.rm = T),
                      plot.err.dead.IQR = IQR(plot.err.dead, na.rm = T)
                      )
  output
}



# allometry data is in cwc.params; actual data is in cwc

### Part 1: use the allometric equation from each time point to estimate 
### biomass for every other time point, calculating error as the different 
### between predicted and observed biomass.


