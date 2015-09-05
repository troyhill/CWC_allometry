### Script assesses LUMCON allometry equations - seasonality, spatial variation

library(zoo)
library(plyr)

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150826.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150827.R")

# function to apply allometry to other time points (same plot). Generates a residual (observed - predicted)
predictBiomass <- function(plotData = cwc, monthYear, plot, quadrat = 0.25, 
                           removeTrainingData = "FALSE", returnData = "TRUE", start_nls = 0.03) {
  # allometryData:  object with allometry parameters (cwc.params)
  # plotData:       object with plot data on plant heights and masses (cwc)
  # summaryData:    data with total biomass per plot, per time interval (maybe not needed)
  # monthYear:      month and year to use for estimating biomass. monthYear should be of the form, e.g., "May-13" for May 2013
  # plot:           which plot's data will be used
  # removeTrainingData: determines whether training data should also be modeled
  # returnData:     indicates whether predicted is returned by the function
#   
#   plotData      <- cwc
#   summaryData   <- CWC.plots # plot.totals has combined live/dead, CWC.plots distinguishes live/dead
   
  countsAsTrue <- c("TRUE", "true", "T", "True")
  plotSize     <- quadrat^2
  
  # isolate plot x's data (and exclude month(s) used to parameterize model, if desired)
  if (removeTrainingData %in% countsAsTrue) {
    newData <- plotData[(as.character(plotData$site) %in% plot) & (!plotData$monthYear %in% monthYear), ]
  } else {
    newData <- plotData[(as.character(plotData$site) %in% plot), ]
  }
  
  test.live <- nrow(newData[(newData$type %in% "LIVE") & (newData$monthYear %in% monthYear), ]) > 3
  test.dead <- nrow(newData[(newData$type %in% "DEAD") & (newData$monthYear %in% monthYear), ]) > 3
  
  # generate allometry params (allows multiple months to be combined)  
  # why doesn't this work?!?!
  if (test.live) {
    y.live <- newData$mass[(newData$monthYear %in% monthYear) & (newData$type %in% "LIVE")]
    x.live <- newData$hgt[(newData$monthYear %in% monthYear) & (newData$type %in% "LIVE")]
    y.live2 <- y.live[!is.na(y.live) & !is.na(x.live)]
    x.live2 <- x.live[!is.na(y.live) & !is.na(x.live)]
    live.coefs <- coef(model <- nls(y.live2 ~ I(a * exp(b * x.live2)), start = list(a = start_nls, 
                                                          b = start_nls)))
  } else {
    live.coefs <- c(NA, NA)
  }
  # force through zero?
#   plot(y.live2 ~ x.live2)
#   lines(1:150, y = live.coefs[1] * exp(live.coefs[2] * 1:150))
   
  if (test.dead) {
    y.dead <- newData$mass[(newData$monthYear %in% monthYear) & (newData$type %in% "DEAD")]
    x.dead <- newData$hgt[(newData$monthYear %in% monthYear) & (newData$type %in% "DEAD")]
    y.dead2 <- y.dead[!is.na(y.dead) & !is.na(x.dead)]
    x.dead2 <- x.dead[!is.na(y.dead) & !is.na(x.dead)]
    dead.coefs <- coef(model <- nls(y.dead2 ~ I(a * exp(b * x.dead2)), start = list(a = start_nls, 
        b = start_nls)))
  } else {
    dead.coefs <- c(NA, NA)
  }
#   plot(y.dead2 ~ x.dead2)
#   lines(1:150, y = dead.coefs[1] * exp(dead.coefs[2] * 1:150))

  a.live <- live.coefs[1]
  b.live <- live.coefs[2]
  a.dead <- dead.coefs[1]    # coef
  b.dead <- dead.coefs[2]

  # apply allometry from t1 to other time points
  newData$predicted <- NA
  newData$predicted[newData$type %in% "LIVE"] <- a.live * exp(b.live * newData$hgt[newData$type %in% "LIVE"])
  newData$predicted[newData$type %in% "DEAD"] <- a.dead * exp(b.dead * newData$hgt[newData$type %in% "DEAD"])
  # error per stem (difference as percent of observed)
  newData$stem.err <- (newData$predicted - newData$mass) / newData$mass
  
  # for each month, calculate error metrics
  for (i in 1:length(unique(newData$monthYear))) {
    targetTime <- unique(newData$monthYear)[i]
    subData    <- newData[newData$monthYear %in% targetTime, ]
    
    
    # this if statement accommodates months where one type of biomass is absent
    if (test.live) {
      # error per stem
      live.stem.err  <- median(subData$stem.err[subData$type %in% "LIVE"], na.rm = T)
      # error per plot
      live.mass.obs  <- sum(subData$mass[subData$type %in% "LIVE"], na.rm = T) / plotSize
      live.mass.pred <- sum(subData$predicted[subData$type %in% "LIVE"], na.rm = T) / plotSize
      live.err.pct   <- (live.mass.pred - live.mass.obs) / live.mass.obs
    } else {
      live.stem.err  <- NA 
      live.mass.obs  <- NA
      live.mass.pred <- NA
      live.err.pct   <- NA
    }
    
    # do the same for dead biomass
    if (test.dead) {
      dead.stem.err  <- median(subData$stem.err[subData$type %in% "DEAD"], na.rm = T)
      dead.mass.obs  <- sum(subData$mass[subData$type %in% "DEAD"], na.rm = T) / plotSize
      dead.mass.pred <- sum(subData$predicted[subData$type %in% "DEAD"], na.rm = T) / plotSize
      dead.err.pct   <- (dead.mass.pred - dead.mass.obs) / dead.mass.obs
    } else {
      dead.stem.err  <- NA 
      dead.mass.obs  <- NA
      dead.mass.pred <- NA
      dead.err.pct   <- NA
    }
  
    
  # output median, IQR of errors
  intData <- data.frame(plot            = as.character(plot),
                      trainingData      = as.character(paste(monthYear, collapse = ",")),
                      monthYear         = as.character(unique(newData$monthYear)[i]),
                      stem.err.live     = live.stem.err, # median percent errors on per-stem basis
                      stem.err.dead     = dead.stem.err,
                      stem.err.live.IQR = IQR(subData$stem.err[subData$type %in% "LIVE"], na.rm = T),
                      stem.err.dead.IQR = IQR(subData$stem.err[subData$type %in% "DEAD"], na.rm = T),
                      plot.err.live     = live.err.pct, # error as a percent of observed (positive nos mean predicted biomass is greater)
                      plot.err.dead     = dead.err.pct,
                      obs.biomass.live  = live.mass.obs, # observed biomass
                      obs.biomass.dead  = dead.mass.obs,
                      pred.biomass.live = live.mass.pred, # predicted biomass
                      pred.biomass.dead = dead.mass.pred
                      )
  
  if (i != 1) {
    output <- rbind(output, intData)
  } else {
    output <- intData
  }
  
  }
  
  if(returnData %in% countsAsTrue) {
    output <- list(summary = output, data = newData)
  } else {
    output
  }
}


test <- predictBiomass(monthYear = "Jul-14", plot = "LUM1")
test[[1]]
test <- predictBiomass(monthYear = "Jan-14", plot = "LUM1")
test[[1]]


# use all time points to define a plot's model. Positive errors = over-predict biomass
unique(cwc$monthYear)
test <- predictBiomass(monthYear = c(unique(cwc$monthYear)), plot = "LUM1")
test[[1]]


# allometry data is in cwc.params; actual data is in cwc

### Part 1: use the allometric equation from each time point to estimate 
### biomass for every other time point, calculating error as the different 
### between predicted and observed biomass.

# plot <- "LUM1"
# startVal <- 0.03

for (j in 1:length(unique(cwc$site))) {
  for (i in 1:length(unique(cwc$monthYear))) {
    temp <- predictBiomass(monthYear = c(unique(cwc$monthYear)[i]), plot = unique(cwc$site)[j], returnData = "FALSE")
  
    if (i != 1) {
      intData <- rbind(intData, temp)
    } else {
      intData <- temp
    }
  }
  temp <- predictBiomass(monthYear = c(unique(cwc$monthYear[cwc$site %in% plot])), plot = unique(cwc$site)[j], returnData = "FALSE")
  intData <- rbind(intData, temp)
  
  if (j != 1) {
    finalData <- rbind(finalData, intData)
  } else {
    finalData <- intData
  }
}


finalData$monthYear     <- as.yearmon(finalData$monthYear, "%B-%y")
finalData$plot          <- as.character(finalData$plot)
finalData$trainingData  <- as.character(finalData$trainingData)
finalData$trainingData[nchar(finalData$trainingData) > 6] <- "all"
finalData$biomass.error <- finalData$pred.biomass.live - finalData$obs.biomass.live # magnitude of biomass estimation error (g m2 yr)
finalData$dead.error    <- finalData$pred.biomass.dead - finalData$obs.biomass.dead    # magnitude of dead biomass estimation error (g m2 yr)

# finalData_allPlots <- finalData
finalData <- finalData[finalData$plot %in% paste0("LUM", 1:3), ]

# set legend in two columns

# full-scale plot
ggplot(data = finalData, aes(x = as.factor(monthYear), y = plot.err.live, col = trainingData)) + geom_point() + 
  theme_bw() + labs(x = "", y = "Estimation error (decimal fraction)") + scale_colour_discrete(name = "Training data") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + facet_grid(plot ~ .) + ylim(-2, 45)

ggplot(data = finalData, aes(x = as.factor(monthYear), y = plot.err.live, 
  col = trainingData)) + geom_point() + 
  theme_bw() + labs(x = "", y = "Estimation error (decimal fraction)") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  ylim(-10, 10)  + facet_grid(plot ~ .) + guides(colour = guide_legend(title = "Training data", override.aes = list(alpha = 1), ncol = 2))
# ggsave("C:/RDATA/SPAL_allometry/allom_errors_zoom1.png", width = 11, height= 5, units = "in", dpi = 300)



ggplot(data = finalData, aes(x = as.factor(monthYear), y = plot.err.live, colour = trainingData)) + 
  geom_point(alpha = 0.7) + 
  theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  labs(x = "", y = "Estimation error (decimal fraction)") + 
  guides(colour = guide_legend(title = "Training data", override.aes = list(alpha = 1), ncol = 2)) +
  ylim(-1, 1)  + facet_grid(plot ~ .) 
# ggsave("C:/RDATA/SPAL_allometry/allom_errors_zoom2.png", width = 11, height= 7, units = "in", dpi = 300)



ggplot(data = finalData[finalData$trainingData %in% "all", ], 
       aes(x = as.factor(monthYear), y = plot.err.live)) + # geom_point() + 
  theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + labs(x = "", y = "Estimation error (decimal fraction)") + 
  scale_colour_discrete(name = "Training data") + geom_text(aes(label = round(plot.err.live, 1)), hjust = 0.5, vjust = 0) +
  facet_grid(plot ~ .) + ylim(-0.5, 4.2)
# ggsave("C:/RDATA/SPAL_allometry/allom_errors_allData.png", width = 11, height= 8, units = "in", dpi = 300)


# how does estimation error change as function of observed biomass?
# for individual months, the relationship is very noisy and probably nonexistent
ggplot(data = finalData, aes(x = obs.biomass.live, y = plot.err.live, colour = as.factor(trainingData))) + # geom_point() + 
  theme_bw() + labs(x = expression("Observed biomass (g "%.%m^-2~")"), y = "Estimation error (decimal fraction)") + 
  guides(colour = guide_legend(title = "Training data", override.aes = list(alpha = 1), ncol = 2)) +
  geom_point() + ylim(-0.5, 4.2)
# ggsave("C:/RDATA/SPAL_allometry/allom_errors_vs_Obs.png", width = 10, height= 6, units = "in", dpi = 300)


# for 'all' data, the relationship is one of relatively low error at high biomasses, 
# and high (and highly variable) errors at low biomass.
# This is consistent with allometric equation which seems to show upward bias at low masses.
ggplot(data = finalData[finalData$trainingData %in% "all", ], aes(x = obs.biomass.live, y = plot.err.live)) + # geom_point() + 
  theme_bw() + labs(x = expression("Observed biomass (g "%.%m^-2~")"), y = "Estimation error (decimal fraction)") + 
  scale_colour_discrete(name = "Training data") + geom_point() +
  facet_grid(plot ~ .)
# ggsave("C:/RDATA/SPAL_allometry/allom_errors_vs_Obs_all.png", width = 8, height= 5, units = "in", dpi = 300)
# pooled data only
plot(x = finalData$obs.biomass.live[finalData$trainingData %in% "all"], y = finalData$plot.err.live[finalData$trainingData %in% "all"],
     cex = 0.7, pch = 19, las = 1, ylab = "Estimation error (decimal fraction)", 
     xlab = expression("Observed biomass (g "%.%m^-2~")"))

# compare with magnitude of error?
ggplot(data = finalData, aes(x = obs.biomass.live, y = biomass.error)) + # geom_point() + 
  theme_bw() + labs(x = expression("Observed biomass (g "%.%m^-2~")"), y = expression("Estimation error (g "%.%m^-2~")")) + 
  scale_colour_discrete(name = "Training data") + geom_point() +
  facet_grid(plot ~ .)

ggplot(data = finalData, aes(x = obs.biomass.live, y = biomass.error, colour = as.factor(trainingData))) + # geom_point() + 
  theme_bw() + labs(x = expression("Observed biomass (g "%.%m^-2~")"), 
                    y = expression("Estimation error (g "%.%m^-2~")")) + 
  guides(colour = guide_legend(title = "Training data", override.aes = list(alpha = 1), ncol = 2)) +
  scale_colour_discrete(name = "Training data") + geom_point() 
# some really large errors produced by Jan/Feb's data

ggplot(data = finalData[finalData$trainingData %in% "all", ], aes(x = obs.biomass.live, y = biomass.error)) + # geom_point() + 
  theme_bw() + labs(x = expression("Observed biomass (g "%.%m^-2~")"), y = expression("Estimation error (g "%.%m^-2~")")) + 
  scale_colour_discrete(name = "Training data") + geom_point() +
  facet_grid(plot ~ .)



### More could be done to investigate outliers (e.g., Dec 2013; facet by trainingData) etc. 
# but the general conclusion would remain: there are massive uncertainties introduced by using 
# a single month's data to parameterize allometry equations
### allometry data 

### which month minimizes error?
# will have to re-do this after examining Jan-14 live and Jun-14 dead
# get average biomass.error and plot.err by trainingData
errSum <- ddply(finalData, .(trainingData), summarise,
                 # magnitude of errors (g/m2)
                 live.mag = mean(biomass.error, na.rm = T),
                 dead.mag = mean(dead.error, na.rm = T),
                 live.pct = mean(plot.err.live, na.rm = T),
                 dead.pct = mean(plot.err.dead, na.rm = T)
  )

errSE <- ddply(finalData, .(trainingData), summarise,
                # magnitude of errors (g/m2)
                live.mag.se = se(biomass.error),
                dead.mag.se = se(dead.error),
                live.pct.se = se(plot.err.live),
                dead.pct.se = se(plot.err.dead)
)
allDat <- join_all(list(errSum, errSE))

ord <- c(as.character(sort(as.yearmon(unique(errSum$trainingData)[2:28], "%B-%y"))), unique(errSum$trainingData)[1])

min

# combine for ggplot2
errMelt <- melt(errSum, id.vars = "trainingData")
seMelt  <- melt(errSE,  id.vars = "trainingData")
errMelt$se <- seMelt$value

# percentages
ggplot(errMelt[grep(".pct", errMelt$variable), ], aes(y = value, x = trainingData)) + geom_point() + facet_grid(variable ~ .) +
  theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  labs(x = "Training data", y = "Estimation error (decimal fraction)") + 
  ylim(-0.5, 2) + geom_errorbar(aes(x = trainingData, ymin = value - se , ymax = value + se), width = 0)
# ggsave("C:/RDATA/SPAL_allometry/PctErrorByTrainingData.png", width = 8, height= 5, units = "in", dpi = 300)

# magnitude
ggplot(errMelt[grep(".mag", errMelt$variable), ], aes(y = value, x = trainingData)) + geom_point() + facet_grid(variable ~ .) +
  theme_bw() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  labs(x = "Training data", y = expression("Estimation error (g "%.%m^-2~")")) + 
  ylim(-250, 1500) + geom_errorbar(aes(x = trainingData, ymin = value - se , ymax = value + se), width = 0)
# ggsave("C:/RDATA/SPAL_allometry/massErrorByTrainingData.png", width = 8, height= 5, units = "in", dpi = 300)


### still need to look at seasonal training data

### what is the effect of this error on biomass estimates? 
