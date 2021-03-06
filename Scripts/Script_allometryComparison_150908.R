### Script assesses LUMCON allometry equations - seasonality, spatial variation

library(zoo)
library(plyr)
library(knitr)

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150826.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150827.R")


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
    temp <- predictBiomass(monthYear = c(unique(cwc$monthYear)[i]), plot = unique(cwc$site)[j], returnData = "FALSE", cutoff = 7)
  
    if (i != 1) {
      intData <- rbind(intData, temp)
    } else {
      intData <- temp
    }
  }
  temp <- predictBiomass(monthYear = c(unique(cwc$monthYear[cwc$site %in% plot])), plot = unique(cwc$site)[j], returnData = "FALSE", cutoff = 7)
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
  ylim(-100, 200) + geom_errorbar(aes(x = trainingData, ymin = value - se , ymax = value + se), width = 0)
# ggsave("C:/RDATA/SPAL_allometry/massErrorByTrainingData.png", width = 8, height= 5, units = "in", dpi = 300)


kable(allDat, digits = 2)

### still need to look at seasonal training data












### what is the effect of this error on biomass estimates? 
for (i in 1:length(unique(cwc$site))) {
  intData <- predictBiomass(monthYear = unique(cwc$monthYear), 
                          plot = unique(cwc$site)[i], returnData = "TRUE")[[2]]
  if (i != 1) {
    nappEst2 <- rbind(nappEst2, intData)
  } else {
    nappEst2 <- intData
  }
}

nappEst2.plots <- ddply(nappEst2, .(site, time, type, marsh, monthYear), summarise,
                       mass         =  sum(mass, na.rm = T) / plotSize,
                       pred.mass    =  sum(predicted, na.rm = T) / plotSize
)

### make sure that an absence of biomass has not been interpreted as an absence of sampling
nrow(nappEst2.plots) # 207
for (i in 1:length(unique(nappEst2.plots$site))) {
  for (j in 1:length(unique(nappEst2.plots$time))) {
    plot    <- unique(nappEst2.plots$site)[i]
    time    <- unique(nappEst2.plots$time)[j]
    subData <- nappEst2.plots[(nappEst2.plots$site %in% plot) & (nappEst2.plots$time == time), ]
    if (nrow(subData) == 1) { 
      fillData                     <- subData[1, ]
      fillData$type                <- ifelse(fillData$type %in% "LIVE", "DEAD", "LIVE")
      fillData[, 5:ncol(fillData)] <- as.numeric(0)
      nappEst2.plots <- rbind(nappEst2.plots, fillData)
    }
  }
  rownames(nappEst2.plots) <- 1:nrow(nappEst2.plots)
}
nrow(nappEst.plots) # 5 larger (212)
tail(nappEst.plots) 
nappEst2.plots$year <- substr(as.character(nappEst2.plots$time), 5, 9)

nappEst2proc <- ddply(nappEst2.plots, .(marsh, site, time, year), summarise,
                      live      = mass[type %in% "LIVE"],
                      dead      = mass[type %in% "DEAD"],
                      live.pred = pred.mass[type %in% "LIVE"],
                      dead.pred = pred.mass[type %in% "DEAD"]
)


# observed data
napp.obs2 <- nappCalc(nappEst2proc[nappEst2proc$year %in% c(2013, 2014), ], summarize = "TRUE")
napp.obs2$summary <- marshName(napp.obs2$summary)

# predicted data [problem]
napp.pred2 <- nappCalc(nappEst2proc[nappEst2proc$year %in% c(2013, 2014), ], liveCol = "live.pred", 
                      deadCol = "dead.pred", summarize = "TRUE")
napp.pred2$summary <- marshName(napp.pred2$summary)



# another couple ddplys for good luck...
dd.napp.pred2   <- ddply(napp.pred2$summary, .(marsh, year), summarise,
                        smalley.pred  = mean(napp.smalley, na.rm = T),
                        MH.pred       = mean(napp.MH, na.rm = T),
                        VTS.pred      = mean(napp.VTS, na.rm = T),
                        psc.live.pred = mean(napp.psc.a, na.rm = T),
                        psc.tot.pred  = mean(napp.psc.b, na.rm = T)
)
dd.napp.se.pred2   <- ddply(napp.pred2$summary, .(marsh, year), summarise,
                           smalley.se.pred  = se(napp.smalley),
                           MH.se.pred       = se(napp.MH),
                           VTS.se.pred      = se(napp.VTS),
                           psc.live.se.pred = se(napp.psc.a),
                           psc.tot.se.pred  = se(napp.psc.b)
)

m.napp.pred2     <- melt(dd.napp.pred2, id.vars = c("marsh", "year"))
m.napp.se.pred2  <- melt(dd.napp.se.pred2, id.vars = c("marsh", "year"))
m.napp.pred2$value.se  <- m.napp.se.pred2$value

nappAll2 <- rbind(m.napp, m.napp.pred2)
nappAll2$class <- "Observed"
nappAll2$class[grep(".pred", nappAll2$variable)] <- "Predicted"
nappAll2$variable[grep(".pred", nappAll2$variable)] <- gsub(".pred", "", nappAll2$variable[grep(".pred", nappAll2$variable)])


ggplot(nappAll2[(!nappAll2$year %in% "2015"), ], aes(x = class, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(marsh ~ year) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_compare_ObsPred_siteAlloms.png", width = 8, height= 6, units = "in", dpi = 300)











### Use a single allometric equation, built from all plots and all sampling points
# re-build "finalData", but include coefficients
nappEst <- predictBiomass(monthYear = unique(cwc$monthYear), 
               plot = unique(cwc$site), returnData = "TRUE",
               coefReturn = "TRUE")

nappEst[[1]]$monthYear     <- as.yearmon(nappEst[[1]]$monthYear, "%B-%y")
nappEst[[1]]$plot          <- as.character(nappEst[[1]]$plot)
nappEst[[1]]$trainingData  <- as.character(nappEst[[1]]$trainingData)
nappEst[[1]]$trainingData[nchar(nappEst[[1]]$trainingData) > 6] <- "all"
nappEst[[1]]$biomass.error <- nappEst[[1]]$pred.biomass.live - nappEst[[1]]$obs.biomass.live # magnitude of biomass estimation error (g m2 yr)
nappEst[[1]]$dead.error    <- nappEst[[1]]$pred.biomass.dead - nappEst[[1]]$obs.biomass.dead    # magnitude of dead biomass estimation error (g m2 yr)

nappEst[[3]]

# # check manually
# y.live <- cwc$mass[cwc$type %in% "LIVE"]
# x.live <- cwc$hgt[cwc$type %in% "LIVE"]
# y.live2 <- y.live[!is.na(y.live) & !is.na(x.live)]
# x.live2 <- x.live[!is.na(y.live) & !is.na(x.live)]
# # same result
# live.coefs <- coef(model <- nls(y.live2 ~ I(a * exp(b * x.live2)), start = list(a = 0.03, b = 0.03)))
#
# # do the same for dead material
# y.dead <- cwc$mass[cwc$type %in% "DEAD"]
# x.dead <- cwc$hgt[cwc$type %in% "DEAD"]
# y.dead2 <- y.dead[!is.na(y.dead) & !is.na(x.dead)]
# x.dead2 <- x.dead[!is.na(y.dead) & !is.na(x.dead)]
# # same result, again
# dead.coefs <- coef(model <- nls(y.dead2 ~ I(a * exp(b * x.dead2)), start = list(a = 0.03, b = 0.03)))


# build dataset similar to napp
head(nappEst[[2]])

nappEst.plots <- ddply(nappEst[[2]], .(site, time, type, marsh, monthYear), summarise,
                   mass         =  sum(mass, na.rm = T) / plotSize,
                   pred.mass    =  sum(predicted, na.rm = T) / plotSize
)

### make sure that an absence of biomass has not been interpreted as an absence of sampling
nrow(nappEst.plots) # 207
for (i in 1:length(unique(nappEst.plots$site))) {
  for (j in 1:length(unique(nappEst.plots$time))) {
    plot    <- unique(nappEst.plots$site)[i]
    time    <- unique(nappEst.plots$time)[j]
    subData <- nappEst.plots[(nappEst.plots$site %in% plot) & (nappEst.plots$time == time), ]
    if (nrow(subData) == 1) { 
      fillData                     <- subData[1, ]
      fillData$type                <- ifelse(fillData$type %in% "LIVE", "DEAD", "LIVE")
      fillData[, 5:ncol(fillData)] <- as.numeric(0)
      nappEst.plots <- rbind(nappEst.plots, fillData)
    }
  }
  rownames(nappEst.plots) <- 1:nrow(nappEst.plots)
}
nrow(nappEst.plots) # 5 larger (212)
tail(nappEst.plots) # verified these 150825
nappEst.plots$year <- substr(as.character(nappEst.plots$time), 5, 9)

nappEst_proc <- ddply(nappEst.plots, .(marsh, site, time, year), summarise,
              live      = mass[type %in% "LIVE"],
              dead      = mass[type %in% "DEAD"],
              live.pred = pred.mass[type %in% "LIVE"],
              dead.pred = pred.mass[type %in% "DEAD"]
)


# observed data
napp.obs <- nappCalc(nappEst_proc[nappEst_proc$year %in% c(2013, 2014), ], summarize = "TRUE")
napp.obs$summary
napp.obs$summary <- marshName(napp.obs$summary)

# predicted data [problem]
napp.pred <- nappCalc(nappEst_proc[nappEst_proc$year %in% c(2013, 2014), ], liveCol = "live.pred", 
                      deadCol = "dead.pred", summarize = "TRUE")
napp.pred$summary
napp.pred$summary <- marshName(napp.pred$summary)




# another couple ddplys for good luck...
dd.napp.pred   <- ddply(napp.pred$summary, .(marsh, year), summarise,
                   smalley.pred  = mean(napp.smalley, na.rm = T),
                   MH.pred       = mean(napp.MH, na.rm = T),
                   VTS.pred      = mean(napp.VTS, na.rm = T),
                   psc.live.pred = mean(napp.psc.a, na.rm = T),
                   psc.tot.pred  = mean(napp.psc.b, na.rm = T)
)
dd.napp.se.pred   <- ddply(napp.pred$summary, .(marsh, year), summarise,
                      smalley.se.pred  = se(napp.smalley),
                      MH.se.pred       = se(napp.MH),
                      VTS.se.pred      = se(napp.VTS),
                      psc.live.se.pred = se(napp.psc.a),
                      psc.tot.se.pred  = se(napp.psc.b)
)

m.napp.pred     <- melt(dd.napp.pred, id.vars = c("marsh", "year"))
m.napp.se.pred  <- melt(dd.napp.se.pred, id.vars = c("marsh", "year"))
m.napp.pred$value.se  <- m.napp.se.pred$value

nappAll <- rbind(m.napp, m.napp.pred)
# head(nappAll)
nappAll$class <- "Observed"
nappAll$class[grep(".pred", nappAll$variable)] <- "Predicted"
nappAll$variable[grep(".pred", nappAll$variable)] <- gsub(".pred", "", nappAll$variable[grep(".pred", nappAll$variable)])



ggplot(nappAll[(!nappAll$year %in% "2015"), ], aes(x = class, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(marsh ~ year) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_compare_ObsPred.png", width = 8, height= 6, units = "in", dpi = 300)





### Really what I want to do is apply an equation generated by, e.g., LUM1, to the other plots
### and compare the predictions
# LUM1 2014
lum1coefs <- predictBiomass(monthYear = grep("14", unique(cwc$monthYear), value = TRUE), 
                          plot = "LUM1", returnData = "TRUE",
                          coefReturn = "TRUE")[[3]]

cwc.pr <- cwc
cwc.pr$LUM1.pred <- NA
cwc.pr$LUM1.pred[cwc.pr$type %in% "LIVE"] <- lum1coefs$coef.live * exp(lum1coefs$exp.live * cwc.pr$hgt[cwc.pr$type %in% "LIVE"])
cwc.pr$LUM1.pred[cwc.pr$type %in% "DEAD"] <- lum1coefs$coef.dead * exp(lum1coefs$exp.dead * cwc.pr$hgt[cwc.pr$type %in% "DEAD"])

lum1Est <- ddply(cwc.pr, .(site, time, type, marsh, monthYear), summarise,
                       mass         =  sum(mass, na.rm = T) / plotSize,
                       pred.mass    =  sum(LUM1.pred, na.rm = T) / plotSize
)

### make sure that an absence of biomass has not been interpreted as an absence of sampling
nrow(lum1Est) # 207
for (i in 1:length(unique(lum1Est$site))) {
  for (j in 1:length(unique(lum1Est$time))) {
    plot    <- unique(lum1Est$site)[i]
    time    <- unique(lum1Est$time)[j]
    subData <- lum1Est[(lum1Est$site %in% plot) & (lum1Est$time == time), ]
    if (nrow(subData) == 1) { 
      fillData                     <- subData[1, ]
      fillData$type                <- ifelse(fillData$type %in% "LIVE", "DEAD", "LIVE")
      fillData[, 5:ncol(fillData)] <- as.numeric(0)
      lum1Est <- rbind(lum1Est, fillData)
    }
  }
  rownames(lum1Est) <- 1:nrow(lum1Est)
}

lum1Est$year <- substr(as.character(lum1Est$time), 5, 9)

lum1EstProc <- ddply(lum1Est, .(marsh, site, time, year), summarise,
                      live      = mass[type %in% "LIVE"],
                      dead      = mass[type %in% "DEAD"],
                      live.pred = pred.mass[type %in% "LIVE"],
                      dead.pred = pred.mass[type %in% "DEAD"]
)


# observed data
napp.obs <- nappCalc(lum1EstProc[lum1EstProc$year %in% c(2013, 2014), ], summarize = "TRUE")
napp.obs$summary <- marshName(napp.obs$summary)

# predicted data [problem]
napp.pred <- nappCalc(lum1EstProc[lum1EstProc$year %in% c(2013, 2014), ], liveCol = "live.pred", 
                      deadCol = "dead.pred", summarize = "TRUE")
napp.pred$summary <- marshName(napp.pred$summary)


# another couple ddplys for good luck...
dd.napp.pred   <- ddply(napp.pred$summary, .(marsh, year), summarise,
                        smalley.pred  = mean(napp.smalley, na.rm = T),
                        MH.pred       = mean(napp.MH, na.rm = T),
                        VTS.pred      = mean(napp.VTS, na.rm = T),
                        psc.live.pred = mean(napp.psc.a, na.rm = T),
                        psc.tot.pred  = mean(napp.psc.b, na.rm = T)
)
dd.napp.se.pred   <- ddply(napp.pred$summary, .(marsh, year), summarise,
                           smalley.se.pred  = se(napp.smalley),
                           MH.se.pred       = se(napp.MH),
                           VTS.se.pred      = se(napp.VTS),
                           psc.live.se.pred = se(napp.psc.a),
                           psc.tot.se.pred  = se(napp.psc.b)
)

m.napp.pred     <- melt(dd.napp.pred, id.vars = c("marsh", "year"))
m.napp.se.pred  <- melt(dd.napp.se.pred, id.vars = c("marsh", "year"))
m.napp.pred$value.se  <- m.napp.se.pred$value

nappAll <- rbind(m.napp, m.napp.pred)
nappAll$class <- "Observed"
nappAll$class[grep(".pred", nappAll$variable)] <- "Predicted"
nappAll$variable[grep(".pred", nappAll$variable)] <- gsub(".pred", "", nappAll$variable[grep(".pred", nappAll$variable)])



ggplot(nappAll[(!nappAll$year %in% "2015"), ], aes(x = class, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(marsh ~ year) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_ObsPred_LUM1data.png", width = 8, height= 6, units = "in", dpi = 300)
