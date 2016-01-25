### Script assesses LUMCON allometry equations - seasonality, spatial variation
# - Does allometry differ between sites, seasons? 
#     - seasons: plot seasonally-aggregated data together, see if there are differences
#     - sites: use LUM allometry to predict TB biomass, calculate error
# - Compare allometry parameters to literature values
#     - Thursby et al., Morris datasets?

library(zoo)
library(plyr)
library(knitr)
library(multcomp) # for glht
library(sandwich) # for HC3 estimator
library(car) # for Levene's test

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150918.R") # runs CWC_functions.R
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_ancillaryParams_150921.R") # sources no files





#####
### seasonal allometry equations at plot level
#####

### add season column to cwc dataset
cwc$season <- as.character(NA)

seasons <- list(
  sprg = c("Mar", "Apr", "May"),
  sumr = c("Jun", "Jul", "Aug"),
  fall = c("Sep", "Oct", "Nov"),
  wint = c("Dec", "Jan", "Feb")
)

seasons_alt <- list(
  grow = c("Mar", "Apr", "May", "Jun"),
  peak = c("Jul", "Aug", "Sep", "Oct"),
  sene = c("Nov", "Dec", "Jan", "Feb")
)


year <- c(13:15)

# Seasonal allometry by region, season, and year
# Get allometry after forming pooled, regional datasets
# this reduces sample size to 1-3; n = 1 per year
startVals <- 0.2
for (h in 1:length(unique(cwc$marsh))) {
  targSite <- as.character(unique(cwc$marsh)[h])
  for (i in 1:length(seasons)) {
    for (j in 1:length(year)) {
      # account for winter spanning two years
      if (i == 4) {
        targetDates <- paste0(seasons[[i]], "-", c(year[j], year[j] + 1, year[j] + 1))
      } else {
        targetDates <- paste0(seasons[[i]], "-", year[j])
      }  
      cwc$season[(cwc$monthYear %in% targetDates) & (cwc$marsh %in% targSite)] <- paste(names(seasons)[i], year[j])
        
      subData          <- cwc[(cwc$monthYear %in% targetDates) & (cwc$marsh %in% targSite), ]
      y.live           <- subData$mass[subData$type %in% "LIVE"]
      x.live           <- subData$hgt[subData$type %in% "LIVE"]
      if (length(y.live) > 5) {
        liveCoefs      <- coef(model <- nls(y.live ~ I(a * x.live^b), start = list(a = startVals, b = startVals)))
        predicted           <- liveCoefs[1] * x.live^liveCoefs[2]
        squared_error       <- (predicted - y.live)^2
        liveVals   <- data.frame(coef.live = liveCoefs[1],
                                 exp.live  = liveCoefs[2],
                                 MSE.live  = sum(squared_error, na.rm = T),
                                 r.live    = sqrt(1 - (deviance(model)  / sum((y.live[!is.na(y.live)] - mean(y.live, na.rm = T))^2)))
        )
      } else {
        liveVals   <- data.frame(coef.live = NA,
                                 exp.live  = NA,
                                 MSE.live  = NA,
                                 r.live    = NA
        )
      }
      
      y.dead           <- subData$mass[subData$type %in% "DEAD"]
      x.dead           <- subData$hgt[subData$type %in% "DEAD"]
      if (length(y.dead) > 5) {
        deadCoefs      <- coef(model <- nls(y.dead ~ I(a * x.dead^b), start = list(a = startVals, b = startVals)))
        predicted           <- deadCoefs[1] * x.dead^deadCoefs[2]
        squared_error       <- (predicted - y.dead)^2
        deadVals   <- data.frame(coef.dead = deadCoefs[1],
                                 exp.dead  = deadCoefs[2],
                                 MSE.dead  = sum(squared_error, na.rm = T),
                                 r.dead    = sqrt(1 - (deviance(model)  / sum((y.dead[!is.na(y.dead)] - mean(y.dead, na.rm = T))^2)))
        )
      } else {
        deadVals   <- data.frame(coef.dead = NA,
                                 exp.dead  = NA,
                                 MSE.dead  = NA,
                                 r.dead    = NA
        )
      }
      
      paramsSite_temp        <- cbind(liveVals, deadVals)
      paramsSite_temp$season <- paste(names(seasons)[i], year[j])
      paramsSite_temp$site   <- targSite
      
      #     paramsSite_temp <- getAllometryParams(dataset = subData, sitesIncluded = "LUM", combinePlots = "TRUE", returnData = "TRUE")[[2]]
      if ((i != 1) | (j != 1)) {
        paramsSite[[1]] <- rbind(paramsSite[[1]], paramsSite_temp)
        paramsSite[[length(paramsSite) + 1]] <- subData
      } else if ((i == 1) & (j != 1)) {
        paramsSite <- rbind(paramsSite, paramsSite_temp)
        paramsSite[[length(paramsSite) + 1]] <- subData
      } else if ((i == 1) & (j == 1)) {
        paramsSite <- list(parameters = paramsSite_temp)
        paramsSite[[length(paramsSite) + 1]] <- subData
        #       names(paramsSite[[length(paramsSite)]]) <- paste0(names(seasons)[i], year[j])
      }
    }
  }
  if (h != 1) {
    params[[1]] <- rbind(params[[1]], paramsSite[[1]])
#     params[length(params):(length(params) + length(paramsSite) - 1)] <- paramsSite[2:length(paramsSite)]
    rownames(params[[1]]) <- 1:nrow(params[[1]])
  } else if (h == 1) {
    params <- paramsSite
  }
}
# params[[1]]


### Allometry for each season-year-region combination (combining sites)
# this produces the same result as lines 50:127!
for (i in 1:length(unique(cwc$season))) {
  for(j in 1:length(unique(cwc$marsh))) {
    a <- getParams(cwc[(cwc$season %in% unique(cwc$season)[i]) & (cwc$marsh %in% unique(cwc$marsh)[j]), ], 
                   timeName = paste0(unique(cwc$marsh)[j], "-", unique(cwc$season)[i]),
                   dataName = "seasonYearRegion")
    a$marsh <- unique(cwc$marsh)[j]
    a$seas  <- substr(unique(cwc$season)[i], 1, 4)
    if ((i == 1) & (j == 1)) {
      seasonYearRegionParams <- a
    } else {
      seasonYearRegionParams <- rbind(seasonYearRegionParams, a)
    } 
  }
}

ddply(seasonYearRegionParams, .(marsh, seas), summarise,
      coefL = mean(coef.live, na.rm = T),
      expL  = mean(exp.live, na.rm = T),
      coefD = mean(coef.dead, na.rm = T),
      expD  = mean(exp.dead, na.rm = T),
      expL.se = se(exp.live),
      expD.se = se(exp.dead)
      )

### Allometry for each season-region combination (combining sites, years)
for (i in 1:length(unique(substr(cwc$season, 1, 4)))) {
  for(j in 1:length(unique(cwc$marsh))) {
    a <- getParams(cwc[(substr(cwc$season, 1, 4) %in% unique(substr(cwc$season, 1, 4))[i]) & (cwc$marsh %in% unique(cwc$marsh)[j]), ], 
                   timeName = paste0(unique(cwc$marsh)[j], "-", unique(substr(cwc$season, 1, 4))[i]),
                   dataName = "seasonRegion") 
    if ((i == 1) & (j == 1)) {
      seasonRegionParams <- a
    } else {
      seasonRegionParams <- rbind(seasonRegionParams, a)
    } 
  }
}



### Allometry for each season-year combination (combining sites, regions)
for (i in 1:length(unique(cwc$season))) {
  a <- getParams(cwc[(cwc$season %in% unique(cwc$season)[i]), ],
                 timeName = unique(cwc$season)[i],
                 dataName = "seasonYear")
  if (i == 1) {
    seasonYearParams <- a
  } else {
    seasonYearParams <- rbind(seasonYearParams, a)
  }
}


### Allometry for each season (combining sites, regions, and years)
for (i in 1:length(unique(substr(cwc$season, 1, 4)))) {
  season_name <- unique(substr(cwc$season, 1, 4))[i]
  a <- getParams(cwc[(substr(cwc$season, 1, 4) %in% season_name), ], timeName = season_name,
                 returnPlot = TRUE, plotTitle = season_name,
                 dataName = "season")  # could add '& (cwc$marsh %in% "LUM")' if you want to do this for a single region
  if (i == 1) {
    seasonalParams <- a
  } else {
    seasonalParams <- rbind(seasonalParams, a)
  }
}

### Allometry for each year (combining sites, regions, and seasons)
for (i in 1:length(unique(substr(cwc$season, 6, 7)))) {
  targ <- unique(substr(cwc$season, 6, 7))[i]
  a <- getParams(cwc[(substr(cwc$season, 6, 7) %in% targ), ], timeName = targ,
                 returnPlot = TRUE, plotTitle = targ,
                 dataName = "year")  # could add '& (cwc$marsh %in% "LUM")' if you want to do this for a single region
  if (i == 1) {
    yearParams <- a
  } else {
    yearParams <- rbind(yearParams, a)
  }
}


### Allometry for each site (combining seasons and years)
for (i in 1:length(unique(cwc$site))) {
  targ <- unique(cwc$site)[i]
  a <- getParams(cwc[cwc$site %in% targ, ], timeName = targ,
                 returnPlot = TRUE, plotTitle = targ,
                 dataName = "year")  # could add '& (cwc$marsh %in% "LUM")' if you want to do this for a single region
  if (i == 1) {
    siteParams <- a
  } else {
    siteParams <- rbind(siteParams, a)
  }
}


### Allometry for each region (combining sites, seasons, and years)
for (i in 1:length(unique(cwc$marsh))) {
  targ <- unique(cwc$marsh)[i]
  a <- getParams(cwc[cwc$marsh %in% targ, ], timeName = targ,
                 returnPlot = TRUE, plotTitle = targ,
                 dataName = "region")  # could add '& (cwc$marsh %in% "LUM")' if you want to do this for a single region
  if (i == 1) {
    regionParams <- a
  } else {
    regionParams <- rbind(regionParams, a)
  }
}



### Allometry entire dataset (combining sites, regions, seasons, and years)
allParams <- getParams(cwc, timeName = "all",
                 returnPlot = TRUE, plotTitle = "all",
                 dataName = "all")  # could add '& (cwc$marsh %in% "LUM")' if you want to do this for a single region

# # comparison with lm method
# sub <- cwc
# lm1 <- lm(I(log10(sub[, "mass"][sub[, "type"] %in% "LIVE"])) ~ I(log10(sub[, "hgt"][sub[, "type"] %in% "LIVE"])))
# 10^coef(lm1)[1]
# coef(lm1)[2]
# plot(sub[, "mass"][sub[, "type"] %in% "LIVE"] ~ sub[, "hgt"][sub[, "type"] %in% "LIVE"])
#   x <- sub[, "hgt"][sub[, "type"] %in% "LIVE"]
#   y.pred <- 10^coef(lm1)[1] * x^(coef(lm1)[2])
#   y.pred2 <- allParams$coef.live * x^allParams$exp.live
# lines(x = x[order(y.pred)], y = y.pred[order(y.pred)], col = "red", lwd = 2)
# lines(x = x[order(y.pred2)], y = y.pred2[order(y.pred2)], col = "green", lwd = 2)
# nlm1 <- nls(mass ~ a * hgt^b, data = sub[sub$type %in% "LIVE", ],
#             start = list(a = 0.2, b = 0.2))
# AIC(lm1)
# AIC(nlm1)
#anova(lm1, nlm1)
# NLS models perform much better than log transformation on the basis of R2, but not AIC



### summarize region-wide seasonality using exponential terms
params[[1]]$seas <- substr(params[[1]]$season, 1, 4)

params[[1]]      <- params[[1]][!is.na(params[[1]]$site), ]

seaParams <- ddply(params[[1]], .(site, seas), summarise,
                   live.mean   = mean(exp.live, na.rm = T),
                   live.se     = se(exp.live),
                   dead.mean   = mean(exp.dead, na.rm = T),
                   dead.se     = se(exp.dead),
                   r.live.mean = mean(r.live, na.rm = T),
                   r.dead.mean = mean(r.dead, na.rm = T),                 
                   r.live.se   = se(r.live),
                   r.dead.se   = se(r.dead),
                   n           = sum(!is.na(exp.live))
                   )

m.seaParams         <- melt(seaParams[, c(1:3, 5, 7:8)], id.vars = c("site", "seas"))
m.seaParams.se      <- melt(seaParams[, c(1:2, 4, 6, 9:10)], id.vars = c("site", "seas"))
m.seaParams$se.name <- m.seaParams.se$variable
m.seaParams$se      <- m.seaParams.se$value

ggplot(m.seaParams[m.seaParams$variable %in% c("live.mean", "dead.mean"), ], aes(x = seas, y = value, colour = site)) + geom_point() + 
  theme_bw() + theme(legend.title = element_blank()) + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0) + 
  scale_colour_discrete(labels = c("LUMCON", "Bay La Fleur", "Lake Barre")) + 
  labs(x = "Season (split by site)", y = "Allometry exponent") +
  facet_grid(variable ~ site)
# ggsave("C:/RDATA/SPAL_allometry/allomCompare_150921.png", width = 8, height= 6, units = "in", dpi = 300)
# evidence of seasonality (at LUM and TB-A)


ggplot(m.seaParams[m.seaParams$variable %in% c("live.mean", "dead.mean"), ], aes(x = site, y = value)) + geom_point() + 
  theme_bw() + theme(legend.title = element_blank()) + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0) + 
  labs(x = "", y = "Allometry exponent") +
  scale_x_discrete(labels = c("LUMCON", "Bay La
Fleur", "Lake
Barre")) + 
  labs(x = "Site (split by season)", y = "Allometry exponent") +
  facet_grid(seas ~ variable)
#  ggsave("C:/RDATA/SPAL_allometry/allomCompareRegions_150922.png", width = 5, height= 6, units = "in", dpi = 300)
# Differences don't seem strong enough to support regional differences. 




### same analysis as above, but this doesn't pool data for each region-season. 
### plot-level equations are used (with lower MSE and higher predictive power in allometric equations, and 
### larger sample sizes for aggregation)
### seasonal allometry by plot
startVals <- 0.2
for (h in 1:length(levels(cwc$site))) {
  targSite <- as.character(levels(cwc$site)[h]) # c("TB4")[h] 
  for (i in 1:length(seasons)) {
    for (j in 1:length(year)) {
      # account for winter spanning two years
      if (i == 4) {
        targetDates <- paste0(seasons[[i]], "-", c(year[j], year[j] + 1, year[j] + 1))
      } else {
        targetDates <- paste0(seasons[[i]], "-", year[j])
      }  
      #       cwc$season[(cwc$monthYear %in% targetDates) & (cwc$site %in% targSite)] <- paste(names(seasons)[i], year[j])
      subData          <- cwc[(cwc$monthYear %in% targetDates) & (cwc$site %in% targSite), ]
      y.live           <- subData$mass[subData$type %in% "LIVE"]
      x.live           <- subData$hgt[subData$type %in% "LIVE"]
      if (length(y.live) > 5) {
        liveCoefs      <- tryCatch(coef(model <- nls(y.live ~ I(a * x.live^b), start = list(a = startVals, b = startVals))),
                                   error = function(e) c(NA, NA))
        predicted           <- liveCoefs[1] * x.live^liveCoefs[2]
        squared_error       <- (predicted - y.live)^2
        liveVals   <- data.frame(coef.live = liveCoefs[1],
                                 exp.live  = liveCoefs[2],
                                 MSE.live  = ifelse(!is.na(liveCoefs[1]), sum(squared_error, na.rm = T), NA),
                                 r.live    = ifelse(!is.na(liveCoefs[1]), 
                                                    sqrt(1 - (deviance(model)  / sum((y.live[!is.na(y.live)] - mean(y.live, na.rm = T))^2))), 
                                                    NA)
        )
      } else {
        liveVals   <- data.frame(coef.live = NA,
                                 exp.live  = NA,
                                 MSE.live  = NA,
                                 r.live    = NA
        )
      }
      
      y.dead           <- subData$mass[subData$type %in% "DEAD"]
      x.dead           <- subData$hgt[subData$type %in% "DEAD"]
      if (length(y.dead) > 5) {
        deadCoefs      <- tryCatch(coef(model <- nls(y.dead ~ I(a * x.dead^b), start = list(a = startVals, b = startVals))),
                                   error = function(e) c(NA, NA))
        predicted           <- deadCoefs[1] * x.dead^deadCoefs[2]
        squared_error       <- (predicted - y.dead)^2
        deadVals   <- data.frame(coef.dead = deadCoefs[1],
                                 exp.dead  = deadCoefs[2],
                                 MSE.dead  = ifelse(!is.na(deadCoefs[1]), sum(squared_error, na.rm = T), NA),
                                 r.dead    = ifelse(!is.na(deadCoefs[1]), 
                                                    sqrt(1 - (deviance(model)  / sum((y.dead[!is.na(y.dead)] - mean(y.dead, na.rm = T))^2))), 
                                                    NA)
        )
      } else {
        deadVals   <- data.frame(coef.dead = NA,
                                 exp.dead  = NA,
                                 MSE.dead  = NA,
                                 r.dead    = NA
        )
      }
      
      paramsSite_temp        <- cbind(liveVals, deadVals)
      paramsSite_temp$season <- paste(names(seasons)[i], year[j])
      paramsSite_temp$site   <- targSite
      
      #     paramsSite_temp <- getAllometryParams(dataset = subData, sitesIncluded = "LUM", combinePlots = "TRUE", returnData = "TRUE")[[2]]
      if ((i != 1) | (j != 1)) {
        paramsSite[[1]] <- rbind(paramsSite[[1]], paramsSite_temp)
        paramsSite[[length(paramsSite) + 1]] <- subData
      } else if ((i == 1) & (j != 1)) {
        paramsSite <- rbind(paramsSite, paramsSite_temp)
        paramsSite[[length(paramsSite) + 1]] <- subData
      } else if ((i == 1) & (j == 1)) {
        paramsSite <- list(parameters = paramsSite_temp)
        paramsSite[[length(paramsSite) + 1]] <- subData
        #       names(paramsSite[[length(paramsSite)]]) <- paste0(names(seasons)[i], year[j])
      }
    }
  }
  if (h != 1) {
    plotParams[[1]] <- rbind(plotParams[[1]], paramsSite[[1]])
    rownames(plotParams[[1]]) <- 1:nrow(plotParams[[1]])
#     plotParams[length(plotParams):(length(plotParams) + length(paramsSite) - 1)] <- paramsSite[2:length(paramsSite)]
  } else if (h == 1) {
    plotParams <- paramsSite
  }
}
sum(!is.na(plotParams[[1]]$exp.live)) # adjust startVals and maximize this number # max 59
# rm(plotParams)

plotParams[[1]]$seas <- substr(plotParams[[1]]$season, 1, 4)
plotParams[[1]]      <- marshName(plotParams[[1]])
seaPlotParams <- ddply(plotParams[[1]], .(marsh, seas), summarise,
                   live.mean   = mean(exp.live, na.rm = T),
                   live.se     = se(exp.live),
                   dead.mean   = mean(exp.dead, na.rm = T),
                   dead.se     = se(exp.dead),
                   r.live.mean = mean(r.live, na.rm = T),
                   r.dead.mean = mean(r.dead, na.rm = T),                 
                   r.live.se   = se(r.live),
                   r.dead.se   = se(r.dead),
                   n.live      = sum(!is.na(exp.live))
)
names(seaPlotParams)

m.seaPlots         <- melt(seaPlotParams[, c(1:3, 5, 7:8)], id.vars = c("marsh", "seas"))
m.seaPlots.se      <- melt(seaPlotParams[, c(1:2, 4, 6, 9:10)], id.vars = c("marsh", "seas"))
m.seaPlots$se.name <- m.seaPlots.se$variable
m.seaPlots$se      <- m.seaPlots.se$value





ggplot(m.seaPlots[m.seaPlots$variable %in% c("live.mean", "dead.mean"), ], aes(x = seas, y = value, colour = marsh)) + geom_point() + 
  theme_bw() + theme(legend.title = element_blank()) + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0) + 
  scale_colour_discrete(labels = c("LUMCON", "Bay La Fleur", "Lake Barre")) + 
  labs(x = "Season (split by region)", y = "Allometry exponent") +
  facet_grid(variable ~ marsh)
# ggsave("C:/RDATA/SPAL_allometry/allomComparePlots_150922.png", width = 8, height= 6, units = "in", dpi = 300)
# evidence of seasonality (in all regions)


ggplot(m.seaPlots[m.seaPlots$variable %in% c("live.mean", "dead.mean"), ], aes(x = marsh, y = value)) + geom_point() + 
  theme_bw() + theme(legend.title = element_blank()) + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0) + 
  labs(x = "", y = "Allometry exponent") +
  scale_x_discrete(labels = c("LUMCON", "Bay La
Fleur", "Lake
Barre")) + 
  labs(x = "Region (split by season)", y = "Allometry exponent") +
  facet_grid(seas ~ variable)
# ggsave("C:/RDATA/SPAL_allometry/allomComparePlots2_150922.png", width = 5, height= 6, units = "in", dpi = 300)



# ANOVAs with unbalanced design and unequal variances. Following Herberich et al. 2010
tempData      <- plotParams[[1]][!is.na(plotParams[[1]]$exp.live) & (plotParams[[1]]$marsh %in% "LUM"), c("exp.live", "seas", "site", "marsh")]
tempData[, 2] <- ordered(as.character(tempData$seas), c("sprg", "sumr", "fall", "wint"))
leveneTest(exp.live ~ seas, data = tempData) # p = 0.12
# but,  with small sample sizes these tests have low power to detect violations. 
# boxplots are probably the better means of evaluating variance, or just not to 
# assume homogenous variances, since I can't justify that assumption by reasoning about the data
ddply(tempData, .(seas), summarise, exp.live.se = se(exp.live)) # sd varies by a factor of 2.7; se varies by factor of 3.8. I think this justifies 
n <- table(tempData$seas) # observations per season
plot(exp.live ~ seas, data = tempData,
     xlab = "season", ylab = "allometry exponent (live biomass)")
axis(3, at = 1:4, labels = paste("n = ", n))

tempData$seas <- 
  ordered(as.character(tempData$seas), c("sprg", "sumr", "fall", "wint"))

amod          <- aov(exp.live ~ seas, data = tempData)
# glht sets up multiple contrasts. vcov = vcovHC specifies heteroscedasticity-consisten estimator of covariance matrix (Herberich et al. 2010)
amod_glht     <- glht(amod, mcp(seas = "Tukey"), vcov = vcovHC) # remove ", vcov = vcovHC)" if variances are assumed equal 
coef(amod_glht)
summary(amod_glht)
plot(confint(amod_glht)) # sprg-fall; sumr-fall; sumr-sprg

explabs <- c("A", "B", "C", "ABC")
plot(exp.live ~ seas, data = tempData,
     xlab = "season", ylab = "allometry exponent (live biomass)",
     ylim = c(1, 3.1), staplewex = 0)
# axis(3, at = 1:4, labels = paste0("n = ", n))
text(x = 1:4, y = 3.05, explabs)


# png(filename = "allom_exponent_live2.png", width = 3.6, height = 3, units = "in", res = 400)
par(mar = c(3.5, 4, 0.5, 0.7), fig = c(0,1,0,1))
a <- ddply(tempData, .(seas), summarise, mean = mean(exp.live), se = se(exp.live))
plot(a$mean, type = "p", pch = 19, cex = 0.6,
     xlab = "", ylab = "allometry exponent (live stems)",
     ylim = c(0.5, 3),  yaxt = "n", xaxt = "n", bty = "n", yaxs = "i")
axis(1, at = 0:length(levels(a$seas)), labels = c("", "spring", "summer", "fall", "winter") )
axis(2, las = 1, at = axTicks(2), labels = axTicks(2))
# add error bars
segments(index(a$seas), a$mean + a$se, index(a$seas), a$mean - a$se)
text(x = c(1:3, 3.85), y = 2.85, explabs)
# dev.off()



plot(a$mean, type = "n", pch = 19, cex = 0.6,
     xlab = "", ylab = "allometry exponent (live stems)",
     ylim = c(0.5, 3),  yaxt = "n", xaxt = "n", bty = "n", yaxs = "i")
points(exp.live ~ seas, data = tempData, subset = site %in% "LUM1", cex = 0.6, pch = 19, col = "red")
points(exp.live ~ seas, data = tempData, subset = site %in% "LUM2", cex = 0.6, pch = 19)
points(exp.live ~ seas, data = tempData, subset = site %in% "LUM3", cex = 0.6, pch = 19, col = "cornflowerblue")
text(x = 1.5, y = 1.2, "LUM1", col = "red")
text(x = 1.5, y = 1.05, "LUM2")
text(x = 1.5, y = 0.9, "LUM3", col = "cornflowerblue")
axis(1, at = 0:length(levels(tempData$seas)), labels = c("", "spring", "summer", "fall", "winter") )
axis(2, las = 1, at = axTicks(2), labels = axTicks(2))




tempData      <- plotParams[[1]][!is.na(plotParams[[1]]$exp.dead) & (plotParams[[1]]$marsh %in% "LUM"), c("exp.dead", "seas", "site", "marsh")]
tempData$seas <- ordered(as.character(tempData$seas), c("sprg", "sumr", "fall", "wint"))
summary(glht(aov(exp.dead ~ seas, data = tempData), mcp(seas = "Tukey"), vcov = vcovHC)) # remove ", vcov = vcovHC)" if variances are assumed equal 

# png(filename = "allom_exponent_dead2.png", width = 3.8, height = 3, units = "in", res = 400)
par(mar = c(3.5, 4, 0.5, 0.7), fig = c(0,1,0,1))
a <- ddply(tempData, .(seas), summarise, mean = mean(exp.dead), se = se(exp.dead))
plot(a$mean, type = "p", pch = 19, cex = 0.6,
     xlab = "", ylab = "allometry exponent (dead stems)",
     ylim = c(0.5, 3),  yaxt = "n", xaxt = "n", bty = "n", yaxs = "i")
axis(1, at = 0:length(levels(a$seas)), labels = c("", "spring", "summer", "fall", "winter") )
axis(2, las = 1, at = axTicks(2), labels = axTicks(2))

# add error bars
segments(index(a$seas), a$mean + a$se, index(a$seas), a$mean - a$se)
text(x = c(1:4), y = 2.85, c("AB", "A", "B", "B"))
# dev.off()

plot(a$mean, type = "n", pch = 19, cex = 0.6,
     xlab = "", ylab = "allometry exponent (dead stems)",
     ylim = c(0.5, 3),  yaxt = "n", xaxt = "n", bty = "n", yaxs = "i")
points(exp.dead ~ seas, data = tempData, subset = site %in% "LUM1", cex = 0.6, pch = 19, col = "red")
points(exp.dead ~ seas, data = tempData, subset = site %in% "LUM2", cex = 0.6, pch = 19)
points(exp.dead ~ seas, data = tempData, subset = site %in% "LUM3", cex = 0.6, pch = 19, col = "cornflowerblue")
text(x = 1.5, y = 2.8, "LUM1", col = "red")
text(x = 1.5, y = 2.65, "LUM2")
text(x = 1.5, y = 2.5, "LUM3", col = "cornflowerblue")
axis(1, at = 0:length(levels(a$seas)), labels = c("", "spring", "summer", "fall", "winter") )
axis(2, las = 1, at = axTicks(2), labels = axTicks(2))




# Same analysis as above, but for TB-A/B
tempData      <- plotParams[[1]][!is.na(plotParams[[1]]$exp.live) & (plotParams[[1]]$marsh %in% "TB-B"), c("exp.live", "seas", "site", "marsh")]
tempData[, 2] <- as.factor(tempData[, 2])
tempData[, 3] <- as.factor(tempData[, 3])
leveneTest(exp.live ~ seas, data = tempData) # p = 0.12
# but,  with small sample sizes these tests have low power to detect violations. 
# boxplots are probably the better means of evaluating variance, or just not to 
# assume homogenous variances, since I can't justify that assumption by reasoning about the data
ddply(tempData, .(seas), summarise, exp.live.se = se(exp.live)) # sd varies by a factor of 2.7; se varies by factor of 3.8. I think this justifies 
n <- table(tempData$seas) # observations per season
plot(exp.live ~ seas, data = tempData,
     xlab = "season", ylab = "allometry exponent (live biomass)")
axis(3, at = 1:4, labels = paste("n = ", n))


amod          <- aov(exp.live ~ seas, data = tempData)
# glht sets up multiple contrasts. vcov = vcovHC specifies heteroscedasticity-consisten estimator of covariance matrix (Herberich et al. 2010)
amod_glht     <- glht(amod, mcp(seas = "Tukey"), vcov = vcovHC) # remove ", vcov = vcovHC)" if variances are assumed equal 
coef(amod_glht)
summary(amod_glht)
plot(confint(amod_glht)) # sprg-fall; sumr-fall; sumr-sprg

explabs <- c("A", "B", "C", "ABC") # for TB-A
explabs <- c("AB", "AB", "A", "B") # for TB-B

plot(exp.live ~ seas, data = tempData,
     xlab = "season", ylab = "allometry exponent (live biomass)",
     ylim = c(1, 3.5), staplewex = 0)
axis(3, at = 1:4, labels = paste0("n = ", n))
text(x = 1:4, y = 3.45, explabs)




### plot of all data, showing seasonal allometry
### seasonRegionParams, regionParams
type <- "LIVE"
sub <- cwc
# png(filename = "allom_plot_part1.png", width = 10, height = 6, units = "in", res = 300)
par(fig = c(0,1,0,1), mar = c(4, 4, 0.5, 0.5))
plot(sub[, "mass"][sub[, "type"] %in% type] ~ sub[, "hgt"][sub[, "type"] %in% type], cex = 0.35, pch = 19, col = "gray",
     ylab = "mass (g)", xlab = "height (cm)", xlim = c(0, 150), ylim = c(0, 13), las = 1, tcl = 0.25, tck = 0.01, bty = "n", yaxs = "i", xaxs = "i")
  abline(h = 0)
  abline(v = 0)
  x <- sub[, "hgt"][sub[, "type"] %in% type]
  y.pred <- regionParams$coef.live[1] * x^regionParams$exp.live[1] # all seasons, all data
  y.pred2 <- seasonRegionParams$coef.live[seasonRegionParams$time.period %in% "LUM-sprg"] * x^seasonRegionParams$exp.live[seasonRegionParams$time.period %in% "LUM-sprg"] # sprg
  y.pred3 <- seasonRegionParams$coef.live[seasonRegionParams$time.period %in% "LUM-sumr"] * x^seasonRegionParams$exp.live[seasonRegionParams$time.period %in% "LUM-sumr"] # sumr
  y.pred4 <- seasonRegionParams$coef.live[seasonRegionParams$time.period %in% "LUM-fall"] * x^seasonRegionParams$exp.live[seasonRegionParams$time.period %in% "LUM-fall"] # fall
  y.pred5 <- seasonRegionParams$coef.live[seasonRegionParams$time.period %in% "LUM-wint"] * x^seasonRegionParams$exp.live[seasonRegionParams$time.period %in% "LUM-wint"] # wint
lines(x = x[order(y.pred)], y = y.pred[order(y.pred)], lwd = 2, lty = 2)
lines(x = x[order(y.pred2)], y = y.pred2[order(y.pred2)], col = "red1", lwd = 2)
lines(x = x[order(y.pred3)], y = y.pred3[order(y.pred3)], col = "orange", lwd = 2)
lines(x = x[order(y.pred4)], y = y.pred4[order(y.pred4)], col = "green", lwd = 2)
lines(x = x[order(y.pred5)], y = y.pred5[order(y.pred5)], col = "brown4", lwd = 2)
# dev.off()


type <- "DEAD"
# png(filename = "allom_plot_part4.png", width = 10, height = 6, units = "in", res = 300)
par(fig = c(0,1,0,1), mar = c(4, 4, 0.5, 0.5))
plot(sub[, "mass"][sub[, "type"] %in% type] ~ sub[, "hgt"][sub[, "type"] %in% type], cex = 0.35, pch = 19, col = "gray",
     ylab = "mass (g)", xlab = "height (cm)", xlim = c(0, 150), ylim = c(0, 13), las = 1, tcl = 0.25, tck = 0.01, bty = "n", yaxs = "i", xaxs = "i")
abline(h = 0)
abline(v = 0)
x <- sub[, "hgt"][sub[, "type"] %in% type]
y.pred <- regionParams$coef.dead[1] * x^regionParams$exp.dead[1] # all seasons, all data
y.pred2 <- seasonRegionParams$coef.dead[seasonRegionParams$time.period %in% "LUM-sprg"] * x^seasonRegionParams$exp.dead[seasonRegionParams$time.period %in% "LUM-sprg"] # sprg
y.pred3 <- seasonRegionParams$coef.dead[seasonRegionParams$time.period %in% "LUM-sumr"] * x^seasonRegionParams$exp.dead[seasonRegionParams$time.period %in% "LUM-sumr"] # sumr
y.pred4 <- seasonRegionParams$coef.dead[seasonRegionParams$time.period %in% "LUM-fall"] * x^seasonRegionParams$exp.dead[seasonRegionParams$time.period %in% "LUM-fall"] # fall
y.pred5 <- seasonRegionParams$coef.dead[seasonRegionParams$time.period %in% "LUM-wint"] * x^seasonRegionParams$exp.dead[seasonRegionParams$time.period %in% "LUM-wint"] # wint
lines(x = x[order(y.pred)], y = y.pred[order(y.pred)], lwd = 2, lty = 2)
lines(x = x[order(y.pred2)], y = y.pred2[order(y.pred2)], col = "red1", lwd = 2)
lines(x = x[order(y.pred3)], y = y.pred3[order(y.pred3)], col = "orange", lwd = 2)
lines(x = x[order(y.pred4)], y = y.pred4[order(y.pred4)], col = "green", lwd = 2)
lines(x = x[order(y.pred5)], y = y.pred5[order(y.pred5)], col = "brown4", lwd = 2)

# dev.off()
