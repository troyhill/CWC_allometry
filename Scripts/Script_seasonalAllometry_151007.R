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
  # spring: Mar Apr May
  sprg = c("Mar", "Apr", "May"),
  # summer: Jun Jul Aug
  sumr = c("Jun", "Jul", "Aug"),
  # fall: Sep Oct Nov
  fall = c("Sep", "Oct", "Nov"),
  # winter: Dec Jan Feb
  wint = c("Dec", "Jan", "Feb")
)

year <- c(13:15)

# Seasonal allometry by "region"
# Get allometry after forming pooled, regional datasets
# this reduces sample size to 1-3; n = 1 per year
startVals <- 0.2
### tdh 151006: why is this producing rows with NA as site?
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
tempData      <- plotParams[[1]][!is.na(plotParams[[1]]$exp.live) & (plotParams[[1]]$marsh %in% "LUM"), c("exp.live", "seas", "marsh")]
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
plot(confint(amod_glht))

# assuming equal variances: sprg-fall; wint-sprg
# unequal variances: sprg-fall; sumr-fall; sumr-sprg






