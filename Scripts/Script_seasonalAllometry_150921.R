### Script assesses LUMCON allometry equations - seasonality, spatial variation

library(zoo)
library(plyr)
library(knitr)

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150918.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150918.R")






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

# all sites are combined here
startVals <- 0.3
for (h in 1:length(unique(cwc$marsh))) {
  targSite <- as.character(unique(cwc$site)[h])
  for (i in 1:length(seasons)) {
    for (j in 1:length(year)) {
      # account for winter spanning two years
      if (i == 4) {
        targetDates <- paste0(seasons[[i]], "-", c(year[j], year[j] + 1, year[j] + 1))
      } else {
        targetDates <- paste0(seasons[[i]], "-", year[j])
      }  
      cwc$season[(cwc$monthYear %in% targetDates) & (cwc$site %in% targSite)] <- paste(names(seasons)[i], year[j])
      subData          <- cwc[(cwc$monthYear %in% targetDates) & (cwc$site %in% targSite), ]
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
}
paramsSite[[1]]

# predict biomass using spring 2014 data pooled from all LUM plots
# seasonalData <- list(
#   sum13  = predictBiomass(monthYear = paste0(seasons$sumr, c("-13")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#   fall13 = predictBiomass(monthYear = paste0(seasons$sumr, c("-13")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#   wint13 = predictBiomass(monthYear = paste0(seasons$wint, c("-13", "-14", "-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
# 
#   spr14  = predictBiomass(monthYear = paste0(seasons$sprg, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#   sum14  = predictBiomass(monthYear = paste0(seasons$sumr, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#   fall14 = predictBiomass(monthYear = paste0(seasons$fall, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#   wint14 = predictBiomass(monthYear = paste0(seasons$wint, c("-14", "-15", "-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#     
#   spr15 = predictBiomass(monthYear = paste0(seasons$sumr, c("-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05),
#   sum15 = predictBiomass(monthYear = paste0(seasons$sumr, c("-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
# )


sum13  <- predictBiomass(monthYear = paste0(seasons$sumr, c("-13")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
fall13 <- predictBiomass(monthYear = paste0(seasons$sumr, c("-13")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
wint13 <- predictBiomass(monthYear = paste0(seasons$wint, c("-13", "-14", "-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)

spr14  <- predictBiomass(monthYear = paste0(seasons$sprg, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
sum14  <- predictBiomass(monthYear = paste0(seasons$sumr, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
fall14 <- predictBiomass(monthYear = paste0(seasons$fall, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
wint14 <- predictBiomass(monthYear = paste0(seasons$wint, c("-14", "-15", "-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)

spr15  <- predictBiomass(monthYear = paste0(seasons$sumr, c("-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
sum15  <- predictBiomass(monthYear = paste0(seasons$sumr, c("-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)

# combine seasonal estimates into a single dataframe
sum13[[1]]$trainingSeason  <- "sum13"
fall13[[1]]$trainingSeason <- "fall13"
wint13[[1]]$trainingSeason <- "wint13"

spr14[[1]]$trainingSeason  <- "spr14"
sum14[[1]]$trainingSeason  <- "sum14"
fall14[[1]]$trainingSeason <- "fall14"
wint14[[1]]$trainingSeason <- "wint14"

spr15[[1]]$trainingSeason  <- "spr15"
sum15[[1]]$trainingSeason  <- "sum15"

seasDat <- rbind(sum13[[1]], rbind(fall13[[1]], rbind(wint13[[1]], rbind(spr14[[1]], rbind(sum14[[1]], rbind(fall14[[1]], rbind(wint14[[1]], rbind(spr15[[1]], sum15[[1]]))))))))
seasDat <- marshName(seasDat, siteCol = "plot")


# convert from long to wide before running nappCalc
napp <- ddply(seasDat, .(trainingSeason, marsh, plot, monthYear, year), summarise,
              live.obs = mass.obs[type %in% "LIVE"],
              dead.obs = mass.obs[type %in% "DEAD"],
              live.pred = mass.pred[type %in% "LIVE"],
              dead.pred = mass.pred[type %in% "DEAD"]
)
head(napp)

for(i in 1:length(unique(napp$trainingSeason))) {
  subD <- napp[(napp$trainingSeason %in% unique(napp$trainingSeason)[i]) & (!napp$year %in% 2015), ]
  temp <- nappCalc(subD, liveCol = "live.obs", deadCol = "dead.obs", siteCol = "plot", timeCol = "monthYear", summarize = "TRUE")[[2]]
  temp.pred <- nappCalc(subD, liveCol = "live.pred", deadCol = "dead.pred", siteCol = "plot", timeCol = "monthYear", summarize = "TRUE")[[2]]
  names(temp.pred)[grep("napp", names(temp.pred))] <- paste0("pred.", names(temp.pred)[grep("napp", names(temp.pred))])
  temp <- cbind(temp, temp.pred[, !names(temp.pred) %in% names(temp)])
  
  temp$trainingSeason <- unique(napp$trainingSeason)[i]
  if (!i == 1) {
    nappSeas <- rbind(nappSeas, temp)
  } else {
    nappSeas <- temp
  }
}
head(nappSeas)
nappSeas <- marshName(nappSeas)

dd.napp   <- ddply(nappSeas, .(trainingSeason, year, marsh), numcolwise(mean, na.rm = T))
dd.napp.se   <- ddply(nappSeas, .(trainingSeason, year, marsh), numcolwise(se))

m.napp     <- melt(dd.napp, id.vars = c("marsh", "year", "trainingSeason"))
m.napp.se  <- melt(dd.napp.se, id.vars = c("marsh", "year", "trainingSeason"))
m.napp$value.se  <- m.napp.se$value


m.napp$class <- "Observed"
m.napp$class[grep("pred.", m.napp$variable)] <- "Predicted"
m.napp$variable[grep("pred.", m.napp$variable)] <- gsub("pred.", "", m.napp$variable[grep("pred.", as.character(m.napp$variable))])


ggplot(m.napp[(m.napp$year %in% "2014") & (m.napp$marsh %in% "LUM"), ], aes(x = class, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~ trainingSeason, ncol = 3, nrow = 3) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP14_LUM_byTrainingData.png", width = 8, height= 6, units = "in", dpi = 300)

ggplot(m.napp[(m.napp$year %in% "2013") & (m.napp$marsh %in% "LUM"), ], aes(x = class, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_wrap(~ trainingSeason, ncol = 3, nrow = 3) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP13_LUM_byTrainingData.png", width = 8, height= 6, units = "in", dpi = 300)








### exploratory plots
ggplot(spr14[[2]], aes(y = stem.err, x = mass, col = season)) + geom_point( alpha = 0.65) +
  theme_bw() + facet_grid(type ~ .) + labs(title = "Training data: Spring 2014") +
  ylim(-2, 7) 

ggplot(sum14[[2]], aes(y = stem.err, x = mass, col = season)) + geom_point( alpha = 0.65) +
  theme_bw() + facet_grid(type ~ .) + labs(title = "Training data: Summer 2014") +
  ylim(-2, 7) 

ggplot(fall14[[2]], aes(y = stem.err, x = mass, col = season)) + geom_point( alpha = 0.65) +
  theme_bw() + facet_grid(type ~ .) + labs(title = "Training data: Fall 2014") +
  ylim(-2, 7) 

ggplot(wint14[[2]], aes(y = stem.err, x = mass, col = season)) + geom_point( alpha = 0.65) +
  theme_bw() + facet_grid(type ~ .) + labs(title = "Training data: Winter 2014") +
  ylim(-2, 7) 





