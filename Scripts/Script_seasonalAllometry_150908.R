### Script assesses LUMCON allometry equations - seasonality, spatial variation

library(zoo)
library(plyr)
library(knitr)

# load custom functions
# source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version

# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150826.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150827.R")






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
for (i in 1:length(seasons)) {
  for (j in 1:length(year)) {
    # account for winter spanning two years
    if (i == 4) {
      targetDates <- paste0(seasons[[i]], "-", c(year[j], year[j] + 1, year[j] + 1))
    } else {
      targetDates <- paste0(seasons[[i]], "-", year[j])
    }  
    cwc$season[(cwc$monthYear %in% targetDates) & (cwc$marsh %in% "LUM")] <- paste(names(seasons)[i], year[j])
    subData          <- cwc[(cwc$monthYear %in% targetDates) & (cwc$marsh %in% "LUM"), ]
    y.live           <- subData$mass[subData$type %in% "LIVE"]
    x.live           <- subData$hgt[subData$type %in% "LIVE"]
    if (length(y.live) > 5) {
      liveCoefs      <- coef(model <- nls(y.live ~ I(a * x.live^b), start = list(a = 0.1, b = 0.1)))
      predicted           <- liveCoefs[1] * x.live^liveCoefs[2]
      squared_error       <- (predicted - y.live)^2
      liveVals   <- data.frame(coef.live = liveCoefs[1],
                               exp.live  = liveCoefs[2],
                               MSE.live  = sum(squared_error, na.rm = T),
                               r.live    = sqrt(1 - (deviance(model)  / sum((y.live[!is.na(y.live)] - mean(y, na.rm = T))^2)))
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
      deadCoefs      <- coef(model <- nls(y.dead ~ I(a * x.dead^b), start = list(a = 0.1, b = 0.1)))
      predicted           <- deadCoefs[1] * x.dead^deadCoefs[2]
      squared_error       <- (predicted - y.dead)^2
      deadVals   <- data.frame(coef.dead = deadCoefs[1],
                               exp.dead  = deadCoefs[2],
                               MSE.dead  = sum(squared_error, na.rm = T),
                               r.dead    = sqrt(1 - (deviance(model)  / sum((y.dead[!is.na(y.dead)] - mean(y, na.rm = T))^2)))
      )
    } else {
      deadVals   <- data.frame(coef.dead = NA,
                               exp.dead  = NA,
                               MSE.dead  = NA,
                               r.dead    = NA
      )
    }
    
    paramsSite_temp <- cbind(liveVals, deadVals)
    paramsSite_temp$season <- paste(names(seasons)[i], year[j])
    
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
paramsSite[[1]]

# predict biomass using spring 2014 data pooled from all LUM plots
spr14 <- predictBiomass(monthYear = paste0(seasons$sprg, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
sum14  <- predictBiomass(monthYear = paste0(seasons$sumr, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
fall14 <- predictBiomass(monthYear = paste0(seasons$fall, c("-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
wint14 <- predictBiomass(monthYear = paste0(seasons$wint, c("-14", "-15", "-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)

sum13 <- predictBiomass(monthYear = paste0(seasons$sumr, c("-13")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
fall13 <- predictBiomass(monthYear = paste0(seasons$sumr, c("-13")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
wint13 <- predictBiomass(monthYear = paste0(seasons$wint, c("-13", "-14", "-14")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)

spr15 <- predictBiomass(monthYear = paste0(seasons$sumr, c("-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
sum15 <- predictBiomass(monthYear = paste0(seasons$sumr, c("-15")), plot = c("LUM1", "LUM2", "LUM3"), start_nls = 0.05)
head(sum15[[1]])

# estimate NAPP? compare plot biomass estimates?
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


# predictBiomass[[1]] is supposed to be a plot-scale summary, but doesn't work when multiple times or plots are used in training data
