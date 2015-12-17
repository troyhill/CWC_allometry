### Script explores relationships between LUM NAPP, MSL, climate, bay nutrients
### Using May-September data

source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_151006.R")
source("C:/RDATA/functions/NCDC_LCD_txt.R")
source("C:/RDATA/wll/lumcon/Script_WaterLevelCorrection_151117.R")
months_used    <- c(month.abb[6:8])
months_usedNumeric <- which(month.abb %in% months_used)

# summary NAPP data
napp.sum <- cbind(dd.napp[1:3, ], dd.napp.se[1:3, 3:ncol(dd.napp.se)])

# summarize nutrient data (bay only)
str(nut.prep) # see also nut.se

### nutrient data is too incomplete to be useful. Only 2014-May 2015
nut.sum <- nut.prep # cbind(nut.prep, nut.se)
nut.sum$year <- substr(as.character(nut.sum$moYr), 5, 9)
nut.sum <- ddply(nut.sum[(substr(as.character(nut.sum$moYr), 1, 3) %in% months_used) & (nut.sum$marsh %in% "LUM"), ],
                 .(marsh, year), colwise(mean))


### NCDC temperature/rainfall
ncdc_NOLA$year <- substr(as.character(ncdc_NOLA$moYr), 1, 4)
grSeas <- ddply(ncdc_NOLA[as.numeric(substr(as.character(ncdc_NOLA$moYr), 6, 7)) %in% months_usedNumeric, ], 
                .(year, location), summarise,
                maxTemp  = mean(meanDailyMaxTemp),
                meanTemp = mean(meanDailyMeanTemp),
                meanTemp.se = se(meanDailyMeanTemp),
                meanCDD  = mean(meanDailyCDD),
                ppt      = sum(totalPrecip)
)
grSeas


### MSL, hydroperiod data
msl <- ddply(pf_mnthly[pf_mnthly$Month %in% c(months_usedNumeric), 2:7], .(Year), colwise(mean, na.rm = T))
msl <- msl[msl$Year %in% 2013:2015, ]


comp <- cbind(grSeas, napp.sum, msl[, 4:6])
plot(comp[, c(4, 7, 10, 14, 22:24)])

summary(lm1 <- lm(smalley ~ meanTemp, data = comp))
summary(lm2 <- lm(smalley ~ ppt, data = comp))
summary(lm3 <- lm(smalley ~ MHW, data = comp))

# png(filename = "NAPP_temp_ErrBars.png", width = 5, height = 4, units = "in", res = 300)
par(fig = c(0,1,0,1), mar = c(4, 4.5, 0.5, 0.5))
plot(smalley ~ meanTemp, data = comp, 
     pch = 19, cex = 0.6, las = 1,
     ylim = c(1700, 4300), xlim = c(27.5, 30),
     ylab = expression("NAPP (Smalley: g"%.%m^-2%.%~yr^-1~")"),
     xlab = "mean daily temperature (June-Aug; degrees C)")
segments(comp$meanTemp + comp$meanTemp.se, comp$smalley, comp$meanTemp - comp$meanTemp.se, comp$smalley)
segments(comp$meanTemp, comp$smalley + comp$smalley.se, comp$meanTemp, comp$smalley - comp$smalley.se)
abline(lm1, col = "red")
# dev.off()
