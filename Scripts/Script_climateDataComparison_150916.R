### This script combines and analyzes NCDC and NAPP data 

# to download NOAA tide data
library(VulnToolkit)

# custom functions
simpleCap <- Vectorize(function(x) {
  # capitalize first letter of each word in a string
  # "Vectorize" makes it handle a vector input 
  s <- strsplit(tolower(x), " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
)



# run dependent scripts
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150826.R")
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_NAPPCalc_150827.R")
source("C:/RDATA/functions/NCDC_pdf_fxn_150915.R") # or see https://gist.github.com/troyhill/d730135b099bdd054c98

# load climate data
fileList <- paste0("C:/RDATA/NCDC/NewOrleans/NOLA_annual_", 2013:2014, ".pdf") # doesn't work for pre-2012 pdfs
nola <- NCDC_combine(pdfFiles = fileList)
nola$moYr <- simpleCap(nola$moYr)
names(nola)[12:13] <- substr(names(nola)[12:13], 1, 10)


# load sea level data
psmsl.stations(country = "USA")
gi      <- psmsl(station = 526, interval = "monthly")
gi$mm   <- gi$msl_mm - min(gi$msl_mm, na.rm = T)
gi$yr   <- gi$year
gi$year <- floor(gi$year)
plot(gi$mm ~ gi$yr)
gi.recent <- gi[gi$year > 2012, ]
gi.recent$moYr <- paste0(month.abb, "-", substr(gi.recent$year, 3, 4))

### build summary data to append to napp2$summary
### This isn't very interesting with two datapoints...
ncdcMeans <- ddply(nola[, c(1, 2, 12, 15:17, 19)], .(year), colwise(mean, na.rm = T))
ncdcMax   <- ddply(nola[, c(1:2, 4, 5, 10:13, 15:17, 19)], .(year), colwise(max, na.rm = T))
ncdcMin   <- ddply(nola[, c(1:2, 4:5, 10:13, 15:17, 19)], .(year), colwise(min, na.rm = T))
ncdcSum   <- ddply(nola[, c(10:12, 15:17, 19)], .(year), colwise(sum, na.rm = T))
msl_gi    <- ddply(gi, .(year), summarise, mm = mean(mm, na.rm = T))


ncdcMeans <- join_all(list(napp2$summary, ncdcMeans), by = "year")
ncdcMax   <- join_all(list(napp2$summary, ncdcMax), by = "year")
ncdcMin   <- join_all(list(napp2$summary, ncdcMin), by = "year")
ncdcSum   <- join_all(list(napp2$summary, ncdcSum), by = "year")
msl_gi    <- join_all(list(napp2$summary, msl_gi), by = "year")

# find average growing-season MSL
for (i in 1:nrow(msl_gi)) {
  targYear              <- as.numeric(as.character(msl_gi$year[i]))
  if (targYear == 2015) {
    msl_gi$summerMSL[i] <- NA
  } else {
    msl_gi$summerMSL[i] <- mean(gi$mm[(gi$year == targYear) & (((gi$yr - gi$year) > 0.370) & ((gi$yr - gi$year) < 0.795))], na.rm = T)
  }
}

yVar <- "napp.smalley"
ggplot(ncdcSum, aes(y = eval(parse(text = yVar)), x = PPT_TOTAL_, colour = site)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)
ggplot(ncdcSum, aes(y = eval(parse(text = yVar)), x = DAYS_OVER_0.01, colour = site)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)
ggplot(ncdcMax, aes(y = eval(parse(text = yVar)), x = COOLING_DEGREE_DAYS, colour = site)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)
ggplot(msl_gi, aes(y = eval(parse(text = yVar)), x = summerMSL, colour = site)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)


dd_msl    <- ddply(msl_gi, .(year, marsh), colwise(mean, na.rm = T))
dd_msl_se <- ddply(msl_gi, .(year, marsh), colwise(se))
names(dd_msl_se)[c(4:8)] <- paste0("se.", names(dd_msl_se)[c(4:8)])
dd_msl <- cbind(dd_msl, dd_msl_se[c(4:8)])

# find average growing-season MSL
for (i in 1:nrow(dd_msl)) {
  targYear              <- as.numeric(as.character(dd_msl$year[i]))
  if (targYear == 2015) {
    dd_msl$summerMSL[i] <- NA
  } else {
    dd_msl$summerMSL[i] <- mean(gi$mm[(gi$year == targYear) & (((gi$yr - gi$year) > 0.370) & ((gi$yr - gi$year) < 0.795))], na.rm = T)
  }
}

# pattern holds for LUM and TB-B but not for TB-A
sub <- dd_msl[(dd_msl$marsh %in% "LUM") & (!dd_msl$year %in% 2015), ]
plot(sub$napp.smalley ~ sub$summerMSL, ylim = c(0, 2500), pch = 19, cex = 0.7, 
     ylab = "NAPP", xlab = "growing season MSL (mm)", las = 1)
segments(sub$summerMSL, sub$napp.smalley + sub$se.napp.smalley, sub$summerMSL, sub$napp.smalley - sub$se.napp.smalley)

# biomass and productivity data are in napp2$intervalData and napp2$summary (peak productivity and timing)
### First, explore monthly correlations (should really have a three-month average or sum for the meteorological variables)
plotData      <- napp2$intervalData
plotData$moYr <- paste0(substr(plotData$time, 1, 3), "-", substr(plotData$time, 7, 8))
plotData2     <- join_all(list(plotData, nola, gi.recent[, c("mm", "moYr")]), by = "moYr")

site <- "LUM2" # plotData2[plotData2$site %in% site]
ggplot(plotData2, aes(y = live, x = LOWEST_DAILY_MINIMUM, colour = site)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)
ggplot(plotData2, aes(y = live, x = HIGHEST_DAILY_MAXIMUM, colour = site)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)

ggplot(plotData2, aes(y = live.inc, x = PPT_TOTAL_)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)
ggplot(plotData2, aes(y = live.inc, x = MEAN_DAILY_MAXIMUM)) + geom_point() + theme_bw() + geom_smooth(method="lm", fill=NA)

ggplot(plotData2, aes(y = live, x = mm)) + geom_point() + theme_bw()  + 
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)
ggplot(plotData2, aes(y = live.inc, x = mm)) + geom_point() + theme_bw() + 
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1)

par(mar = c(4, 4, 0.3, 1), fig = c(0,1,0,1))
plot(gi.recent$mm ~ gi.recent$yr, pch = 19, cex = 0.6, type = "l",
     las = 1, 
     ylab = "monthly mean sea level (mm)", xlab = "")
  polygon(x =  c(2013.375, 2013.792, 2013.792, 2013.375), y = c(0, 0, 1e3, 1e3), col = "firebrick2")
  polygon(x =  c(2014.375, 2014.792, 2014.792, 2014.375), y = c(0, 0, 1e3, 1e3), col = "firebrick2")
  lines(gi.recent$mm ~ gi.recent$yr)


### add MSL in previous three months to site.tot
site.tot
# 
# for (i in 1:nrow(site.tot)) {
#   targYear              <- as.numeric(substr(site.tot$time[i], 5, 9))
#   samplingMo            <- substr(site.tot$time[i], 1, 3)
#   
#   month.abb[month.abb %in% samplingMo]
#   sm <- grep(samplingMo, month.abb)
#   tm <- ((sm + 12) - 1:3)
#   targetMos             <- c(paste0(month.abb, "-", targYear - 2001), paste0(month.abb, "-", targYear - 2000))[tm]
#   
#   if (targYear == 2015) {
#     site.tot$MSL_lag[i] <- NA
#   } else {
#     site.tot$MSL_lag[i] <- mean(gi.recent$mm[(gi.recent$year == targYear) & 
#                                 (gi.recent$moYr %in% targetMos)], 
#                                 na.rm = T)
#   }
# }

getMSLlag <- function(data = site.tot, mslData = gi.recent, months = 3,
                      dataTimeCol = "time", mslTimeCol = "moYr") {
  
  for (i in 1:nrow(data)) {
    targYear              <- as.numeric(substr(data[, dataTimeCol][i], 5, 9))
    samplingMo            <- substr(data[, dataTimeCol][i], 1, 3)
    sm <- grep(samplingMo, month.abb)
    tm <- ((sm + 12) - 1:3)
    targetMos             <- c(paste0(month.abb, "-", targYear - 2001), paste0(month.abb, "-", targYear - 2000))[tm]
    
    if (targYear == 2015) {
      data[, "MSL_lag"][i] <- NA
    } else {
      data[, "MSL_lag"][i] <- mean(mslData[, "mm"][(mslData[, "year"] == targYear) & 
          (mslData[, mslTimeCol] %in% targetMos)], na.rm = T)
    }
  }
}

getMSLlag()

plot(site.tot$biomass[site.tot$marsh %in% "LUM"] ~ site.tot$MSL_lag[site.tot$marsh %in% "LUM"], pch = 19, cex = 0.6)
plot(site.tot$biomass[site.tot$marsh %in% "TB-A"] ~ site.tot$MSL_lag[site.tot$marsh %in% "TB-A"], pch = 19, cex = 0.6)
plot(site.tot$biomass[site.tot$marsh %in% "TB-B"] ~ site.tot$MSL_lag[site.tot$marsh %in% "TB-B"], pch = 19, cex = 0.6)
